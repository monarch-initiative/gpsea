import abc
import io
import json
import pickle
import os
import typing


from gpsea.model import (
    ProteinMetadata,
    TranscriptAnnotation,
    TranscriptCoordinates,
    TranscriptInfoAware,
    VariantCoordinates,
)
from gpsea.io import GpseaJSONDecoder, GpseaJSONEncoder

from ._api import (
    FunctionalAnnotator,
    ProteinMetadataService,
    TranscriptCoordinateService,
)

T = typing.TypeVar("T")


class Cache(typing.Generic[T], metaclass=abc.ABCMeta):
    """
    Cache stores items locally to prevent having to consult a remote resource.

    The itemes are persisted across runs.
    """

    @abc.abstractmethod
    def load_item(self, identifier: str) -> typing.Optional[T]:
        pass

    @abc.abstractmethod
    def store_item(self, identifier: str, item: T):
        pass


class FilesystemCache(typing.Generic[T], Cache[T], metaclass=abc.ABCMeta):

    def __init__(
        self,
        data_dir: str,
    ):
        if not os.path.isdir(data_dir):
            raise ValueError(f"datadir {data_dir} must be an existing directory")
        self._datadir = data_dir

    @property
    def data_dir(self) -> str:
        return self._datadir

    def load_item(self, identifier: str) -> typing.Optional[T]:
        name = self._prepare_resource_path(identifier=identifier)
        if os.path.isfile(name):
            with open(name, "rb") as fh:
                return self._deserialize(fh)  # type: ignore
        else:
            return None

    def store_item(self, identifier: str, item: T):
        name = self._prepare_resource_path(identifier=identifier)
        with open(name, "wb") as fh:
            self._serialize(item, fh)  # type: ignore

    @abc.abstractmethod
    def _prepare_resource_path(self, identifier: str) -> str:
        pass

    @abc.abstractmethod
    def _serialize(self, item: T, fh: io.BytesIO):
        pass

    @abc.abstractmethod
    def _deserialize(self, fh: io.BytesIO) -> T:
        pass


class PicklingCache(typing.Generic[T], FilesystemCache[T]):

    def _prepare_resource_path(self, identifier: str) -> str:
        return os.path.join(self._datadir, f"{identifier}.pickle")

    def _serialize(self, item: T, fh: io.BytesIO):
        pickle.dump(item, fh)

    def _deserialize(self, fh: io.BytesIO) -> T:
        return pickle.load(fh)


class JsonCache(typing.Generic[T], FilesystemCache[T]):

    def __init__(
        self,
        data_dir: str,
        indent: typing.Optional[typing.Union[int, str]] = None,
    ):
        super().__init__(data_dir=data_dir)
        self._indent = indent

    def _prepare_resource_path(self, identifier: str) -> str:
        return os.path.join(self._datadir, f"{identifier}.json")

    def _serialize(
        self,
        item: T,
        fh: io.BytesIO,
    ):
        text_fh = JsonCache._wrap_with_text(fh)
        json.dump(
            item,
            text_fh,
            cls=GpseaJSONEncoder,
            indent=self._indent,
        )

    def _deserialize(
        self,
        fh: io.BytesIO,
    ) -> T:
        text_fh = JsonCache._wrap_with_text(fh)
        return json.load(fp=text_fh, cls=GpseaJSONDecoder)

    @staticmethod
    def _wrap_with_text(fh: io.BytesIO) -> typing.TextIO:
        return io.TextIOWrapper(fh)


class CachingProteinMetadataService(ProteinMetadataService):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        cache: Cache[ProteinMetadata],
        fallback: ProteinMetadataService,
    ):
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, ProteinMetadataService)
        self._fallback = fallback

    def annotate(self, protein_id: str) -> ProteinMetadata:
        assert isinstance(protein_id, str)
        item = self._cache.load_item(protein_id)
        if item is None:  # cache miss
            item = self._fallback.annotate(protein_id)
            self._cache.store_item(protein_id, item)

        return item


class CachingTranscriptCoordinateService(TranscriptCoordinateService):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        cache: Cache[TranscriptCoordinates],
        fallback: TranscriptCoordinateService,
    ):
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, TranscriptCoordinateService)
        self._fallback = fallback

    def fetch(
        self,
        tx: typing.Union[str, TranscriptInfoAware],
    ) -> TranscriptCoordinates:
        tx_id = self._parse_tx(tx)
        item = self._cache.load_item(tx_id)
        if item is None:  # cache miss
            item = self._fallback.fetch(tx_id)
            self._cache.store_item(tx_id, item)

        return item


class CachingFunctionalAnnotator(FunctionalAnnotator):
    """A class that retrieves a Variant object if it exists or will run the fallback Fuctional Annotator if it does not exist.

    Methods:
        annotate(variant_coordinates:VariantCoordinates): Gets data and returns a Variant object for given variant coordinates
    """

    def __init__(
        self,
        cache: Cache[typing.Sequence[TranscriptAnnotation]],
        fallback: FunctionalAnnotator,
    ):
        assert isinstance(cache, Cache)
        self._cache = cache

        assert isinstance(fallback, FunctionalAnnotator)
        self._fallback = fallback

    @staticmethod
    def _create_cache_key(vc: VariantCoordinates) -> str:
        vk = vc.variant_key
        if len(vk) <= 50:
            return f"{vk}"
        else:
            # long INDELs in sequence notation
            return f"{vc.chrom}_{vc.start}_{vc.end}_{vc.variant_class}"

    def annotate(
        self,
        variant_coordinates: VariantCoordinates,
    ) -> typing.Sequence[TranscriptAnnotation]:
        assert isinstance(variant_coordinates, VariantCoordinates)

        cache_key = CachingFunctionalAnnotator._create_cache_key(variant_coordinates)
        annotations = self._cache.load_item(cache_key)
        if annotations is None:  # cache miss
            annotations = self._fallback.annotate(variant_coordinates)
            self._cache.store_item(cache_key, annotations)

        return annotations
