import abc
import typing

from ._variant_data import Variant, VariantCoordinates

T = typing.TypeVar('T')


class VariantCoordinateFinder(typing.Generic[T], metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def find_coordinates(self, item: T) -> VariantCoordinates:
        """
        Determine :class:`VariantCoordinates` from an `item` of some sort.
        """
        pass


class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        pass
