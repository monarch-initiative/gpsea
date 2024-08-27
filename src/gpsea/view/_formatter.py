import abc
import typing

from gpsea.model import Variant

T = typing.TypeVar('T')

class Formatter(typing.Generic[T], metaclass=abc.ABCMeta):
    
    @abc.abstractmethod
    def format_as_string(self, item: T) -> str:
        """ Inputs an item and outputs a human readable string that can be used
        to more easily read the item in tables or other visualizers.

        Args:
            item (T): an element to be formatted

        Returns:
            str: a human readable string 
        """
        pass
    
class VariantFormatter(Formatter[Variant]):
    """
    A class that can be used to format a `Variant` to a human readable string
    """
    def __init__(self, tx_id: typing.Optional[str] = None) -> None:
        self._tx_id = tx_id
        
    def format_as_string(self, item: Variant) -> str:
        """ 
        Args:
            item (Variant): An object of class `Variant` representing a variant.

        Returns:
            str: A human readable string for the variant.
        """
        if self._tx_id is not None:
            transcript = item.get_tx_anno_by_tx_id(self._tx_id)
        else:
            transcript = item.get_preferred_tx_annotation()
        if transcript is not None and transcript.hgvs_cdna is not None:
            if len(transcript.hgvs_cdna) > 50:
                return "Long HGVS"
            return transcript.hgvs_cdna
        else:
            return item.variant_info.variant_key
