import abc
import typing

T = typing.TypeVar('T')

class Formatter(typing.Generic[T], metaclass=abc.ABCMeta):
    
    @abc.abstractmethod
    def format_as_string(self, item: T, tx_id: typing.Optional[str]) -> str:
        """ Inputs an item and outputs a human readable string that can be used
        to more easily read the item in tables or other visualizers.

        Args:
            item (T): a genophenocorr class
            tx_id (typing.Optional[str]): The desired transcript ID, if not the MANE transcript

        Returns:
            str: a human readable string 
        """
        pass
    
