import abc
import typing


class Summarizable(metaclass=abc.ABCMeta):
    """
    A mixin for entities that can summarize themselves into a provided IO handle.
    """

    @abc.abstractmethod
    def summarize(
        self,
        other: typing.Optional[str] = None,
    ):
        """
        Summarize the item while also considering `other` (default `None`).
        """
        pass
