import abc
import io
import typing


class Summarizable(metaclass=abc.ABCMeta):
    """
    A mixin for entities that can summarize themselves into a provided IO handle.
    """

    @abc.abstractmethod
    def summarize(
        self,
        out: typing.TextIO,
    ):
        """
        Summarize the item into the provided IO handle.

        :param out: an IO handle to write into.
        """
        pass

    def summary(self) -> str:
        """
        Get the summary.
        """
        buf = io.StringIO()
        self.summarize(buf)
        return buf.getvalue()
