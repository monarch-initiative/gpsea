import abc
import io
import typing

from gpsea.util import open_text_io_handle_for_writing


class GpseaReport(metaclass=abc.ABCMeta):
    """
    `GpseaReport` summarizes an aspect of an analysis for the user.
    """

    @abc.abstractmethod
    def write(self, fh: typing.Union[io.IOBase, str]):
        """
        Write the report into the provided path or file handle.

        :param fh: a `str` with path
          or :class:`io.IOBase` with the file-like object for writing the report into.
        """
        pass
    

class HtmlGpseaReport(GpseaReport):
    """
    A report where the content is formatted as a HTML `str`.
    """
    # NOT PART OF THE PUBLIC API
    
    def __init__(
        self,
        html: str,
    ):
        assert isinstance(html, str)
        self._html = html
    
    @property
    def html(self) -> str:
        """
        Get a `str` with the HTML report.
        """
        return self._html

    def write(self, fh: typing.Union[io.IOBase, str]):
        should_close = isinstance(fh, str)
        try:
            fout = open_text_io_handle_for_writing(fh)
            fout.write(self._html)
        except Exception:
            if should_close:
                fout.close()

    def _repr_html_(self) -> str:
        return self._html
