import abc
import io
import typing

from gpsea.util import open_text_io_handle_for_writing


from jinja2 import Environment, PackageLoader


def pluralize(count: int, stem: str) -> str:
    if isinstance(count, int):
        if count == 1:
            return f"{count} {stem}"
        else:
            if stem.endswith("s"):
                return f"{count} {stem}es"
            else:
                return f"{count} {stem}s"
    else:
        raise ValueError(f"{count} must be an `int`!")


def was_were(count: int) -> str:
    if isinstance(count, int):
        if count == 1:
            return "1 was"
        else:
            return f"{count} were"
    else:
        raise ValueError(f"{count} must be an `int`!")


class BaseViewer(metaclass=abc.ABCMeta):

    def __init__(self):
        self._environment = Environment(loader=PackageLoader("gpsea.view", "templates"))
        self._environment.filters["pluralize"] = pluralize
        self._environment.filters["was_were"] = was_were


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
        fout = None
        try:
            fout = open_text_io_handle_for_writing(fh)
            fout.write(self._html)
        finally:
            if should_close and fout is not None:
                fout.close()

    def _repr_html_(self) -> str:
        return self._html
