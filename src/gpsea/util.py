import io
import typing


def open_text_io_handle_for_reading(
    file: typing.Union[io.IOBase, str],
    encoding: str = "utf-8",
) -> io.TextIOWrapper:
    """
    Open a `io.TextIO` file handle based on `file`.

    :param file: a `str` or :class:`io.IOBase` to read from. If `str`, then the input is interpreted as a path
      to a local file.
      If `fh` is an IO wrapper, the function ensures we get a text wrapper that uses given encoding.
    :param encoding: encoding used to decode the input (`utf-8` by default).
    :return: the :class:`io.TextIOWrapper` wrapper.
    """
    if isinstance(file, str):
        # A path to local file
        return open(file, "r", encoding=encoding)
    elif isinstance(file, io.IOBase):
        if isinstance(file, (io.TextIOWrapper, io.TextIOBase)):
            return file
        elif isinstance(file, (io.BytesIO, io.BufferedIOBase)):
            return io.TextIOWrapper(file, encoding=encoding)

    raise ValueError(f"Unsupported type {type(file)}")
