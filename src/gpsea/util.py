import io
import platform
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


def open_text_io_handle_for_writing(
    file: typing.Union[io.IOBase, str],
    encoding: str = "utf-8",
) -> io.TextIOWrapper:
    """
    Open a `io.TextIO` file handle based on `file`.

    :param file: a `str` or :class:`io.IOBase` to write to. If `str`, then the input is interpreted as a path
      to a local file.
      If `fh` is an IO wrapper, the function ensures we get a text wrapper that uses given encoding.
    :param encoding: encoding used to encode the output (`utf-8` by default).
    :return: the :class:`io.TextIOWrapper` wrapper.
    """
    if isinstance(file, str):
        # A path to local file
        return open(file, "w", encoding=encoding)
    elif isinstance(file, io.IOBase):
        if isinstance(file, (io.TextIOWrapper, io.TextIOBase)):
            return file
        elif isinstance(file, (io.BytesIO, io.BufferedIOBase)):
            return io.TextIOWrapper(file, encoding=encoding)

    raise ValueError(f"Unsupported type {type(file)}")


def open_resource(package: str, path: str) -> typing.IO[str]:
    """
    Get an IO handle for reading a package resource on the current platform.

    The returned callable expects two arguments:
    :param package: a `str` with package name (e.g. `gpsea.model.genome`)
    :param path: a `str` with the file name that is located in the package (e.g. `GCF_000001405.25_GRCh37.p13_assembly_report.tsv`)
    """
    major, minor, patch = platform.python_version_tuple()

    if major == "3":
        import importlib.resources

        minor = int(minor)
        if minor < 9:
            # Versions 3.7, 3.8
            return importlib.resources.open_text(package, path)

        else:
            return importlib.resources.files(package).joinpath(path).open()

    else:
        raise ValueError(f"Untested Python version v{major}.{minor}.{patch}")
