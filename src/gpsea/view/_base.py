import abc

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
