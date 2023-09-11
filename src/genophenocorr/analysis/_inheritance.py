import abc
import typing

from genophenocorr.model import Variant


class InheritanceCompatibilityChecker(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def check_compatibility(self, patient: typing.Sequence[Variant]) -> bool:
        pass


class AutosomalDominantChecker(InheritanceCompatibilityChecker):

    def check_compatibility(self, patient: typing.Sequence[Variant]) -> bool:
        pass


class AutosomalRecessiveChecker(InheritanceCompatibilityChecker):

    def check_compatibility(self, patient: typing.Sequence[Variant]) -> bool:
        pass


class MitochondrialChecker(InheritanceCompatibilityChecker):

    def check_compatibility(self, patient: typing.Sequence[Variant]) -> bool:
        pass
