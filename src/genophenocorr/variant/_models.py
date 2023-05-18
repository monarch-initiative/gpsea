import abc
from ._variant_data import Variant, VariantCoordinates


class VariantCoordinateFinder(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def find_coordinates(self, item) -> VariantCoordinates:
        pass

class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        pass