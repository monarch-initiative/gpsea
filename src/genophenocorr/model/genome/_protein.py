import hpotk.util

from ._genome import Region


class ProteinCoordinates:

    def __init__(self, identifier: str,
                 region: Region):
        self._id = hpotk.util.validate_instance(identifier, str, 'identifier')
        self._region = hpotk.util.validate_instance(region, Region, 'region')

    @property
    def identifier(self) -> str:
        return self._id

    @property
    def region(self) -> Region:
        return self._region

    def __eq__(self, other):
        return (isinstance(other, ProteinCoordinates)
                and self.identifier == other.identifier
                and self.region == other.region)

    def __str__(self):
        return f'ProteinCoordinates(identifier={self.identifier}, region={self.region})'

    def __repr__(self):
        return str(self)
