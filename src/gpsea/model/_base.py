import enum
import typing

import hpotk


class Sex(enum.Enum):
    """
    `Sex` represents typical “phenotypic sex”, as would be determined by a midwife or physician at birth.
    
    The definition is aligned with `Phenopacket Schema <https://phenopacket-schema.readthedocs.io/en/2.0.0/sex.html>`_
    """
    
    UNKNOWN_SEX = 0
    """
    Not assessed or not available. Maps to ``NCIT:C17998``.
    """
    
    FEMALE = 1
    """
    Female sex. Maps to ``NCIT:C46113``.
    """

    MALE = 2
    """
    Male sex. Maps to ``NCIT:C46112``.
    """

    def is_provided(self) -> bool:
        """
        Return `True` if the sex is a known value, such as `FEMALE` or `MALE`.
        """
        return self != Sex.UNKNOWN_SEX

    def is_unknown(self) -> bool:
        """
        Return `True` if this is an `UNKNOWN_SEX`.
        """
        return self == Sex.UNKNOWN_SEX

    def is_female(self) -> bool:
        """
        Return `True` if the sex represents a `FEMALE`.
        """
        return self == Sex.FEMALE

    def is_male(self) -> bool:
        """
        Return `True` if the sex represents a `MALE`.
        """
        return self == Sex.MALE


class SampleLabels:
    """
    A data model for subject identifiers.

    The subject has a mandatory :attr:`label` and an optional :attr:`meta_label`.

    The identifiers support natural ordering, equality tests, and are hashable.
    """

    def __init__(self, label: str,
                 meta_label: typing.Optional[str] = None):
        self._label = hpotk.util.validate_instance(label, str, 'label')
        self._meta_label = hpotk.util.validate_optional_instance(meta_label, str, 'meta_label')

    @property
    def label(self) -> str:
        return self._label

    @property
    def meta_label(self) -> typing.Optional[str]:
        return self._meta_label

    def label_summary(self) -> str:
        """
        Summarize `label` and `meta_label` into a `str` where the sub-parts are inserted as ``<label>[<meta_label>]``.
        """
        return self._label if self._meta_label is None else f'{self._label}[{self._meta_label}]'

    def __eq__(self, other):
        return isinstance(other, SampleLabels) and self._label == other.label and self._meta_label == other._meta_label

    def __lt__(self, other):
        if isinstance(other, SampleLabels):
            if self._label < other._label:
                return True
            elif self._label == other._label:
                if self._meta_label is None or other._meta_label is None:
                    if self._meta_label == other._meta_label:
                        return False
                    else:
                        return True if self._meta_label is None else False  # `None` is less
                else:
                    return self._meta_label < other._meta_label
            return False
        else:
            return NotImplemented

    def __hash__(self):
        return hash((self._label, self._meta_label))

    def __str__(self):
        return self.label_summary()

    def __repr__(self):
        return f'SampleLabels(label={self._label}, meta_label={self._meta_label})'
