import math

from dataclasses import dataclass


@dataclass(frozen=True)
class Survival:
    """
    Information regarding individual's survival.
    """

    value: float
    """
    The survival value, expressed as a `float`.

    The `value` must be finite and non-`NaN`.
    """

    is_censored: bool
    """
    `True` if the survival has been censored and `False` otherwise.
    """

    def __post_init__(self):
        assert math.isfinite(
            self.value
        ), f"`value` must be finite and non-NaN, but was {self.value}"
