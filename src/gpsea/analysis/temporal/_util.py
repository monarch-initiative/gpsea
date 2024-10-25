import typing

from scipy.stats import CensoredData

from ._api import Survival


def prepare_censored_data(
    survivals: typing.Iterable[Survival],
) -> CensoredData:
    uncensored = []
    right_censored = []
    for survival in survivals:
        if survival.is_censored:
            right_censored.append(survival.value)
        else:
            uncensored.append(survival.value)
    return CensoredData(
        uncensored=uncensored,
        right=right_censored,
    )
