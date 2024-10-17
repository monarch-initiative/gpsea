import typing

from gpsea.model import Patient
from ._api import Endpoint
from ._base import Survival


class Death(Endpoint):

    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        # If the patient is alive we use the current age as `value` and `censored=True`
        # If the patient is deceased, we use the age at death as a `value` and `censored=False`
        if patient.age_at_death is not None:
            return Survival(
                value=patient.age_at_death.days,
                is_censored=False,
            )
        raise NotImplementedError


class PhenotypicFeatureOnset(Endpoint):

    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        # If the individual has the target HPO feature:
        #  - if we know about the feature onset, then we use the onset as value and `censored=False`
        #  - else `None`
        # If the individual does not have the feature:
        #   - if the individual has current age, we use the current age and `censored=True` !! CHECK !!
        #   - else `None`
        raise NotImplementedError


class DiseaseOnset(Endpoint):

    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        # If the individual has the target disease:
        #   - we know about the disease onset, then we use the onset as value and `censored=False`
        #   - else `None`
        # If the individual does not have the disease:
        #   - if the individual has current age, we use the current age and `censored=True` !! CHECK !!
        #   - else `None`
        raise NotImplementedError
