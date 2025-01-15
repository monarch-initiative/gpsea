.. _custom-components:

#################
Custom components
#################

GPSEA aims to stay useful in the long run. Therefore, we took a great care
to adhere to "good software development practices" and we framed GPSEA functionalities
as a set of loosely coupled components. As a rule of thumb, each component corresponds
to a Python abstract base class which is then extended by the builtin components
and can also be extended by the future components, to serve both common or exotic use cases.

The abstract base classes define the component API.
Per guidelines in Python's :mod:`abc` module, the abstract classes use :class:`abc.ABCMeta` as a metaclass
and the class API consists of methods annotated with the :func:`abc.abstractmethod` decorator.
These decorated methods must be overridden in the subclasses.

The following sections provide guidance for customizing the most commonly used GPSEA components.


.. _custom-phenotype-scorer:

****************
Phenotype scorer
****************

:class:`~gpsea.analysis.pscore.PhenotypeScorer` computes a phenotype score for an individual.
The phenotype score is a `float` with range :math:`(-\infty, \infty)` where `NaN` indicates
that a score cannot be computed (e.g. the lab measurement value was not ascertained for the individual).

Here we show an example of a toy phenotype scorer
for using body mass index (BMI) as a phenotype score.

>>> import typing
>>> from gpsea.model import Patient
>>> from gpsea.analysis.pscore import PhenotypeScorer
>>> class BmiScorer(PhenotypeScorer):  # ❶
...    
...     def __init__( # ❷
...         self,
...         id2bmi: typing.Mapping[str, float],
...     ):
...         self._id2bmi = id2bmi
...     
...     @property
...     def name(self) -> str: # ❸
...         return "BMI phenotype scorer"
...    
...     @property
...     def description(self) -> str: # ❹
...         return "Body mass index used as a phenotype score"
...    
...     @property
...     def variable_name(self) -> str: # ❺
...         return "BMI"
...     
...     def score(self, patient: Patient) -> float: # ❻
...         try:
...             return self._id2bmi[patient.labels.label]
...         except KeyError:
...             return float('nan')

❶ The ``BmiScorer`` must extend :class:`~gpsea.analysis.pscore.PhenotypeScorer`
to be used as a phenotype scorer.
❷ The scorer needs a ``dict`` with `label` → `BMI` for the analyzed individuals.
We assume the user will pre-compute the corresponding ``dict``.

Then, the scorer must expose several properties, including ❸ ``name``, ❹ ``description``,
and the ❺ ``variable_name`` it operates on.
The properties provide bookkeeping metadata to use in e.g. visualizations.
Try to choose short and concise names.

The most important part of the scorer is the ❻ `score` method
which retrieves the BMI for an individual or returns `NaN` if the value is not available
and the individual should be omitted from the analysis.

