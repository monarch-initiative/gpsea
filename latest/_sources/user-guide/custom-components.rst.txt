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

.. _custom-variant-predicate:

*****************
Variant predicate
*****************

The purpose of a :class:`~gpsea.analysis.predicate.VariantPredicate` is to test
if a variant meets a certain criterion and GPSEA ships with an array
of builtin predicates (see :mod:`gpsea.analysis.predicate` module).
However, chances are a custom predicate will be needed in future,
so we show how to how to extend
the :class:`~gpsea.analysis.predicate.VariantPredicate` class
to create one's own predicate.

Specifically, we show how to create a predicate to test if the variant affects a glycine residue
of the transcript of interest.

>>> from gpsea.model import Variant, VariantEffect
>>> from gpsea.analysis.predicate import VariantPredicate
>>> class AffectsGlycinePredicate(VariantPredicate): # ❶
...     def __init__( # ❷
...         self,
...         tx_id: str,
...     ):
...         self._tx_id = tx_id
...         self._aa_code = "Gly"
... 
...     @property
...     def name(self) -> str: # ❸
...         return "Affects Glycine"
...    
...     @property
...     def description(self) -> str: # ❹
...         return "affects a glycine residue"
...    
...     @property
...     def variable_name(self) -> str: # ❺
...         return "affected aminoacid residue"
...     
...     def test(self, variant: Variant) -> bool: # ❻
...         tx_ann = variant.get_tx_anno_by_tx_id(self._tx_id)
...         if tx_ann is not None:
...            hgvsp = tx_ann.hgvsp
...            if hgvsp is not None:
...                return hgvsp.startswith(f"p.{self._aa_code}")
...         return False
...     
...     def __eq__(self, value: object) -> bool: # ➐
...         return isinstance(value, AffectsGlycinePredicate) and self._tx_id == value._tx_id
...     
...     def __hash__(self) -> int: # ❽
...         return hash((self._tx_id,))
...     
...     def __repr__(self) -> str: # ❾
...         return str(self)
...     
...     def __str__(self) -> str: # ➓
...         return f"AffectsGlycinePredicate(tx_id={self._tx_id})"

❶ The ``AffectsGlycinePredicate`` must extend :class:`~gpsea.analysis.predicate.VariantPredicate`.
❷ We ask the user to provide the transcript accession `str` and we set the target aminoacid code to glycine ``Gly``.
Like in the :ref:`custom-phenotype-scorer` above, ❸❹❺ provide metadata required for the bookkeeping.
The ❻ ``test`` method includes the most interesting part - we retrieve the :class:`~gpsea.model.TranscriptAnnotation`
with the functional annotation data for the transcript of interest, and we test if the HGVS protein indicates
that the reference aminoacid is glycine.
Last, we override ➐ ``__eq__()`` and ❽ ``__hash__()`` (required) as well as ❾ ``__repr__()`` and ➓ ``__str__()`` (recommended).
