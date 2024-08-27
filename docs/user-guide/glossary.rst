.. _glossary:

========
Glossary
========

The glossary summarizes several frequently used concepts.

.. _true-path-rule:

True path rule
~~~~~~~~~~~~~~

The true path rule of ontologies states that an item (e.g. an individual) annotated with an ontology term
(e.g. `Focal-onset seizure <https://hpo.jax.org/browse/term/HP:0007359>`_)
is implicitly annotated with all its *ancestor* terms
(`Seizure <https://hpo.jax.org/browse/term/HP:0001250>`_,
`Abnormal nervous system physiology <https://hpo.jax.org/browse/term/HP:0012638>`_, ...).
Conversely, exclusion of a term (e.g. `Abnormal ventricular septum morphology <https://hpo.jax.org/browse/term/HP:0010438>`_)
implies exclusion of all its *descendants*
(`Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_,
`Ventricular septal aneurysm <https://hpo.jax.org/browse/term/HP:0030957>`_, ...).


.. _length-of-the-reference-allele:

Length of the reference allele
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The length of the REF allele corresponds to the length of the genomic region affected by the variant.

Let's show a few examples.

>>> from gpsea.model import VariantCoordinates
>>> from gpsea.model.genome import GRCh38
>>> chr1 = GRCh38.contig_by_name("chr1")

The length of the reference allele of a missense variant is 1
because the variant affects a 1-bp region spanned by the ``C`` nucleotide:

>>> missense = VariantCoordinates.from_vcf_literal(chr1, 1001, 'C', 'T')
>>> len(missense)
1

The length of a "small" deletion is the same as the length of the ref allele `str`:
(``'CCC'`` in the example below):

>>> deletion = VariantCoordinates.from_vcf_literal(chr1, 1001, 'CCC', 'C')
>>> len(deletion)
3

This is because the literal notation spells out the alleles.
However, this simple rule does not apply in symbolic notation.
Here, the REF length corresponds to the length of the allele region.

For instance, for the following structural variant::

    #CHROM   POS    ID   REF  ALT     QUAL   FILTER   INFO
    1        1001   .    C    <DEL>   6      PASS     SVTYPE=DEL;END=1003;SVLEN=-3

the length of the REF allele is `3`:

>>> sv_deletion = VariantCoordinates.from_vcf_symbolic(
...     chr1, pos=1001, end=1003,
...     ref='C', alt='<DEL>', svlen=-3,
... )
>>> len(sv_deletion)
3

because the deletion removes 3 base pairs at the coordinates :math:`[1001, 1003]`.


.. _change-length-of-an-allele:

Change length of an allele
~~~~~~~~~~~~~~~~~~~~~~~~~~

Change length is defined as the difference between lengths of the *alt* and *ref* alleles.
SNVs lead to change length of zero, deletions and insertions/duplications lead to negative
and positive change lengths, respectively. See the table below for examples.

==================  ==================  ===============
 Reference allele    Alternate allele    Change length
==================  ==================  ===============
 C                   T                    0
 CG                  AT                   0
 ACTG                A                   -3
 A                   AAT                  2
==================  ==================  ===============
