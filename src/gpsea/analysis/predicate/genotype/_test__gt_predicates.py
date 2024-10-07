import typing
import pytest

from ._gt_predicates import (
    build_count_to_cat,
    ModeOfInheritancePredicate,
    ModeOfInheritanceInfo,
)


@pytest.mark.parametrize(
    "partitions, ac2cat_name",
    [
        (
            ((0,), (1,), (2,)),
            {
                (2, 0): "A/A",
                (1, 1): "A/B",
                (0, 2): "B/B",
            },
        ),
        (
            ((0, 1), (2,)),
            {
                (2, 0): "A/A OR A/B",
                (1, 1): "A/A OR A/B",
                (0, 2): "B/B",
            },
        ),
        (
            ((0,), (1, 2)),
            {
                (2, 0): "A/A",
                (1, 1): "A/B OR B/B",
                (0, 2): "A/B OR B/B",
            },
        ),
        (
            ((1,), (2,)),
            {
                (1, 1): "A/B",
                (0, 2): "B/B",
            },
        ),
    ],
)
def test_build_count_to_cat(
    partitions: typing.Sequence[typing.Sequence[int]],
    ac2cat_name: typing.Mapping[typing.Tuple[int, int], str],
):
    names = ("A", "B")
    ac2cat = build_count_to_cat(
        names=names,
        partitions=partitions,
    )

    for ac, expected_cat_name in ac2cat_name.items():
        categorization = ac2cat[ac]
        assert categorization.category.name == expected_cat_name


@pytest.mark.parametrize(
    "partitions, moi_data, ac2cat_name",
    [
        (
            [(0,), (1,)],
            ModeOfInheritanceInfo.autosomal_dominant(),
            {
                0: "No allele",
                1: "Monoallelic",
            },
        ),
        (
            [(0,), (1,), (2,)],
            ModeOfInheritanceInfo.autosomal_recessive(),
            {
                0: "No allele",
                1: "Monoallelic",
                2: "Biallelic",
            },
        ),
        (
            [(1,), (2,)],
            ModeOfInheritanceInfo.autosomal_recessive(),
            {
                1: "Monoallelic",
                2: "Biallelic",
            },
        ),
    ],
)
def test_prepare_count_to_cat(
    partitions: typing.Sequence[typing.Sequence[int]],
    moi_data: ModeOfInheritanceInfo,
    ac2cat_name: typing.Mapping[int, str],
):
    ac2cat = ModeOfInheritancePredicate.prepare_count2cat(
        partitions=partitions,
        mode_of_inheritance_data=moi_data,
    )

    for ac, expected_cat_name in ac2cat_name.items():
        categorization = ac2cat[ac]
        assert categorization.category.name == expected_cat_name
