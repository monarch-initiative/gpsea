import typing
import pytest

from ._gt_classifiers import _build_count_to_cat, _build_ac_to_cat, _qc_partitions


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
    ac2cat = _build_count_to_cat(
        a_label="A",
        b_label="B",
        partitions=partitions,
    )

    for ac, expected_cat_name in ac2cat_name.items():
        categorization = ac2cat[ac]
        assert categorization.category.name == expected_cat_name


@pytest.mark.parametrize(
    "partitions, ac2cat_name",
    [
        (
            ((0,), (1,), (2,)),
            {
                0: "0 alleles",
                1: "1 allele",
                2: "2 alleles",
            },
        ),
        (
            ((0, 1), (2,)),
            {
                0: "0 alleles OR 1 allele",
                1: "0 alleles OR 1 allele",
                2: "2 alleles",
            },
        ),
        (
            ((0,), (1, 2)),
            {
                0: "0 alleles",
                1: "1 allele OR 2 alleles",
                2: "1 allele OR 2 alleles",
            },
        ),
        (
            ((1,), (2,)),
            {
                1: "1 allele",
                2: "2 alleles",
            },
        ),
    ],
)
def test__build_ac_to_cat(
    partitions: typing.Sequence[typing.Sequence[int]],
    ac2cat_name: typing.Mapping[int, str],
):
    ac2cat = _build_ac_to_cat(
        partitions=partitions,
    )

    for ac, expected_cat_name in ac2cat_name.items():
        categorization = ac2cat[ac]
        assert categorization.category.name == expected_cat_name


@pytest.mark.parametrize(
    "partitions, expected",
    [
        (
            ((0, 1), (1,)),
            "element 1 was present 2!=1 times",
        ),
        (((0,), (0,), (2,)),
            "partition (0,) was present 2!=1 times"),
        (
            ((0, 1), (1, 2,)),
            "element 1 was present 2!=1 times",
        ),
        (
            ((0, 1), (1, 0)),
            "element 0 was present 2!=1 times, element 1 was present 2!=1 times",
        ),
        (
            ((0, 1, 1), (1, 1, 2,)),
            "element 1 was present 4!=1 times",
        ),
        (
            ((0.1,), (1,), (2,)),
            "Each partition index must be a non-negative int",
        ),
        (
            ((0, 1, 2),),
            "At least 2 partitions must be provided",
        ),
    ],
)
def test__qc_partitions__explodes_when_expected(
    partitions: typing.Sequence[typing.Sequence[int]],
    expected: str,
):
    with pytest.raises(ValueError) as e:
        _qc_partitions(partitions=partitions)

    assert e.value.args == (expected,)
