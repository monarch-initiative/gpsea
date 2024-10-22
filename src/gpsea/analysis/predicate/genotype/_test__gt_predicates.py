import typing
import pytest

from ._gt_predicates import _build_count_to_cat


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
        a_label="A", b_label="B",
        partitions=partitions,
    )

    for ac, expected_cat_name in ac2cat_name.items():
        categorization = ac2cat[ac]
        assert categorization.category.name == expected_cat_name
