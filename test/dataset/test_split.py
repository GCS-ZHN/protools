import pytest

from protools.dataset import split
from test import tools


def test_random_split():
    splits = split.random_split(range(100), [30, 30, 40])
    assert len(splits[0]) == 30
    assert len(splits[1]) == 30
    assert len(splits[2]) == 40


def test_random_split_with_seed():
    splits = split.random_split(range(2), [1, 1], seed=42)
    assert splits[0][0] == 1
    assert splits[1][0] == 0


def test_random_split_with_over_length():
    with pytest.raises(
        ValueError,
        match='Sum of lengths exceeds the size of the dataset.'):
        split.random_split(range(2), [5, 5])


def test_random_split_with_ilegal_length():
    with pytest.raises(
        ValueError,
        match='positive integer'):
        split.random_split(range(2), [1, -1])


def test_grouped_split_with_sorted():
    splits = split.grouped_split(
        group_sizes=[4, 3, 2, 1],
        lengths=[8, 2],
        mode='sorted'
    )
    assert splits[0] == [0, 1, 3]
    assert splits[1] == [2]


def test_grouped_split_with_random():
   splits_all = [split.grouped_split(
        group_sizes=[4, 3, 2, 1],
        lengths=[8, 2],
        mode='random'
    ) for _ in range(10)]
   assert not tools.all_equal(splits_all)


def test_grouped_split_with_partial_random():
   splits_all = [split.grouped_split(
        group_sizes=[4, 3, 2, 1],
        lengths=[8, 2],
        mode='partial_random',
        chunk_size=3
    ) for _ in range(10)]
   assert not tools.all_equal(splits_all)
   assert tools.all_equal(s[1] for s in splits_all)


def test_grouped_split_with_min_chunk():
   splits_all = [split.grouped_split(
        group_sizes=[4, 3, 2, 1],
        lengths=[8, 2],
        mode='partial_random',
        chunk_size=1
    ) for _ in range(10)]
   assert tools.all_equal(splits_all)
   sorted_splits = split.grouped_split(
        group_sizes=[4, 3, 2, 1],
        lengths=[8, 2],
        mode='sorted')
   assert sorted_splits == splits_all[0]
