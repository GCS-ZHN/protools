import numpy as np
from typing import Iterable, Optional, List


def random_split(
    data: Iterable,
    lengths: Iterable[int],
    seed: Optional[int] = None
) -> List[List]:
    """
    Randomly split a dataset into non-overlapping new datasets of given lengths.

    Parameters
    ----------
    data : Iterable
        The dataset to be split.
    lengths : Iterable[int]
        A list of lengths for each split. The sum of lengths must not exceed the size of the dataset.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    list
        A list of datasets, each corresponding to a split of the original dataset.
    """
    rng = np.random.default_rng(seed)
    data = list(data)
    total_length = sum(lengths)
    
    if total_length > len(data):
        raise ValueError("Sum of lengths exceeds the size of the dataset.")
    
    indices = rng.permutation(len(data))[:total_length]
    splits = []
    current_index = 0
    
    for length in lengths:
        split_indices = indices[current_index:current_index + length]
        splits.append([data[i] for i in split_indices])
        current_index += length
    
    return splits


def grouped_split(
    group_sizes: Iterable[int],
    lengths: Iterable[int],
    mode: str = 'partial_random',
    seed: Optional[int] = None,
    chunk_size: int = 3
) -> List[List[int]]:
    """
    Split a grouped dataset into non-overlapping new datasets of given lengths.
    Algorithm will try to minimize the difference between result lengths and
    given lengths.

    Parameters
    ----------
    group_sizes : Iterable[int]
        A list of sizes for each group in the dataset.

    lengths : Iterable[int]
        A list of lengths for each split. The sum of lengths must be equal to the total size of input dataset.

    mode : str, optional
        The mode of splitting. Options are 'partial_random', 'random', 'sorted'.
        Default is 'partial_random'.
        - 'partial_random': dataset will sorted by size and shuffled for each chunk.
        - 'random': dataset will be shuffled randomly.
        - 'sorted': dataset will be sorted by size.

    seed : int, optional
        Random seed for reproducibility.

    chunk_size : int, optional
        The size of chunks for 'partial_random' mode. Default is 3.

    Returns
    -------
    list
        A list of datasets, each corresponding to a split of the original dataset.
    """

    group_sizes = list(group_sizes)
    groups = np.arange(len(group_sizes), dtype=int)
    rng = np.random.default_rng(seed)
    lengths = list(lengths)
    assert sum(lengths) == sum(group_sizes), \
        "Sum of lengths must be equal to the total size of the dataset."
    
    assert mode in ('partial_random', 'random', 'sorted'), \
        "Mode must be one of 'partial_random', 'random', 'sorted'."
    if mode in ('partial_random', 'sorted'):
        groups = np.argsort(group_sizes)[::-1]
    else:
        groups = rng.permutation(len(group_sizes))

    if mode == 'partial_random':
        if chunk_size <= 0:
            raise ValueError("chunk_size must be a positive integer.")
        if chunk_size >= len(groups):
            raise ValueError("chunk_size must be smaller than the number of groups.")

        groups = np.concatenate([
            rng.permutation(groups[i:i + chunk_size])
            for i in range(0, len(groups), chunk_size)
        ])

    splits = [[] for _ in lengths]
    split_sizes = [0 for _ in lengths]

    for group in groups:
        group_size = group_sizes[group]
        best_split = np.argmin([
            (length - (split_size + group_size)) ** 2
            if split_size + group_size <= length else float('inf')
            for split_size, length in zip(split_sizes, lengths)
        ])
        splits[best_split].append(group)
        split_sizes[best_split] += group_size
    return splits
        