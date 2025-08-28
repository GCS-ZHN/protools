# -*- coding: utf-8 -*-
# @Author  : GCS-ZHN
# @Time    : 2025-08-28 17:55:14

"""
This module provides functions to visualize amino acid sequences.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from typing import Iterable, Optional
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from protools.utils import Intervals
from Bio.Data.IUPACData import protein_letters


__all__ = ['aa_heatmap']


def aa_heatmap(
        seqs: Iterable[str],
        ref_seq: Optional[str] = None,
        cmap_start_color: str = '#FFFFFF',
        cmap_end_color: str = '#4B0082',
        ref_box_color: str = 'red',
        display_intervals: Optional[Intervals] = None,
        title: Optional[str] = None,
        ax: Optional[plt.Axes] = None,
        height: float = 10.0) -> plt.Axes:
    """
    Draw a heatmap of amino acid sequences.

    Parameters
    ----------
    seqs : Iterable[str]
        A list of amino acid sequences with the same length to be visualized.
    ref_seq : Optional[str], default=None
        A reference sequence to highlight in the heatmap.
    cmap_start_color : str, default='#FFFFFF'
        The starting color of the colormap.
    cmap_end_color : str, default='#B0082'
        The ending color of the colormap.
    ref_box_color : str, default='red'
        The color of the box to highlight the reference sequence.
    display_intervals : Optional[Intervals], default=None
        Intervals to display. If None, display all positions.
    title : Optional[str], default=None
        Title of the heatmap.
    ax : Optional[plt.Axes], default=None
        Matplotlib Axes to plot on. If None, a new figure and axes will be
        created.
    height : float, default=10.0
        Height of the figure if a new figure is created.
    """

    seqs = list(seqs)
    if not seqs:
        raise ValueError("The input sequence list is empty.")
    
    if not all(len(seq) == len(seqs[0]) for seq in seqs):
        raise ValueError("All sequences must have the same length.")
    
    if ref_seq is not None and len(ref_seq) != len(seqs[0]):
        raise ValueError("Reference sequence must have the same length as the input sequences.")
    
    positions = list(range(1, len(seqs[0]) + 1))
    if display_intervals is not None:
        assert not display_intervals.zero_based, \
                "display_intervals must be one-based, not zero-based."
        positions = [p for p in positions if p in display_intervals]

    counts = np.zeros((len(protein_letters), len(positions)), dtype=int)
    aa2idx = {aa: idx for idx, aa in enumerate(protein_letters)}
    for seq in seqs:
        for i, pos in enumerate(positions):
            aa = seq[pos - 1].upper()
            if aa in aa2idx:
                counts[aa2idx[aa], i] += 1

    probs = counts / len(seqs)
    probs = pd.DataFrame(
        probs,
        index=list(protein_letters),
        columns=positions
    )

    custom_cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', [cmap_start_color, cmap_end_color]
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=(len(positions) * height / 20, height))
    else:
        height = ax.figure.get_size_inches()[1]

    ax = sns.heatmap(
        probs,
        cmap=custom_cmap,
        cbar_kws={'label': 'Frequency', 'shrink': 0.5},
        square=True,
        ax=ax
    )

    if ref_seq:
        for i, pos in enumerate(positions):
            aa = ref_seq[pos - 1].upper()
            if aa in aa2idx:
                ax.add_patch(Rectangle(
                    (i, aa2idx[aa]),
                    1, 1,
                    fill=False,
                    edgecolor=ref_box_color,
                    lw=2 * height / 10
                ))
    ax.set_xlabel('Position')
    ax.set_ylabel('Amino Acid')
    if title:
        ax.set_title(title)
    plt.tight_layout()
    return ax