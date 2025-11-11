import pytest

from protools import seqplot


def test_aa_heatmap():
    seqplot.aa_heatmap(['ACGDE', 'ACGDE'])