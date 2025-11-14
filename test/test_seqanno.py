import pytest

from protools import seqanno


@pytest.mark.parametrize(
    "s1, s2, identity",
    [
        ('QVQ', 'QVQ', 1.0),
        ('QV', 'QL', 0.5)
    ])
def test_calc_seq_identity(s1: str, s2: str, identity: float):
    assert seqanno.calc_seq_identity(s1, s2) == identity


@pytest.mark.xfail(reason="Not implemented")
def test_anno_cdr():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_numbering_seq():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_anno_vj_gene():
    # TODO
    assert False
