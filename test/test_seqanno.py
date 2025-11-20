from typing import Callable, List
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


@pytest.mark.parametrize(
        's1, s2, comparsion, chain_id, exp_mutations',
        [
            ('QVQ', 'QVQ', 'type', 'H', []),
            ('EVQ', 'DVQ', 'property', 'H', []),
            ('QVQE', 'EVDE', 'type', '', ['Q1E', 'Q3D']),
            ('-VQ', 'EVQ', 'type', 'H', ['-H1E']),
            ('QVD', 'EYG', lambda x, y: True, '', []),
            ('QVQ', 'QVQ', lambda x, y: True, '', [])
        ]
)
def test_get_mutations(
    s1: str, s2: str, comparsion: str|Callable[[str, str], bool],
    chain_id: str, exp_mutations: List[str]):
    assert seqanno.get_mutations(s1, s2, comparsion, chain_id) == exp_mutations
