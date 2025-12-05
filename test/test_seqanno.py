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



@pytest.fixture(name='vdomain_data', scope='module')
def fixture_vdomain_data():
    from protools import seqio
    import pandas as pd
    seqs = seqio.read_fasta('data/vdomain.fasta')
    annotations = pd.read_csv('data/vdomain_anno.csv', index_col=[0, 1, 2])
    return seqs, annotations


@pytest.mark.parametrize(
        'seq_id, chain, scheme, is_tcr',
        [
            ('1iqd_B', 'H', 'imgt', False),
            ('1iqd_B', 'H', 'kabat', False),
            ('1iqd_B', 'H', 'chothia', False),
            ('1iqd_B', 'H', 'martin', False),
            ('1dfb_H','H', 'imgt', False),
            ('1dfb_H','H', 'kabat', False),
            ('1dfb_H','H', 'chothia', False),
            ('1dfb_H','H', 'martin', False)
        ]
)
def test_anno_cdr(vdomain_data, seq_id, chain, scheme, is_tcr):
    seqs, exp_annotations = vdomain_data
    annotations = seqanno.anno_cdr(
        seqs[seq_id], chain=chain, scheme=scheme, is_tcr=is_tcr)
    exp_annotations = exp_annotations.loc[(seq_id, chain, scheme)]
    names = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    for name in names:
        assert annotations[chain + name] == exp_annotations[name]


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
            ('-DG', 'EVQ', 'type', 'H', ['-H1E', 'DH1V', 'GH2Q']),
            ('W-E-VQ', 'WGGEVE', 'type', '', ['-2G', 'E2G', '-3E', 'Q4E']),
            ('E-GA', 'EGQ-', 'type', '', ['-2G', 'G2Q', 'A3-']),
            ('QVD', 'EYG', lambda x, y: True, '', []),
            ('QVQ', 'QVQ', lambda x, y: True, '', [])
        ]
)
def test_get_mutations(
    s1: str, s2: str, comparsion: str|Callable[[str, str], bool],
    chain_id: str, exp_mutations: List[str]):
    assert seqanno.get_mutations(s1, s2, comparsion, chain_id) == exp_mutations
