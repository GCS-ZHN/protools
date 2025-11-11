import pytest

from protools import seqio

class TestFasta:
    
    def test_fasta_getitem_str(self):
        fasta = seqio.Fasta()
        fasta['seq1'] = 'ACDEFGHIK'
        assert str(fasta['seq1'].seq) == 'ACDEFGHIK'

    def test_fasta_getitem_int(self):
        fasta = seqio.Fasta()
        fasta['seq1'] = 'ACDEFGHIK'
        fasta['seq2'] = 'LMNPQRSTV'
        assert str(fasta[0].seq) == 'ACDEFGHIK'
        assert str(fasta[1].seq) == 'LMNPQRSTV'

    def test_fasta_getitem_slice(self):
        fasta = seqio.Fasta()
        fasta['seq1'] = 'ACDEFGHIK'
        fasta['seq2'] = 'LMNPQRSTV'
        fasta['seq3'] = 'WYACDEFGH'
        sliced = fasta[1:]
        assert str(sliced[0].seq) == 'LMNPQRSTV'
        assert str(sliced[1].seq) == 'WYACDEFGH'

    def test_fasta_setitem_non_str_key(self):
        fasta = seqio.Fasta()
        with pytest.warns(UserWarning, match='Key 123 is not a string'):
            fasta[123] = 'ACDEFGHIK'
        assert str(fasta['123'].seq) == 'ACDEFGHIK'

    def test_fasta_setitem_str_value(self):
        fasta = seqio.Fasta()
        fasta['seq1'] = 'ACDEFGHIK'
        assert str(fasta['seq1'].seq) == 'ACDEFGHIK'

    def test_fasta_setitem_seq_value(self):
        from Bio.Seq import Seq
        fasta = seqio.Fasta()
        fasta['seq1'] = Seq('ACDEFGHIK')
        assert str(fasta['seq1'].seq) == 'ACDEFGHIK'

    def test_fasta_setitem_seqrecord_err(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        with pytest.raises(AssertionError, match='Mismatch id and SeqRecord'):
            fasta = seqio.Fasta()
            fasta['seq1'] = SeqRecord(Seq('ACDEFGHIK'))
            assert str(fasta['seq1'].seq) == 'ACDEFGHIK'


def test_read_fasta():
    fasta = seqio.read_fasta('data/test.fasta')
    assert len(fasta) == 3
    assert str(fasta['seq1'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta['seq2'].seq) == 'WYACDEFGHIKLMNPQRSTV'
    assert str(fasta['seq3'].seq) == 'RSTVWYACDEFGHIKLMNPQ'


def test_read_fasta_with_gz():
    fasta = seqio.read_fasta('data/test.fasta.gz')
    assert len(fasta) == 3
    assert str(fasta['seq1'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta['seq2'].seq) == 'WYACDEFGHIKLMNPQRSTV'
    assert str(fasta['seq3'].seq) == 'RSTVWYACDEFGHIKLMNPQ'


def test_read_a3m():
    msa = seqio.read_a3m('data/ab1218_0.a3m')
    assert len(msa) == 9933


def test_read_multimer_a3m():
    msa = seqio.read_a3m('data/ab1218.a3m', is_multimer=True)
    assert len(msa) == 3
    assert len(msa[0]) == 9933
    assert len(msa[1]) == 12405
    assert len(msa[2]) == 401


def test_save_fasta(tmp_path):
    fasta = seqio.Fasta()
    fasta['seq1'] = 'ACDEFGHIKLMNPQRSTVWY'
    fasta['seq2'] = 'WYACDEFGHIKLMNPQRSTV'
    out_file = tmp_path / 'out.fasta'
    seqio.save_fasta(fasta.values(), out_file)
    fasta2 = seqio.read_fasta(out_file)
    assert len(fasta2) == 2
    assert str(fasta2['seq1'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta2['seq2'].seq) == 'WYACDEFGHIKLMNPQRSTV'


def test_save_fasta_with_gz(tmp_path):
    fasta = seqio.Fasta()
    fasta['seq1'] = 'ACDEFGHIKLMNPQRSTVWY'
    fasta['seq2'] = 'WYACDEFGHIKLMNPQRSTV'
    out_file = tmp_path / 'out.fasta.gz'
    seqio.save_fasta(fasta.values(), out_file)
    fasta2 = seqio.read_fasta(out_file)
    assert len(fasta2) == 2
    assert str(fasta2['seq1'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta2['seq2'].seq) == 'WYACDEFGHIKLMNPQRSTV'


def test_dataframe2fasta(tmp_path):
    import pandas as pd
    df = pd.DataFrame({
        'id': ['seq1', 'seq2'],
        'sequence': ['ACDEFGHIKLMNPQRSTVWY', 'WYACDEFGHIKLMNPQRSTV']
    })
    out_file = tmp_path / 'out.fasta'
    seqio.df2fasta(df, out_file, 'sequence', id_col='id')
    fasta = seqio.read_fasta(out_file)
    assert len(fasta) == 2
    assert str(fasta['seq1'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta['seq2'].seq) == 'WYACDEFGHIKLMNPQRSTV'


def test_dataframe2fasta_no_id_col(tmp_path):
    import pandas as pd
    df = pd.DataFrame({
        'sequence': ['ACDEFGHIKLMNPQRSTVWY', 'WYACDEFGHIKLMNPQRSTV']
    })
    out_file = tmp_path / 'out.fasta'
    seqio.df2fasta(df, out_file, 'sequence')
    fasta = seqio.read_fasta(out_file)
    assert len(fasta) == 2
    assert str(fasta['0'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta['1'].seq) == 'WYACDEFGHIKLMNPQRSTV'


def test_dataframe2fasta_with_gz(tmp_path):
    import pandas as pd
    df = pd.DataFrame({
        'id': ['seq1', 'seq2'],
        'sequence': ['ACDEFGHIKLMNPQRSTVWY', 'WYACDEFGHIKLMNPQRSTV']
    })
    out_file = tmp_path / 'out.fasta.gz'
    seqio.df2fasta(df, out_file, 'sequence', id_col='id')
    fasta = seqio.read_fasta(out_file)
    assert len(fasta) == 2
    assert str(fasta['seq1'].seq) == 'ACDEFGHIKLMNPQRSTVWY'
    assert str(fasta['seq2'].seq) == 'WYACDEFGHIKLMNPQRSTV'


@pytest.mark.xfail(reason="Not implemented")
def test_read_pdb_seqres():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_read_mmcif_seqres():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_read_seqres():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_temp_fasta():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_create_complex_seq():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_cross_create():
    # TODO
    assert False
