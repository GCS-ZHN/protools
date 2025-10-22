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