from pathlib import Path
import pytest

from protools.cluster import mmseqs
from protools import utils, seqio
try:
    mmseqs.MMseqs2()
    mmseqs_available = True
    helper_msg = ""
except utils.CmdNotFoundError as e:
    mmseqs_available = False
    helper_msg = str(e)


@pytest.mark.skipif(not mmseqs_available, reason=helper_msg)
class TestMMseqs2():

    def setup_method(self):
        self.cmd = mmseqs.MMseqs2()

    def test_cluster(self, tmp_path: Path):
        data = seqio.read_fasta('data/test.fasta')
        result = self.cmd.cluster(
            seqs=data,
            output=tmp_path / 'result.csv',
            tmp_dir=tmp_path
        )
        assert len(result) == 3
        assert (result.columns == ['cluster_id', 'sequence_id']).all()
