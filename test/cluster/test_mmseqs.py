import pytest

from protools.cluster import mmseqs
from protools import utils
try:
    mmseqs.MMseqs2()
    mmseqs_available = True
    helper_msg = ""
except utils.CmdNotFoundError as e:
    mmseqs_available = False
    helper_msg = str(e)


@pytest.mark.skipif(not mmseqs_available, reason=helper_msg)
class TestMMseqs2():

    @pytest.mark.xfail(reason="Not implemented")
    def test_cluster():
        # TODO
        assert False
    