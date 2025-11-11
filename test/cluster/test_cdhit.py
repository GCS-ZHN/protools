import pytest

from protools.cluster import cdhit
from protools import utils
try:
    cdhit.CdHit()
    cdhit_available = True
    helper_msg = ""
except utils.CmdNotFoundError as e:
    cdhit_available = False
    helper_msg = str(e)


@pytest.mark.skipif(not cdhit_available, reason=helper_msg)
class TestCdHit():

    @pytest.mark.xfail(reason="Not implemented")
    def test_call():
        # TODO
        assert False
    