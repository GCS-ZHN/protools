import pytest

from protools import dock, utils


try:
    dock.HDock()
    hdock_available = True
    helper_msg = ""
except utils.CmdNotFoundError as e:
    hdock_available = False
    helper_msg = str(e)


@pytest.mark.skipif(not hdock_available, reason=helper_msg)
class TestHDock():

    @pytest.mark.xfail(reason="Not implemented")
    def test_dock():
        # TODO
        assert False
    
    @pytest.mark.xfail(reason="Not implemented")
    def test_create_complex():
        # TODO
        assert False
    