import pytest

from protools import dock, utils, pdbio


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


@pytest.mark.parametrize(
        'model_file, native_file, chain_map, dockq_label', 
        [
            ('data/4i77.pdb', 'data/4i77.pdb', {'H': 'H', 'L': 'L', 'Z': 'Z'}, 'High')
        ]
)
def test_dockq_score(model_file, native_file, chain_map, dockq_label):
    model = pdbio.get_structure(model_file)[0]
    native = pdbio.get_structure(native_file)[0]
    dockqs = dock.dockq_score(model, native, chain_map)
    assert dockqs['DockQ_label'] == dockq_label
