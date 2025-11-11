import pytest

from protools import pdbanno, utils, pdbio

try:
    pdbanno.TMalign()
    tm_available = True
    helper_msg = ""
except utils.CmdNotFoundError as e:
    tm_available = False
    helper_msg = str(e)

@pytest.mark.skipif(not tm_available, reason=helper_msg)
def test_tmalign():
    tmalign = pdbanno.TMalign()
    v = tmalign("data/gen.pdb", "data/native.pdb", 'data', 'aligned_gen')
    assert v == (0.95422, 0.99)


@pytest.mark.skipif(not tm_available, reason=helper_msg)
def test_batch_tmalign_paral(tmp_path):

    tmalign = pdbanno.TMalign(num_workers=5)
    data_iter = [("data/gen.pdb", "data/native.pdb")] * 5
    res = tmalign.batch_align(data_iter, tmp_path)
    # all tmscore should be 0.95422
    assert all(res['tmscore'] == 0.95422)
    # all rmsd should be 0.99
    assert all(res['rmsd'] == 0.99)
    # dimension of res should be 5
    assert res.shape[0] == 5


@pytest.mark.skipif(not tm_available, reason=helper_msg)
def test_batch_tmalign(tmp_path):
    tmalign = pdbanno.TMalign(num_workers=1)
    data_iter = [("data/gen.pdb", "data/native.pdb")] * 5
    res = tmalign.batch_align(data_iter, tmp_path)
    # all tmscore should be 0.95422
    assert all(res['tmscore'] == 0.95422)
    # all rmsd should be 0.99
    assert all(res['rmsd'] == 0.99)
    # dimension of res should be 5
    assert res.shape[0] == 5


@pytest.mark.skipif(not tm_available, reason=helper_msg)
def test_batch_tmalign_async(tmp_path):
    import asyncio
    loop = asyncio.get_event_loop()
    tmalign = pdbanno.TMalign(num_workers=1)
    data_iter = [("data/gen.pdb", "data/native.pdb")] * 5
    res = loop.run_until_complete(tmalign.async_batch_align(data_iter, tmp_path))
    # all tmscore should be 0.95422
    assert all(res['tmscore'] == 0.95422)
    # all rmsd should be 0.99
    assert all(res['rmsd'] == 0.99)
    # dimension of res should be 5
    assert res.shape[0] == 5



@pytest.mark.parametrize(
        'pdb_file',
        ['data/3inj.pdb', 'data/4I77.pdb']
)
def test_neighbor_water_count(pdb_file: str):
    pdbanno.neighbor_water_count(pdbio.get_structure(pdb_file))


@pytest.mark.parametrize(
        'pdb_file',
        ['data/3inj.pdb', 'data/4I77.pdb']
)
def test_calc_sasa(pdb_file: str):
    pdbanno.calc_sasa(pdbio.get_structure(pdb_file))