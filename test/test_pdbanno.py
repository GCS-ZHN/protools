import pytest

from protools import pdbanno, utils
from shutil import rmtree
from pathlib import Path

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
def test_batch_tmalign_paral():
    target_dir = Path('data/tmalign_batch')
    try:
        tmalign = pdbanno.TMalign(num_workers=5)
        data_iter = [("data/gen.pdb", "data/native.pdb")] * 5
        res = tmalign.batch_align(data_iter, target_dir)
        # all tmscore should be 0.95422
        assert all(res['tmscore'] == 0.95422)
        # all rmsd should be 0.99
        assert all(res['rmsd'] == 0.99)
        # dimension of res should be 5
        assert res.shape[0] == 5

    finally:
        if target_dir.exists():
            rmtree(target_dir)


@pytest.mark.skipif(not tm_available, reason=helper_msg)
def test_batch_tmalign():
    target_dir = Path('data/tmalign_batch')
    try:
        tmalign = pdbanno.TMalign(num_workers=1)
        data_iter = [("data/gen.pdb", "data/native.pdb")] * 5
        res = tmalign.batch_align(data_iter, target_dir)
        # all tmscore should be 0.95422
        assert all(res['tmscore'] == 0.95422)
        # all rmsd should be 0.99
        assert all(res['rmsd'] == 0.99)
        # dimension of res should be 5
        assert res.shape[0] == 5

    finally:
        if target_dir.exists():
            rmtree(target_dir)


@pytest.mark.skipif(not tm_available, reason=helper_msg)
def test_batch_tmalign_async():
    import asyncio
    loop = asyncio.get_event_loop()
    target_dir = Path('data/tmalign_batch')
    try:
        tmalign = pdbanno.TMalign(num_workers=1)
        data_iter = [("data/gen.pdb", "data/native.pdb")] * 5
        res = loop.run_until_complete(tmalign.async_batch_align(data_iter, target_dir))
        # all tmscore should be 0.95422
        assert all(res['tmscore'] == 0.95422)
        # all rmsd should be 0.99
        assert all(res['rmsd'] == 0.99)
        # dimension of res should be 5
        assert res.shape[0] == 5

    finally:
        if target_dir.exists():
            rmtree(target_dir)
