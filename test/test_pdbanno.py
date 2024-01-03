from protools import pdbanno
from shutil import rmtree
from pathlib import Path


def test_tmalign():
    tmalign = pdbanno.TMalign()
    v = tmalign("data/gen.pdb", "data/native.pdb", 'data', 'aligned_gen')
    assert v == (0.95422, 0.99)


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
