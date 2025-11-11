import pytest
from protools import utils
from . import tools


@pytest.mark.parametrize(
        'seq',
        [123, None, 1.23, list()]
)
def test_ensure_seq_string_type_err(seq):
    with pytest.raises(TypeError):
        utils.ensure_seq_string(seq)  # type: ignore

@pytest.mark.parametrize(
        'seq',
        ['ACDEFGHIKLMNPQRSTVWY']
)
def test_ensure_seq_string(seq):
    """
    Test the ensure_seq_string function.
    """
    seq1 = seq
    from Bio.Seq import Seq
    seq2 = Seq(seq1)
    from Bio.SeqRecord import SeqRecord
    seq3 = SeqRecord(Seq(seq1))

    assert utils.ensure_seq_string(seq1) == seq1, "Failed on string input"
    assert utils.ensure_seq_string(seq2) == seq1, "Failed on Seq input"
    assert utils.ensure_seq_string(seq3) == seq1, "Failed on SeqRecord input"


@pytest.mark.parametrize(
        'intervals_pattern, exp_slices',
        [
            ('1-5,7,9-10,12-', [slice(0, 5), slice(6, 7), slice(8, 10), slice(11, None)]),
            ('1-9', [slice(0, 9)]),
            ('1-9, 8', [slice(0, 9)]),
            ('-10', [slice(0,10)]),
            ('5-9,7-11', [slice(4, 11)]),
            ('-', [slice(0, None)])
        ]
)
def test_intervals_parse(intervals_pattern, exp_slices):
    """
    Test the parsing of intervals from a string.
    """
    intervals = utils.Intervals(intervals_pattern)
    slices = list(intervals)
    assert slices == exp_slices


@pytest.mark.parametrize(
        'intervals_pattern, contains, not_contains',
        [
            ('1-5,7,9-10', [1,2,3,4,5,7,9,10, slice(8,10)], [0,6,8, slice(9,11)]),
            ('-', [1, 2, 100000, slice(10)], [0]),
            ('10-', [10,10000], [1,9, slice(8,19)])
        ]
)
def test_intervals_contains(intervals_pattern, contains, not_contains):
    """
    Test the containment functionality of Intervals.
    """
    intervals = utils.Intervals(intervals_pattern)
    for v in contains:
        assert v in intervals
    for v in not_contains:
        assert v not in intervals


@pytest.mark.parametrize(
        'slices, intervals_pattern, indices_len',
        [
            ([slice(0, 5), slice(6, 7), slice(8, 10)], '1-5,7,9-10', None),
            ([slice(0, None)], '-10', 10),
            ([slice(0, 10)],'1-10', None),
            ([slice(-10, -5)], '1-5', 10)
        ]
)
def test_intervals_from_slices(slices, intervals_pattern, indices_len):
    """
    Test the creation of Intervals from slices.
    """
    intervals = utils.Intervals.from_slices(slices, indices_len=indices_len)
    expected_intervals = utils.Intervals(intervals_pattern)
    
    assert intervals == expected_intervals


def test_intervals_add():
    """
    Test the addition of two Intervals objects.
    """
    intervals1 = utils.Intervals('1-5,7,9-10')
    intervals2 = utils.Intervals('2-6,8,10-11')
    assert intervals1 + 1 == intervals2, "Expected intervals1 + 1 to equal intervals2"


@pytest.mark.parametrize(
        'raw, result',
        [
            ('1-5,7,9-10', '6,8,11-'),
            ('2-', '1'),
            ('-', '')
        ]
)
def test_intervals_invert(raw: str, result: str):
    """
    Test the invert functionality of Intervals.
    """
    intervals = utils.Intervals(raw)
    inverted_intervals = ~intervals
    expected_intervals = utils.Intervals(result)
    
    assert inverted_intervals == expected_intervals


@pytest.mark.parametrize(
        'raw, result, low, high',
        [
            ('1-5,7,9-10', '6,8,11-12', 1, 12),
            ('2-5,7,9-10', '1,6,8,11-12', 1, 12),
            ('-5', '6-10', 1, 10)
        ]
)
def test_intervals_invert_with_bounds(raw: str, result: str, low: int, high: int):
    """
    Test the invert functionality of Intervals with specified bounds.
    """
    intervals = utils.Intervals(raw)
    inverted_intervals = intervals.invert(low=low, high=high)
    expected_intervals = utils.Intervals(result)
    
    assert inverted_intervals == expected_intervals


@pytest.mark.parametrize(
        'raw1, raw2, result',
        [
            ('1-5,7,9-10', '1-5,7,9-10', '1-5,7,9-10'),
            ('1-5,7,9-10', '3-8,10-12', '3-5,7,10'),
            ('3-8,10-12', '1-5,7,9-10', '3-5,7,10'),
            ('2', '2-', '2')
        ]
)
def test_intervals_and(raw1, raw2, result):
    """
    Test the intersection functionality of Intervals.
    """
    intervals1 = utils.Intervals(raw1)
    intervals2 = utils.Intervals(raw2)
    intersected_intervals = intervals1 & intervals2
    expected_intervals = utils.Intervals(result)
    
    assert intersected_intervals == expected_intervals


@pytest.mark.parametrize(
        'raw1, raw2, result',
        [
            ('1-5,7,9-10', '1-5,7,9-10', '1-5,7,9-10'),
            ('1-5,7,9-10', '3-8,10-12', '1-8,9-12'),
            ('3-8,10-12', '1-5,7,9-10', '1-8,9-12'),
            ('1', '2-', '-')
        ]
)
def test_intervals_or(raw1, raw2, result):
    """
    Test the union functionality of Intervals.
    """
    intervals1 = utils.Intervals(raw1)
    intervals2 = utils.Intervals(raw2)
    union_intervals = intervals1 | intervals2
    expected_intervals = utils.Intervals(result)
    
    assert union_intervals == expected_intervals


@pytest.mark.parametrize(
        'test_content',
        [
            b"Test content for compression",
            b"hello world",
            b"bodaojdoa",
            "Test content for compression",
            "hello world",
            "bodaojdoa"
        ]
)
def test_auto_compression_read(tmp_path, test_content):
    """
    Test the extract_compression context manager.
    """
    import gzip
    gz_file = tmp_path / "test.gz"
    mode = 'b' if isinstance(test_content, bytes) else 't'
    with gzip.open(gz_file, 'w' + mode) as f:
        f.write(test_content)
    
    with utils.auto_compression(gz_file, 'r' + mode) as f:
        content = f.read()
        assert content == test_content


@pytest.mark.parametrize(
        'raw_file, mode',
        [
            ('data/4zqk.cif', 'b'),
            ('data/4zqk.cif', 't'),
            ('LICENSE', 't')
            ]
)
def test_auto_compression_write(tmp_path, raw_file, mode):
    with open(raw_file, 'r' + mode) as f:
        original_content = f.read()
    with utils.auto_compression(tmp_path / "output.gz", 'w' + mode) as f:
        f.write(original_content)

    with utils.auto_compression(tmp_path / 'output.gz', 'r' + mode) as f:
        content = f.read()
        assert original_content == content


@pytest.mark.parametrize(
        'comressed_file, exp_md5sum, mode',
        [
            ('data/4zqk.cif.gz', '5363bcc30d0fdedf1c5a9218e4732304', 'b'),
            ('data/4zqk.cif.gz', '5363bcc30d0fdedf1c5a9218e4732304', 't')
        ]
)
def test_auto_compression_read(tmp_path, comressed_file, exp_md5sum, mode):
    with utils.auto_compression(comressed_file, mode='r' + mode) as f:
        content = f.read()
    with open(tmp_path / 'output', 'w' + mode) as f:
        f.write(content)
    
    assert tools.md5_equal(tmp_path / 'output', exp_md5sum)


@pytest.mark.xfail(reason="Not implemented")
def test_ensure_path():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_ensure_fileio():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_max_retry():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_catch_error():
    # TODO
    assert False

class TestCmdWrapperBase():

    @pytest.mark.xfail(reason="Not implemented")
    def test_find_command(self):
        # TODO
        assert False

    @pytest.mark.xfail(reason="Not implemented")
    def test_call():
        # TODO
        assert False

    @pytest.mark.xfail(reason="Not implemented")
    def test_async_call():
        # TODO
        assert False


@pytest.mark.xfail(reason="Not implemented")
def test_require_package():
    # TODO
    assert False


@pytest.mark.xfail(reason="Not implemented")
def test_local_cwd():
    # TODO
    assert False


    
    
    




