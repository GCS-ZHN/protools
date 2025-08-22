from protools import utils


def test_intervals_parse():
    """
    Test the parsing of intervals from a string.
    """
    intervals = utils.Intervals('1-5,7,9-10,12-')
    slices = list(intervals)
    expected_slices = [
        slice(0, 5),  # 1-5
        slice(6, 7),  # 7
        slice(8, 10), # 9-10
        slice(11, None)  # 12-
    ]
    assert slices == expected_slices, f"Expected {expected_slices}, got {slices}"


def test_intervals_contains():
    """
    Test the containment functionality of Intervals.
    """
    intervals = utils.Intervals('1-5,7,9-10')
    assert 3 in intervals, "Expected 3 to be in intervals"
    assert 6 not in intervals, "Expected 6 to not be in intervals"
    assert 1 in intervals, "Expected 1 to be in intervals"
    assert 5 in intervals, "Expected 5 to be in intervals"
    assert 7 in intervals, "Expected 7 to be in intervals"


def test_intervals_from_slices():
    """
    Test the creation of Intervals from slices.
    """
    slices = [slice(0, 5), slice(6, 7), slice(8, 10)]
    intervals = utils.Intervals.from_slices(slices)
    expected_intervals = utils.Intervals('1-5,7,9-10')
    
    assert intervals == expected_intervals, f"Expected {expected_intervals}, got {intervals}"


def test_intervals_add():
    """
    Test the addition of two Intervals objects.
    """
    intervals1 = utils.Intervals('1-5,7,9-10')
    intervals2 = utils.Intervals('2-6,8,10-11')
    assert intervals1 + 1 == intervals2, "Expected intervals1 + 1 to equal intervals2"


def test_intervals_invert():
    """
    Test the invert functionality of Intervals.
    """
    intervals = utils.Intervals('1-5,7,9-10')
    inverted_intervals = ~intervals
    expected_intervals = utils.Intervals('6,8,11-')
    
    assert inverted_intervals == expected_intervals, f"Expected {expected_intervals}, got {inverted_intervals}"


def test_intervals_invert_with_bounds():
    """
    Test the invert functionality of Intervals with specified bounds.
    """
    intervals = utils.Intervals('1-5,7,9-10')
    inverted_intervals = intervals.invert(low=1, high=12)
    expected_intervals = utils.Intervals('6,8,11-12')
    
    assert inverted_intervals == expected_intervals, f"Expected {expected_intervals}, got {inverted_intervals}"


def test_intervals_and():
    """
    Test the intersection functionality of Intervals.
    """
    intervals1 = utils.Intervals('1-5,7,9-10')
    intervals2 = utils.Intervals('3-8,10-12')
    intersected_intervals = intervals1 & intervals2
    expected_intervals = utils.Intervals('3-5,7,10')
    
    assert intersected_intervals == expected_intervals, f"Expected {expected_intervals}, got {intersected_intervals}"


def test_intervals_or():
    """
    Test the union functionality of Intervals.
    """
    intervals1 = utils.Intervals('1-5,7,9-10')
    intervals2 = utils.Intervals('3-8,10-12')
    union_intervals = intervals1 | intervals2
    expected_intervals = utils.Intervals('1-8,9-12')
    
    assert union_intervals == expected_intervals, f"Expected {expected_intervals}, got {union_intervals}"