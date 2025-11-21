import pytest

from protools import aa


def test_validate_seq():
    aa.validate_seq('ACDEFGHIKLMNPQRSTVWY')
    with pytest.raises(ValueError):
        aa.validate_seq('ACDEFGHIKLMNPQRSTVWYX')
    aa.validate_seq('ACDEFGHIKLMNPQRSTVWY--', extra_symbols='-')


@pytest.mark.parametrize(
        'a1, a2, comparsion, is_equal',
        [
            ('A', 'A', 'type', True),
            ('A', 'C', 'type', False),
            ('D', 'E', 'property', True),
            ('R', 'H', 'property', True),
            ('D', 'C', 'property', False),
            ('-', 'A', 'type', False),
            ('-', '-', 'type', True),
            ('-', 'X', 'property', False),
            ('X', 'X', 'property', True),
            ('X', 'X', 'type', True),
            ('A', 'a', 'type', True),
            ('A', 'D', lambda x, y: True, True),
            ('A', 'A', lambda x, y: False, False)
        ]
)
def test_aa_equal(a1, a2, comparsion, is_equal):
    assert aa.aa_equal(a1, a2, comparsion) == is_equal