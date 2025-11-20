import pytest

from protools import aa


def test_validate_seq():
    aa.validate_seq('ACDEFGHIKLMNPQRSTVWY')
    with pytest.raises(ValueError):
        aa.validate_seq('ACDEFGHIKLMNPQRSTVWYX')
    aa.validate_seq('ACDEFGHIKLMNPQRSTVWY--', extra_symbols='-')
