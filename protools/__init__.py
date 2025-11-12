"""
Protein tool collections
"""

from protools.seqio import Fasta
from protools.utils import Intervals

try:
    from ._version import get_versions # type: ignore
    __version__ = get_versions()['version']
    del get_versions
except ImportError:
    pass