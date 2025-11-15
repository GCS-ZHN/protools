import hashlib

from pathlib import Path
from typing import Iterable


def md5sum(filename: Path) -> str:
    """Compute the md5sum of a file."""
    with open(filename, 'rb') as f:
        md5 = hashlib.md5()
        while True:
            data = f.read(65536)
            if not data:
                break
            md5.update(data)
        return md5.hexdigest()


def md5_equal(file: Path, md5: str) -> bool:
    """Check if the md5sum of a file matches the given md5sum."""
    return md5sum(file).startswith(md5)


def md5_equal_file(file1: Path, file2: Path) -> bool:
    """Check if the md5sums of two files are equal."""
    return md5sum(file1) == md5sum(file2)


def all_equal(data: Iterable) -> bool:
    prev = None
    for i, v in enumerate(data):
        if i == 0:
            prev = v
        if prev != v:
            return False
    return True
