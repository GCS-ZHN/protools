import hashlib

from pathlib import Path


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