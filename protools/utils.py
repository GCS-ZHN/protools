from pathlib import Path
from typing import Union

FilePath = Union[Path, str]


def ensure_path(path: FilePath) -> Path:
    if not isinstance(path, Path):
        path = Path(path).expanduser().resolve().absolute()
    return path