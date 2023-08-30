import time
import functools

from pathlib import Path
from typing import Union

FilePath = Union[Path, str]


def ensure_path(path: FilePath) -> Path:
    if not isinstance(path, Path):
        path = Path(path).expanduser().resolve().absolute()
    return path


def max_retry(max=3, err_types=(RuntimeError,), interval=1):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            current_err = None
            for _ in range(max):
                try:
                    return func(*args, **kwargs)
                except err_types as e:
                    print(f'Error: {e}')
                    current_err = e
                    time.sleep(interval)
            raise RuntimeError(f'Failed after {max} retries') from current_err
        return wrapper
    return decorator