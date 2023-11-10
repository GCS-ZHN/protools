import time
import os
import shutil
import subprocess
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


def catch_error(logger, err_types=(RuntimeError,)):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except err_types as e:
                logger.error(f'{e.__class__.__name__}: {e}')
        return wrapper
    return decorator


class CmdNotFoundError(RuntimeError):
    
    def __init__(self, cmd: str):
        super().__init__(f'Command {cmd} not found.')


class CmdWrapperBase(object):
    def __init__(self, cmd: str):
        self.cmd = self.find_command(cmd)

    def _check_mod(self, cmd) -> str:
        if not os.path.exists(cmd):
            raise FileNotFoundError(f'{cmd} not existed.')
        if not os.access(cmd, os.X_OK):
            raise PermissionError(f'{cmd} not executable.')
        return cmd

    def find_command(self, cmd) -> str:
        """
        Tries to find the full path to a command.

        Parameters
        ----------
        cmd : str
            The command to be found.

        Returns
        -------
        str
            The full path to the command.
        """
        try:
            if os.path.sep in cmd:
                return self._check_mod(cmd)
            
            if f'{cmd.upper()}_PATH' in os.environ:
                return self._check_mod(os.environ[f'{cmd.upper()}_PATH'])

        except (FileNotFoundError, PermissionError) as e:
            raise CmdNotFoundError(cmd) from e
        
        cmd_path = shutil.which(cmd)
        if cmd_path is None:
            raise CmdNotFoundError(cmd)
        return cmd_path
    
    def __call__(self, *args, **kwargs):
        cmds = [self.cmd]
        cmds.extend(map(str, args))
        for k, v in kwargs.items():
            if len(k) == 1:
                cmds.append(f'-{k}')
            else:
                cmds.append(f'--{k}')
            cmds.append(str(v))
        return subprocess.run(cmds)


def require_package(package_name: str, install_cmd: str = None):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                __import__(package_name)
            except ImportError as e:
                msg = f'Package {package_name} is required.'
                if install_cmd:
                    msg += f" You can install it by `{install_cmd}`."
                raise RuntimeError(msg) from e
            return func(*args, **kwargs)
        return wrapper
    return decorator
