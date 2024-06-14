import time
import os
import shutil
import subprocess
import asyncio
import functools
import logging

from pathlib import Path
from typing import Optional
from io import IOBase
from .typedef import FilePathType, FilePathOrIOType
from collections import namedtuple
from contextlib import contextmanager


AsyncCompletedProcess = namedtuple('AsyncCompletedProcess', ['stdout', 'stderr'])


def ensure_path(path: FilePathType) -> Path:
    """
    Convert the input to a Path object.

    Parameters
    ----------
    path : FilePathType
        The input path.

    Returns
    ----------
    Path
        The Path object with absolute path.
    """
    if not isinstance(path, Path):
        path = Path(path)
    return path.expanduser().resolve().absolute()


def ensure_fileio(path_or_io: FilePathOrIOType, mode: str = 'r') -> IOBase:
    """
    Ensure the input is a IO object. If the input is a path,
    open the path as a IO object.

    Parameters
    ----------
    path_or_io : FilePathOrIOType
        The input path or FileIO object.

    Returns
    ----------
    IOBase
        The IO object.
    bool
        Whether the IO object is created by this function.

    Notes
    ----------
    If the input is a IO object, it will be returned directly.
    You should close the IO object by yourself.

    Raises
    ----------
    ValueError
        If the input is a closed IO object or the mode is not matched.
    """
    if isinstance(path_or_io, (str, Path)):
        return ensure_path(path_or_io).open(mode), True
    elif isinstance(path_or_io, IOBase):
        if path_or_io.closed:
            raise ValueError('The input handler is closed.')
        if path_or_io.mode != mode:
            raise ValueError(f'The input handler is not in "{mode}" mode.')
        return path_or_io, False
    else:
        raise TypeError(f'Unsupported type: {type(path_or_io)}')


def max_retry(max: int = 3, err_types: Exception =(RuntimeError,), interval: int = 1):
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


def catch_error(logger: logging.Logger, err_types: Exception=(RuntimeError,)):
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
    """
    Raised when the command is not found.
    """
    def __init__(self, cmd: str, install_help: Optional[str] = None):
        message = f'Command {cmd} not found.'
        if install_help:
            message += ' ' + install_help.strip()
        super().__init__(message)


class CmdWrapperBase(object):
    """
    A base class for wrapping a command.

    Parameters
    ----------
    cmd : str
        The command to be wrapped.
    short_mode : bool, optional
        Whether the command is in short mode.
        If it is in short mode, the command will be
        formatted as `-key v` instead of `--key v`.
    num_workers : int, optional
        The number of workers for asynchronous calls.
    install_help : str, optional
        The command to install the command.
    bool2flag : bool, optional
        Whether the boolean value will be converted to flag.
        If it is True, the boolean argument will be converted
        to a flag without value.
    """
    def __init__(self, 
                 cmd: str,
                 short_mode: bool = False,
                 num_workers: int = 1,
                 install_help: Optional[str] = None,
                 list2nargs: bool = True,
                 bool2flag: bool = False):
        self.cmd = self.find_command(cmd, install_help)
        self.short_mode = short_mode
        self.semaphore = asyncio.Semaphore(num_workers)
        self.bool2flag = bool2flag
        self.list2nargs = list2nargs

    @classmethod
    def _check_mod(cls, cmd) -> str:
        if not os.path.exists(cmd):
            raise FileNotFoundError(f'{cmd} not existed.')
        if not os.access(cmd, os.X_OK):
            raise PermissionError(f'{cmd} not executable.')
        return cmd

    @classmethod
    def find_command(cls, cmd: str, install_help: Optional[str] = None) -> str:
        """
        Tries to find the full path to a command.

        Parameters
        ----------
        cmd : str
            The command to be found.
        install_help : str, optional
            The command to install the command, it will be shown
            when the command is not installed if provided.

        Returns
        -------
        str
            The full path to the command.
        """
        try:
            if os.path.sep in cmd:
                return cls._check_mod(cmd)
            
            if f'{cmd.upper()}_PATH' in os.environ:
                return cls._check_mod(os.environ[f'{cmd.upper()}_PATH'])

        except (FileNotFoundError, PermissionError) as e:
            raise CmdNotFoundError(cmd, install_help) from e
        
        cmd_path = shutil.which(cmd)
        if cmd_path is None:
            raise CmdNotFoundError(cmd, install_help)
        return cmd_path
    
    def _format_args(self, *args, **kwargs) -> list:
        """
        Format the arguments.
        """
        cmds = [self.cmd]
        cmds.extend(map(str, args))
        for k, v in kwargs.items():
            if len(k) == 1 or self.short_mode:
                cmds.append(f'-{k}')
            else:
                cmds.append(f'--{k}')
            if self.bool2flag and isinstance(v, bool):
                if not v:
                    cmds.pop()
                continue
            if self.list2nargs and isinstance(v, (list, tuple)):
                cmds.extend(map(str, v))
                continue
            cmds.append(str(v))
        return cmds
    
    def __call__(
            self, 
            *args, 
            stdout=None,
            stderr=None,
            stdin=None,
            cwd=None,
            env=None,
            **kwargs):
        """
        Run the command.

        Parameters
        ----------
        *args
            The positional arguments.
        stdout
            The stdout of the command. See subprocess.run for more details.
        stderr
            The stderr of the command. See subprocess.run for more details.
        stdin
            The stdin of the command. See subprocess.run for more details.
        cwd
            The current working directory of the command.
        env
            The environment variables of the command.
        **kwargs
            The keyword arguments.
        """
        cmds = self._format_args(*args, **kwargs)
        process = subprocess.run(
            cmds, 
            stdout=stdout, 
            stderr=stderr, 
            stdin=stdin,
            cwd=cwd,
            env=env,
            check=True)
        return process

    async def async_call(
            self, 
            *args,
            stdout=None,
            stderr=None,
            stdin=None,
            **kwargs) -> AsyncCompletedProcess:
        """
        Run the command asynchronously.
        As each command as an subprocess,
        It will run parallelly if multiple commands are called.

        Parameters
        ----------
        *args
            The positional arguments.
        stdout
            The stdout of the command. See subprocess.run for more details.
        stderr
            The stderr of the command. See subprocess.run for more details.
        stdin
            The stdin of the command. See subprocess.run for more details.
        **kwargs
            The keyword arguments.

        Returns
        ---------
        AsyncCompletedProcess
            The result of the command.
        """
        async with self.semaphore:
            cmds = self._format_args(*args, **kwargs)
            process = await asyncio.create_subprocess_exec(
                *cmds, 
                stdout=stdout, 
                stderr=stderr, 
                stdin=stdin)
            async_stdout, async_stderr = await process.communicate()
            return AsyncCompletedProcess(async_stdout, async_stderr)
    
    def __getattr__(self, name: str):
        """
        Create a subcommand wrapper.
        This is useful when the command has subcommands.

        Parameters
        ----------
        name : str
            The name of the subcommand.
            If the name starts with `async_`, the subcommand
            will be wrapped as an asynchronous function.
            And really called subcommand name will be the
            name without `async_`.

        Returns
        ----------
        function or coroutine function
            The wrapper of the subcommand.
        """
        if name.startswith('__'):
            # avoid conflicts with built-in magic methods if it not implemented
            # otherwise, it will cause unknown behavior
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute '{name}'")
        if name.startswith('async_'):
            name = name[6:]
            async def wrapper(*args, **kwargs):
                return await self.async_call(name, *args, **kwargs)
        else:
            def wrapper(*args, **kwargs):
                return self(name, *args, **kwargs)
        return wrapper


def require_package(package_name: str, install_cmd: Optional[str] = None):
    """
    A decorator to check if a package is installed.
    Some functions may require a package to be installed.

    Parameters
    ----------
    package_name : str
        The name of the package.
    install_cmd : str, optional
        The command to install the package, it will be shown
        when the package is not installed if provided.
    """
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


@contextmanager
def local_cwd(path: FilePathType):
    """
    Change the current working directory locally.

    Parameters
    ----------
    path : FilePathType
        The path to the new working directory.
    """
    cwd = Path.cwd()
    path = ensure_path(path)
    if not path.is_dir():
        raise NotADirectoryError(f'{path} is not a directory.')
    if cwd != path:
        os.chdir(path)
    yield path
    if cwd != path:
        os.chdir(cwd)
