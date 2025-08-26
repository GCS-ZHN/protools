import time
import os
import shutil
import subprocess
import asyncio
import functools
import logging
import warnings

from pathlib import Path
from typing import Optional, List, Callable
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


class Intervals(object):
    """
    A class to represent intervals, such as "1-10, 20-30, 40-50".
    It can be used to represent the slices of a sequence. Interval
    will be merged if they are overlapped.


    Parameters
    ----------
    intervals_pattern : str
        The pattern of the intervals, such as "1-10, 20-30, 40-50".
        Negative values are not supported. If start is emtpy, it
        will be regarded as first position. If stop is empty, it
        will be regarded as infinite value. Each interval is separated
        by comma.
    zero_based : bool, optional
        Whether the interval is zero-based. Default is False.
        First position is zero if set it to True, otherwise,
        first position will be one.
    end_inclusive : bool, optional
        Whether the end of the interval is inclusive. Default is True.
        If it is True, the end of the interval will be included in the
        interval, otherwise, it will be excluded.
    interval_hook : Callable[[str], str], optional
        A hook to process each interval pattern. It will be called
        before parsing each interval pattern. It is useful when the
        interval pattern is not standard. For example, the interval
        pattern is "1-10(strong), 20-30(weakly), 40-50" The hook can be
        used to remove the content in the bracket.

    Examples
    ----------
    >>> Intervals('1-10, 20-30, 40-50')
    '1-10,20-30,40-50' (zero_based=False, end_inclusive=True, dynamic=False)
    >>> Intervals.from_slices([slice(5, None), slice(25, 35)])
    '6-' (zero_based=False, end_inclusive=True, dynamic=True)
    """

    def __init__(self,
                 intervals_pattern: str,
                 zero_based: bool = False,
                 end_inclusive: bool = True,
                 interval_hook: Optional[Callable[[str], str]] = None):
        self.dynamic = False
        self.zero_based = zero_based
        self.end_inclusive = end_inclusive
        self._interval_slices = self._parse_intervals(
            intervals_pattern, interval_hook)

    @classmethod
    def from_slices(cls,
                    slices: List[slice],
                    zero_based: bool = False,
                    end_inclusive: bool = True):
        """
        Create Intervals from slice list.
        """
        interval = cls('', 
                       zero_based=zero_based,
                       end_inclusive=end_inclusive)
        slices = [slice(s.start, s.stop) if s.start else slice(0, s.stop) for s in slices]
        interval._interval_slices = cls._merge_intervals(slices)
        interval.dynamic = any([i.stop is None for i in interval._interval_slices])
        return interval

    @staticmethod
    def _merge_intervals(interval_slices: List[slice]) -> List[slice]:
        interval_slices.sort(key=lambda x: x.start)
        merged_interval_slices = []
        for interval in interval_slices:
            if not merged_interval_slices:
                merged_interval_slices.append(interval)
            else:
                last_interval = merged_interval_slices[-1]
                if  last_interval.stop is None or interval.start <= last_interval.stop:
                    if last_interval.stop is None or interval.stop is None:
                        stop = None
                    else:
                        stop = max(last_interval.stop, interval.stop)
                    merged_interval_slices[-1] = slice(last_interval.start, stop)
                else:
                    merged_interval_slices.append(interval)
        return merged_interval_slices

    def _convert_interval(self, start: int, stop: int) -> tuple:
        if not self.zero_based:
            start = 0 if start is None else start - 1
            if stop is not None:
                stop -= 1
        if self.end_inclusive:
            if stop is not None:
                stop += 1
        assert start >= 0, "Invalid interval: negative value is not supported"
        if stop is not None:
            assert stop > start, "Invalid interval: stop should not less than start"
        return start, stop

    def _parse_intervals(self, 
                       intervals_patttern: str,
                       interval_hook: Optional[Callable[[str], str]] = None) -> List[slice]:
        interval_slices = []
        for i in intervals_patttern.split(','):
            if interval_hook:
                i = interval_hook(i)
            i = i.strip()
            if not i:
                continue
            items = i.split('-')
            items = [int(x) if x else None for x in items]
            if len(items) == 1:
                start, stop = items[0], items[0]
            elif len(items) == 2:
                start, stop = items
            else:
                raise ValueError(f"Invalid intervals: {intervals_patttern}")
            start, stop = self._convert_interval(start, stop)
            if stop is None:
                self.dynamic = True
            interval_slices.append(slice(start, stop))
        return self._merge_intervals(interval_slices)

    def size(self) -> int:
        """
        Return intervals' total size.
        """
        if self.dynamic:
            raise ValueError("Dynamic interval does not have length")
        return sum([i.stop - i.start for i in self._interval_slices])

    def create_view_symbol(self, seq: str, mark_symbol: str = '*', other_symbol: str = ' ') -> str:
        """
        Create a view symbol for the intervals in a sequence.

        Parameters
        ----------
        seq : str
            The sequence to create the view symbol.
        mark_symbol : str, optional
            The symbol to mark the intervals. Default is '*'.
        other_symbol : str, optional
            The symbol to mark the other positions. Default is ' '.

        Returns
        -------
        str
            The view symbol of the intervals in the sequence.
        """
        seq_len = len(seq)
        marks = [other_symbol] * seq_len
        for interval in self._interval_slices:
            stop = seq_len if interval.stop is None else min(seq_len, interval.stop)
            for i in range(interval.start, stop):
                marks[i] = mark_symbol
        return ''.join(marks)

    def intersect(self, intervals: 'Intervals') -> 'Intervals':
        """
        Calculate the intersection of two intervals.
        """
        new_interval_slices = []
        dynamic = False
        for i1 in self._interval_slices:
            for i2 in intervals._interval_slices:
                start = max(i1.start, i2.start)
                if i1.stop is None:
                    stop = i2.stop
                elif i2.stop is None:
                    stop = i1.stop
                else:
                    stop = min(i1.stop, i2.stop)
                if stop is None or start < stop:
                    new_interval_slices.append(slice(start, stop))
                if stop is None:
                    dynamic = True
        res = self.__class__(
            '',
            zero_based=self.zero_based,
            end_inclusive=self.end_inclusive)
        res._interval_slices = new_interval_slices
        res.dynamic = dynamic
        return res

    def invert(self, low: int = None, high: int = None) -> 'Intervals':
        """
        Invert the intervals.

        Parameters
        ----------
        low : int, optional
            The lower bound of the interval. should not less than
            first position. If it is None, it will be regarded as
            the first position (0 or 1 depends on `zero_based`).
        high : int, optional
            The upper bound of the interval.

        """
        cls = self.__class__
        new_interval = cls('', 
            zero_based=self.zero_based,
            end_inclusive=self.end_inclusive)
        intervals = []
        start_pos = 0
        for i in self._interval_slices:
            if start_pos is None:
                break
            if i.start > start_pos:
                intervals.append(slice(start_pos, i.start))
            start_pos = i.stop

        if start_pos is not None:
            intervals.append(slice(start_pos, None))

        new_interval._interval_slices = intervals
        low, high = self._convert_interval(low, high)
        bound_interval = cls.from_slices(
            [slice(low, high)],
            zero_based=self.zero_based,
            end_inclusive=self.end_inclusive)
        new_interval = new_interval.intersect(bound_interval)
        return new_interval

    def union(self, other: 'Intervals') -> 'Intervals':
        """
        Calculate the union of two intervals.

        Parameters
        ----------
        other : Intervals
            The other intervals to be united with.
        """
        return self.__class__.from_slices(
            self._interval_slices + other._interval_slices,
            zero_based=self.zero_based,
            end_inclusive=self.end_inclusive)

    def __str__(self):
        patterns = []
        for interval in self._interval_slices:
            start, stop = interval.start, interval.stop
            if not self.zero_based:
                start += 1
                if stop is not None:
                    stop += 1
            if self.end_inclusive and stop is not None:
                stop -= 1
            start = str(start)
            stop = '' if stop is None else str(stop)
            if  start == stop:
                patterns.append(str(start))
            else:
                patterns.append(f"{start}-{stop}")
        return ','.join(patterns)

    def __repr__(self):
        return f"'{self}' (zero_based={self.zero_based}, end_inclusive={self.end_inclusive}, dynamic={self.dynamic})"

    def __eq__(self, other: 'Intervals') -> bool:
        return self._interval_slices == other._interval_slices

    def __add__(self, shift: int) -> 'Intervals':
        assert isinstance(shift, int), "Shift must be an integer"
        new_intervals = self.__class__('',
                                zero_based=self.zero_based,
                                end_inclusive=self.end_inclusive)
        new_interval_slices = []
        dynamic = False
        for i in self._interval_slices:
            new_start = i.start + shift
            if new_start < 0:
                new_start = 0
                warnings.warn(
                    f"Shift {shift} is too large for interval {self}, start is set to 0")
            if i.stop is None:
                new_stop = None
                dynamic = True
            else:
                new_stop = i.stop + shift
                if new_stop < 0:
                    new_stop = 0
                    warnings.warn(
                        f"Shift {shift} is too large for interval {self}, stop is set to 0")
            if new_stop is None or new_start < new_stop:
                new_interval_slices.append(slice(new_start, new_stop))
        new_intervals._interval_slices = new_interval_slices
        new_intervals.dynamic = dynamic
        return new_intervals

    def __sub__(self, shift: int) -> 'Intervals':
        return self.__add__(-shift)
    
    def __iter__(self):
        return iter(self._interval_slices)

    def __contains__(self, value):
        """
        Check if the value is in the intervals.
        """
        if isinstance(value, int):
            if not self.zero_based:
                value -= 1
            for i in self._interval_slices:
                if i.start <= value and (i.stop is None or value < i.stop):
                    return True
            return False
        if isinstance(value, self.__class__):
            cross_intervals = self.intersect(value)
            return cross_intervals == value
        if isinstance(value, slice):
            value = self.__class__.from_slices([value])
            return self.__contains__(value)
        raise TypeError('Element value must be int, slice or Intervals')

    def __invert__(self) -> 'Intervals':
        """
        Default behavior of `~` operator is to invert the intervals.
        """
        return self.invert()

    __and__ = intersect
    __or__ = union
