import os
import re
import shutil
import subprocess
from typing import Iterable
import pandas as pd

class CmdNotFoundError(RuntimeError):
    
    def __init__(self, cmd: str):
        super().__init__(f'Command {cmd} not found.')


class CmdWrapperBase(object):

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


class CdHit(CmdWrapperBase):
    
    def __init__(self, cmd: str = 'cd-hit'):
        self.cmd = self.find_command(cmd)
    
    def __call__(self, 
            input_file: str,
            output_file: str,
            cutoff: float = 0.9,
            num_threads: int = 1,
            word_length: int = 5,
            memory: int = 800,
            *args, **kwargs):
        cmds = [
            self.cmd,
            '-i', input_file,
            '-o', output_file,
            '-c', cutoff,
            '-n', word_length,
            '-M', memory,
            '-T', num_threads
        ]
        cmds.extend(args)
        return subprocess.run(cmds, **kwargs)
    
    @staticmethod
    def _iter_cluster(cluster_file: str) -> Iterable[dict]:
        common_pattern = '\\d+\t(\\d+)aa, >(.+)\\.\\.\\.'
        repres_pattern = re.compile(common_pattern + ' \\*')
        members_pattern = re.compile(common_pattern + ' at ([\\d\\.]+)%')
        with open(cluster_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    cluster_id = int(line.split()[1])
                else:
                    match = repres_pattern.match(line)
                    if match:
                        yield {
                            'cluster_id': cluster_id,
                            'sequence_id': match.group(2),
                            'sequence_length': match.group(1),
                            'representative': True,
                            'similarity': 100.0
                        }
                    else:
                        match = members_pattern.match(line)
                        if match:
                            yield {
                                'cluster_id': cluster_id,
                                'sequence_id': match.group(2),
                                'sequence_length': match.group(1),
                                'representative': False,
                                'similarity': float(match.group(3))
                            }
                        else:
                            raise ValueError(f'Unrecognized line: {line}')

    @staticmethod
    def parse_cluster(cluster_file: str) -> pd.DataFrame:
        return pd.DataFrame(CdHit._iter_cluster(cluster_file))
