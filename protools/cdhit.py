import re
import tempfile
import pandas as pd

from pathlib import Path
from typing import Iterable
from .utils import CmdWrapperBase
from .seqio import read_fasta, save_fasta, temp_fasta

class CdHit(CmdWrapperBase):
    
    def __init__(self, cmd: str = 'cd-hit'):
        super().__init__(cmd)
    
    def __call__(self, 
            input_file: Path,
            output_file: Path,
            cutoff: float = 0.9,
            num_threads: int = 1,
            word_length: int = 5,
            memory: int = 800,
            *args, **kwargs):
        return super().__call__(
            i=input_file,
            o=output_file,
            c=cutoff,
            n=word_length,
            M=memory,
            T=num_threads,
            *args, **kwargs)
    
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
    def parse_cluster(cluster_file: Path) -> pd.DataFrame:
        return pd.DataFrame(CdHit._iter_cluster(cluster_file))


def unique_fasta(fasta_file: Path, output_file: Path, *args, **kwargs):
    cdhit = CdHit()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    temp_output = tempfile.NamedTemporaryFile(
        prefix='cdhit-',
        dir=output_file.parent)
    temp_input, id_map = temp_fasta(fasta_file)
    cdhit(temp_input.name, temp_output.name, *args, **kwargs)
    cluster = cdhit.parse_cluster(temp_output.name + '.clstr')
    cluster['sequence_id'] = cluster['sequence_id'].apply(lambda x: id_map[x])
    cluster.to_csv(output_file.with_suffix('.cluster.csv'), index=False)
    cluster = cluster[cluster['representative']]
    fasta = read_fasta(fasta_file)
    uniqued = (fasta[seq_id] for seq_id in cluster['sequence_id'])
    save_fasta(uniqued, output_file)
    Path(temp_output.name + '.clstr').unlink()


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    unique_parser = subparsers.add_parser('unique')
    unique_parser.add_argument('--fasta_file', '-i', type=Path, required=True)
    unique_parser.add_argument('--output_file', '-o', type=Path, required=True)
    unique_parser.add_argument('--cutoff', '-c', type=float, default=0.9)
    unique_parser.add_argument('--num_threads', '-t', type=int, default=1)

    args = parser.parse_args()

    unique_fasta(
        args.fasta_file,
        args.output_file,
        cutoff=args.cutoff,
        num_threads=args.num_threads)
