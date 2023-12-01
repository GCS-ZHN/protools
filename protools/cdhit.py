import re
import tempfile
from pathlib import Path
from typing import Iterable

import pandas as pd

from .seqio import read_fasta, save_fasta, temp_fasta
from .utils import CmdWrapperBase


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


def parse_cluster(cluster_file: Path) -> pd.DataFrame:
    return pd.DataFrame(_iter_cluster(cluster_file))

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


class CdHit2D(CmdWrapperBase):
    
    def __init__(self, cmd: str = 'cd-hit-2d'):
        super().__init__(cmd)

    def __call__(self,
                 source_file: Path,
                target_file: Path,
                output_file: Path,
                cutoff: float = 0.9,
                num_threads: int = 1,
                word_length: int = 5,
                memory: int = 800,
                *args, **kwargs):
            return super().__call__(
                "-i", source_file,
                "-i2", target_file,
                "-o", output_file,
                "-c", cutoff,
                "-n", word_length,
                "-M", memory,
                "-T", num_threads,
                *args, **kwargs)

def unique_fasta(fasta_file: Path, output_file: Path, *args, **kwargs):
    cdhit = CdHit()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    temp_output = tempfile.NamedTemporaryFile(
        prefix='cdhit-',
        dir=output_file.parent)
    temp_input, id_map = temp_fasta(fasta_file)
    cdhit(temp_input.name, temp_output.name, *args, **kwargs)
    cluster = parse_cluster(temp_output.name + '.clstr')
    cluster['sequence_id'] = cluster['sequence_id'].apply(lambda x: id_map[x])
    cluster.to_csv(output_file.with_suffix('.cluster.csv'), index=False)
    cluster = cluster[cluster['representative']]
    fasta = read_fasta(fasta_file)
    uniqued = (fasta[seq_id] for seq_id in cluster['sequence_id'])
    save_fasta(uniqued, output_file)
    Path(temp_output.name + '.clstr').unlink()


def unique_fasta_2d(source_file: Path, target_file: Path, output_file: Path, *args, **kwargs):
    cdhit = CdHit2D()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    temp_output = tempfile.NamedTemporaryFile(
        prefix='cdhit-',
        dir=output_file.parent)
    temp_source, source_id_map = temp_fasta(source_file, id_prefix='s_')
    temp_target, target_id_map = temp_fasta(target_file, id_prefix='t_')
    cdhit(temp_source.name, temp_target.name, temp_output.name, *args, **kwargs)
    cluster = parse_cluster(temp_output.name + '.clstr')
    cluster['sequence_id'] = cluster.apply(
        lambda x: source_id_map[x['sequence_id']] if x['representative'] else target_id_map[x['sequence_id']],
        axis=1)
    cluster.to_csv(output_file.with_suffix('.cluster.csv'), index=False)
    fasta = read_fasta(source_file)
    group = cluster.groupby('cluster_id')
    
    # only one member each cluster
    uniqued = (fasta[g['sequence_id'].values[0]] for _, g in group if len(g) == 1)

    save_fasta(uniqued, output_file)
    Path(temp_output.name + '.clstr').unlink()


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    common_parser = ArgumentParser(add_help=False)
    common_parser.add_argument('--input_file', '-i', type=Path, required=True)
    common_parser.add_argument('--output_file', '-o', type=Path, required=True)
    common_parser.add_argument('--cutoff', '-c', type=float, default=0.9)
    common_parser.add_argument('--num_threads', '-t', type=int, default=1)
    common_parser.add_argument('--word_length', '-n', type=int, default=5)
    subparsers = parser.add_subparsers(dest='cmd')
    unique_parser = subparsers.add_parser('unique', parents=[common_parser])
    unique2d_parser = subparsers.add_parser('unique2d', parents=[common_parser])
    unique2d_parser.add_argument('--target_file', '-i2', type=Path, required=True)

    args = parser.parse_args()

    if args.cmd == 'unique':
        unique_fasta(
            args.input_file,
            args.output_file,
            cutoff=args.cutoff,
            num_threads=args.num_threads,
            word_length=args.word_length)
        
    elif args.cmd == 'unique2d':
        unique_fasta_2d(
            args.input_file,
            args.target_file,
            args.output_file,
            cutoff=args.cutoff,
            num_threads=args.num_threads,
            word_length=args.word_length)
    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")
