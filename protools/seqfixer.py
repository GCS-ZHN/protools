from .seqio import read_fasta, save_fasta
from .utils import FilePath, ensure_path
from pathlib import Path


def deduplicated_seq(data: FilePath):
    """
    Remove duplicated sequence by remaining the
    first record.
    """
    fasta = read_fasta(data)
    seqs = set()
    def _filter(x):
        if x.seq not in seqs:
            seqs.add(x.seq)
            return True
        return False
        
    save_fasta(
            filter(_filter, fasta.values()),
            data.with_name(f"{data.stem}_deduplicated.fasta"))

    print(f"Total seqs: {len(fasta)}, unique seqs: {len(seqs)}")


def split(data: FilePath, n: int):
    """
    Split sequences into `n` fasta file.
    """
    data = ensure_path(data)
    fasta = read_fasta(data)
    out_dir = data.with_name(data.stem)
    out_dir.mkdir()
    for i in range(n):
        out_fasta = out_dir / f"seq_{i}.fasta"
        seqs = (res[1] for res in enumerate(fasta.values()) if res[0] % n == i)
        save_fasta(seqs, out_fasta)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd')
    rm_parser = subparsers.add_parser("deduplicated")
    rm_parser.add_argument('--data', '-i', type=Path, required=True)

    split_parser = subparsers.add_parser("split")
    split_parser.add_argument('--data', '-i', type=Path, required=True)
    split_parser.add_argument('--num', '-n', type=int, required=True)
    args = parser.parse_args()
    if args.cmd == 'deduplicated':
        deduplicated_seq(
            args.data
        )
    if args.cmd == 'split':
        split(
            args.data,
            args.num
        )