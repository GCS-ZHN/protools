
from .seqio import read_fasta, Fasta
from Bio.SeqUtils import IUPACData
from Bio.SeqRecord import SeqRecord
from typing import Tuple
from pathlib import Path

_REGISTED_FILTERS = {}

def filter_register(name: str):
    def decorator(fn):
        _REGISTED_FILTERS[name] = fn
        return fn
    return decorator


@filter_register('standard_aa')
def standard_aa_filter(items: Tuple[str, SeqRecord]) -> bool:
    """Filter out non-standard amino acids."""
    _, record = items
    return all(aa in IUPACData.protein_letters for aa in record.seq)


def filter_fasta(fasta_file: str, *filter_name: str) -> Fasta:
    """Filter a FASTA file using a filter function."""
    fasta = read_fasta(fasta_file)
    return Fasta(filter(lambda x: all(_REGISTED_FILTERS[fn](x) for fn in filter_name), fasta.items()))


if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Filter a FASTA file.')
    parser.add_argument('--input', '-i', required=True, type=Path, help='FASTA file to filter')
    parser.add_argument('--output', '-o', type=Path, help='Output file if not provided will print to stdout')
    parser.add_argument('--filters', '-f', nargs='+', help='Filters to apply')
    args = parser.parse_args()

    fasta = filter_fasta(args.input, *args.filters)
    if args.output:
        fasta.to_fasta(args.output)
    else:
        fasta.to_fasta(sys.stdout)
