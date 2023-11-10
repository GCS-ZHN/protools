import pandas as pd

from .utils import require_package
from .seqio import read_fasta, df2fasta
from pathlib import Path
try:
    require_abnumber = require_package("abnumber",  "conda install -c bioconda abnumber")
    from abnumber import Chain, ChainParseError
except ImportError:
    pass


@require_abnumber
def remove_constant_region(fasta_file: Path, strict: bool = False, ignore_error: bool = False):
    fasta = read_fasta(fasta_file)

    for seq_id, seq in fasta.items():
        seq = str(seq.seq)

        # below 135 is a heuristic value
        if not strict and len(seq) <= 135:
            yield {"id": seq_id, "seq": seq}
            continue

        try:
            chain = Chain(seq, scheme='imgt')
            yield {"id": seq_id, "seq": str(chain)}
        except ChainParseError as e:
            if ignore_error:
                yield {"id": seq_id, "seq": seq}
            else:
                raise RuntimeError(f"Error parsing {seq_id}: {e}") from e


@require_abnumber
def annotate_chain_type(fasta_file: Path):
    fasta = read_fasta(fasta_file)
    for seq_id, seq in fasta.items():
        seq = str(seq.seq)
        try:
            chain = Chain(seq.seq, scheme='imgt')
            yield {"id": seq_id, "chain_type": chain.chain_type}
        except ChainParseError:
            yield {"id": seq_id, "chain_type": "unknown"}


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")

    remove_constant_region_parser = subparsers.add_parser("rm_const")
    remove_constant_region_parser.add_argument("--fasta", "-i", type=Path, required=True)
    remove_constant_region_parser.add_argument("--strict", "-s", action="store_true")
    remove_constant_region_parser.add_argument("--ignore_error", "-e", action="store_true")
    remove_constant_region_parser.add_argument("--output", "-o", type=Path, required=True)

    annotate_chain_type_parser = subparsers.add_parser("anno_chain")
    annotate_chain_type_parser.add_argument("--fasta", "-i", type=Path, required=True)
    annotate_chain_type_parser.add_argument("--output", "-o", type=Path, required=True)

    args = parser.parse_args()

    if args.subcommand == "rm_const":
        df = pd.DataFrame(remove_constant_region(args.fasta, args.strict, args.ignore_error))
        df2fasta(df, args.output, 'seq', id_col='id')
    elif args.subcommand == "anno_chain":
        df = pd.DataFrame(annotate_chain_type(args.fasta))
        df.to_csv(args.output, index=False)
    else:
        raise RuntimeError(f"Unknown subcommand {args.subcommand}")