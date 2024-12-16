import pandas as pd

from .utils import require_package
from .seqio import read_fasta, df2fasta
from pathlib import Path
try:
    require_anarci = require_package('anarci', 'conda install -c bioconda anarci')
    require_abnumber = require_package("abnumber",  "conda install -c bioconda abnumber")
    from anarci import run_anarci
    from abnumber import Chain, ChainParseError
except ImportError:
    pass

from antpack import SingleChainAnnotator
_HC_ANNOTATOR = SingleChainAnnotator(chains=['H'], scheme='imgt')
_LC_ANNOTATOR = SingleChainAnnotator(chains=['K', 'L'], scheme='imgt')


IMGT_BORDERS = [27, 39, 56, 66, 105, 118, 129]
REGION_NAMES = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']


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


def numbering_seq(seq: str, chain: str) -> tuple:
    if chain == 'H':
        numbering, percent_identity, chain_type, err_message = _HC_ANNOTATOR.analyze_seq(seq)
        assert chain_type == 'H', f"Chain type {chain_type} != H"
    else:
        numbering, percent_identity, chain_type, err_message = _LC_ANNOTATOR.analyze_seq(seq)
        assert chain_type in ['K', 'L'], f"Chain type {chain_type} != K or L"
    return numbering, percent_identity


def anno_cdr(seq: str, chain: str) -> dict:
    """
    Annotate Antibody CDR regions based antpack package

    Parameters
    ----------
    seq : str
        Amino acid sequence of antibody

    chain : str
        Chain type, 'H' for heavy chain, 'K' for kappa light chain, 'L' for lambda light chain

    Returns
    ----------
    dict
        A dictionary containing CDR regions and their sequences
    """

    numbering, percent_identity = numbering_seq(seq, chain)
    if len(numbering) != len(seq):
        raise ValueError(f"Numbering length {len(numbering)} != sequence length {len(seq)}")
    res = [''] * len(IMGT_BORDERS)
    real_borders = [[-1, 0] for _ in range(len(IMGT_BORDERS))]
    region_idx = 0
    for seq_idx, num in enumerate(numbering):
        if num == '-':
            continue
        try:
            if int(num) >= IMGT_BORDERS[region_idx]:
                region_idx += 1
        except ValueError:
            pass
        if real_borders[region_idx][0] == -1:
            real_borders[region_idx][0] = seq_idx
        res[region_idx] += seq[seq_idx]
        real_borders[region_idx][1] = seq_idx + 1
    res = {
        f'{chain}FR1': res[0],
        f'{chain}FR1_slice': slice(*real_borders[0]),
        f'{chain}CDR1': res[1],
        f'{chain}CDR1_slice': slice(*real_borders[1]),
        f'{chain}FR2': res[2],
        f'{chain}FR2_slice': slice(*real_borders[2]),
        f'{chain}CDR2': res[3],
        f'{chain}CDR2_slice': slice(*real_borders[3]),
        f'{chain}FR3': res[4],
        f'{chain}FR3_slice': slice(*real_borders[4]),
        f'{chain}CDR3': res[5],
        f'{chain}CDR3_slice': slice(*real_borders[5]),
        f'{chain}FR4': res[6],
        f'{chain}FR4_slice': slice(*real_borders[6]),
        f'{chain}_percent_identity': percent_identity
    }
    return res


@require_anarci
def anno_tcr_cdr(seq: str) -> dict:
    """
    Annoate TCR CDR info.
    """
    _, numbered, alignment_details, _ = run_anarci([('A', seq)], scheme='imgt')
    assert len(numbered) == len(alignment_details) == 1
    assert numbered[0] is not None, 'No domain found!'
    assert len(numbered[0]) == len(alignment_details[0]) == 1, 'Multiple domain found!'
    numbered_seq, fv_start, fv_end = numbered[0][0]
    alignment_details = alignment_details[0][0]
    chain_type = alignment_details['chain_type']
    if chain_type not in ['B', 'A', 'D', 'G']:
        raise RuntimeError(f'{seq} is not a valid TCR chain')
    
    numbered_seq = filter(lambda x: x[1] != '-', numbered_seq)
    real_borders = [None for _ in range(len(IMGT_BORDERS))]
    region_idx = -1
    for seq_idx, ((aligned_idx, inscode), aa) in enumerate(numbered_seq, fv_start):
        assert seq[seq_idx] == aa, 'AA not matched!'
        if region_idx == -1:
            # find first region
            tmp_idx = 0
            while aligned_idx >= IMGT_BORDERS[tmp_idx]:
                tmp_idx += 1
            region_idx = tmp_idx
        elif aligned_idx >= IMGT_BORDERS[region_idx]:
            region_idx += 1

        if real_borders[region_idx] is None:
            real_borders[region_idx] = [seq_idx, seq_idx + 1]
        else:
            real_borders[region_idx][1] = seq_idx + 1
    
    assert seq_idx == fv_end
    result = {
        'chain_type': chain_type,
        'raw_seq': seq,
        'fv_seq': seq[fv_start: fv_end + 1]
        }
    for region_idx, border in enumerate(real_borders):
        name = REGION_NAMES[region_idx]
        if border is None:
            result[name] = None
        else:
            border = slice(*border)
            result[name] = seq[border]
        result[name + '_slice'] = border
    return result


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    remove_constant_region_parser = subparsers.add_parser("rm_const")
    remove_constant_region_parser.add_argument("--fasta", "-i", type=Path, required=True)
    remove_constant_region_parser.add_argument("--strict", "-s", action="store_true")
    remove_constant_region_parser.add_argument("--ignore_error", "-e", action="store_true")
    remove_constant_region_parser.add_argument("--output", "-o", type=Path, required=True)

    annotate_chain_type_parser = subparsers.add_parser("anno_chain")
    annotate_chain_type_parser.add_argument("--fasta", "-i", type=Path, required=True)
    annotate_chain_type_parser.add_argument("--output", "-o", type=Path, required=True)

    args = parser.parse_args()

    if args.cmd == "rm_const":
        df = pd.DataFrame(remove_constant_region(args.fasta, args.strict, args.ignore_error))
        df2fasta(df, args.output, 'seq', id_col='id')
    elif args.cmd == "anno_chain":
        df = pd.DataFrame(annotate_chain_type(args.fasta))
        df.to_csv(args.output, index=False)
    else:
        raise ValueError(f"Unknown subcommand: {args.cmd}")