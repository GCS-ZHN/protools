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

from antpack import SingleChainAnnotator, VJGeneTool

_ANTIBODY_HC_IDS = ['H']
_ANTIBODY_LC_IDS = ['K', 'L']
_TCR_HC_IDS = ['B', 'D']
_TCR_LC_IDS = ['A', 'G']
_ANTIBODY_HC_ANNOTATOR = SingleChainAnnotator(chains=_ANTIBODY_HC_IDS, scheme='imgt')
_ANTIBODY_LC_ANNOTATOR = SingleChainAnnotator(chains=_ANTIBODY_LC_IDS, scheme='imgt')
_TCR_HC_ANNOTATOR = SingleChainAnnotator(chains=_TCR_HC_IDS, scheme='imgt')
_TCR_LC_ANNOTATOR = SingleChainAnnotator(chains=_TCR_LC_IDS, scheme='imgt')


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


def numbering_seq(seq: str, chain: str, is_tcr: bool = False) -> tuple:
    if chain == 'H':
        annnotator = _TCR_HC_ANNOTATOR if is_tcr else _ANTIBODY_HC_ANNOTATOR
        valid_chains = _TCR_HC_IDS if is_tcr else _ANTIBODY_HC_IDS
    else:
        annnotator = _TCR_LC_ANNOTATOR if is_tcr else _ANTIBODY_LC_ANNOTATOR
        valid_chains = _TCR_LC_IDS if is_tcr else _ANTIBODY_LC_IDS
    
    numbering, percent_identity, chain_type, err_message = annnotator.analyze_seq(seq)
    assert chain_type in valid_chains, f"Chain type {chain_type} not in {valid_chains}, detail: {err_message}"
    return numbering, percent_identity, chain_type, err_message


def anno_cdr(seq: str, chain: str, is_tcr: bool = False) -> dict:
    """
    Annotate Antibody CDR regions based antpack package

    Parameters
    ----------
    seq : str
        Amino acid sequence of antibody

    chain : str
        Chain type, 'H' for heavy chain, 'K' for kappa light chain, 'L' for lambda light chain

    is_tcr: str
        Whether the sequence is TCR or not. Default is False. supported since antpack v0.3.8

    Returns
    ----------
    dict
        A dictionary containing CDR regions and their sequences
    """

    numbering, percent_identity, _, _ = numbering_seq(seq, chain, is_tcr=is_tcr)
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


def anno_vj_gene(seq: str,
                 chain_type: str,
                 species: str,
                 is_tcr: bool = False) -> dict:
    """
    Annotate VJ gene for antibody or TCR sequence.
    
    Parameters
    ----------
    seq : str
        Amino acid sequence of antibody or TCR.
    
    chain_type : str
        Chain type, 'H' for heavy chain, 'K' for kappa light chain, 'L' for lambda light chain,
        'A' for alpha TCR chain, 'B' for beta TCR chain, 'D' for delta TCR chain, 'G' for gamma TCR chain.
    species : str
        Species of the sequence, can be 'human', 'mouse', 'alpaca', 'rabbit', or 'unknown'.
        if is_tcr is True, species must be 'human', 'mouse', or 'unknown'.
    is_tcr : bool
        Whether the sequence is TCR or not. Default is False.
    Returns
    -------
    dict
        A dictionary containing VJ gene information.
    """
    vj_tool = VJGeneTool()
    alignment = numbering_seq(seq, chain_type, is_tcr=is_tcr)
    if is_tcr:
        assert species in ['human', 'mouse', 'unknown'], \
            "Species must be 'human', 'mouse' or 'unknown' for TCR"
    else:
        assert species in ['human', 'mouse', 'alpaca', 'rabbit', 'unknown'], \
            "Species must be 'human', 'mouse', 'alpaca', 'rabbit' or 'unknown' for antibody"
    v_gene, j_gene, v_pident, j_pident, species_matched = vj_tool.assign_vj_genes(
        alignment=alignment,
        sequence=seq,
        species=species,
        mode='identity'
    )
    if species != 'unknown':
        assert species_matched == species, \
            f"Species {species_matched} does not match the input species {species}"
    v_genes = v_gene.split('_')
    j_genes = j_gene.split('_')
    return {
        'v_genes': v_genes if len(v_genes) > 1  else v_genes[0],
        'j_genes': j_genes if len(j_genes) > 1 else j_genes[0],
        'v_pident': v_pident,
        'j_pident': j_pident,
        'species': species_matched,
    }


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