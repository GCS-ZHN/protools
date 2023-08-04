from typing import overload, Dict, Iterable, Tuple


@overload
def save_fasta(seq_dict: Dict[str, str], fasta_file: str):
    ...


@overload
def save_fasta(seq_iter: Iterable[Tuple[str, str]], fasta_file: str):
    ...