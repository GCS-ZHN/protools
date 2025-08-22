#!/usr/bin/env python3
import pandas as pd
import logging

from Bio import ExPASy, SwissProt
from Bio.Seq import Seq
from pathlib import Path
from ..utils import (
    max_retry, FilePathType, ensure_path, catch_error, Intervals)
from multiprocessing import Pool
from itertools import starmap, product
from Bio.SeqFeature import ExactPosition


loger = logging.getLogger(__name__)


def get_uniprot_record(accession_id: str) -> SwissProt.Record:
    """
    Extracts the SwissProt record from the given accession id.

    Parameters
    ----------
    accession_id : str
        The accession id of the protein in UniProt.

    Returns
    ----------
    SwissProt.Record
        The SwissProt record of the protein.
    """
    handle = ExPASy.get_sprot_raw(accession_id)
    record = SwissProt.read(handle)
    return record


def extract_features(record: SwissProt.Record):
    """
    Extracts the features from the SwissProt record.
    """
    loger.debug(f'Record: {record.accessions}')
    seq = Seq(record.sequence)
    loger.debug(f'Sequence: {seq}')
    for feature in record.features:
        feature_type = feature.type
        loger.debug(f'Feature type: {feature_type}')
        pos_start = feature.location.start
        loger.debug(f'Feature start: {pos_start}')
        pos_end = feature.location.end
        pos_str = feature.location.__str__()
        loger.debug(f'Feature end: {pos_end}')
        if isinstance(pos_start, ExactPosition) and isinstance(pos_end, ExactPosition):
            sub_seq = seq[pos_start: pos_end]
        else:
            sub_seq = ''
        loger.debug(f'Feature subseq: {sub_seq}')
        qualifiers = feature.qualifiers
        loger.debug(f'Feature qualifiers: {qualifiers}')
        res = {
            'type': feature_type,
            'pos_range': pos_str,
            'subseq': str(sub_seq)
        }
        res.update(qualifiers)
        yield res


def extract_extracellular_chains(record: SwissProt.Record) -> pd.DataFrame:
    """
    Extracts the extracellular chains from the SwissProt record.

    Parameters
    ----------
    record : SwissProt.Record
        The SwissProt record.
    
    Returns
    ----------
    pd.DataFrame
        The extracellular chains. there are three columns: 
        'seq', 'slice', 'exact'. 'seq' is the sequence of
        the extracellular chain, 'slice' is the slice of the
        extracellular chain in the full sequence, 'exact' is a
        boolean value indicating whether the slice is from the
        exact extracellular region or not.
    """
    non_extracellular_slices = []
    extracellular_slices = []
    signal_slice = None
    for feature in record.features:
        pos_start = feature.location.start
        if not isinstance(pos_start, ExactPosition):
            continue
        pos_end = feature.location.end
        if not isinstance(pos_end, ExactPosition):
            continue
        if feature.type == 'SIGNAL':
            if signal_slice is not None:
                raise ValueError('Multiple signal peptide found')
            signal_slice = slice(pos_start, pos_end)
            non_extracellular_slices.append(signal_slice)
        if feature.type == 'TRANSMEM':
            non_extracellular_slices.append(slice(pos_start, pos_end))
        if feature.type == 'TOPO_DOM':
            if feature.qualifiers['note'] == 'Extracellular':
                extracellular_slices.append(slice(pos_start, pos_end))
            if feature.qualifiers['note'] == 'Cytoplasmic':
                non_extracellular_slices.append(slice(pos_start, pos_end))

    results = []
    if len(extracellular_slices) > 0:
        for extracellular in extracellular_slices:
            extracellular_seq = record.sequence[extracellular]
            results.append({
                'seq': extracellular_seq,
                'slice': extracellular,
                'exact': True,
            })
        return pd.DataFrame(results)

    for extracellular in Intervals.from_slices(
            non_extracellular_slices).invert(high=len(record.sequence)):
        extracellular_seq = record.sequence[extracellular]
        results.append({
            'seq': extracellular_seq,
            'slice': extracellular,
            'exact': False,
        })

    return pd.DataFrame(results)


@catch_error(loger, err_types=(Exception,))
@max_retry(max=5)
def fetch_uniprot_features(accession_id: str, output_path: FilePathType, skip_exist: bool = False):
    main_id = accession_id.split('-')[0]
    isoform = accession_id.split('-')[1] if len(accession_id.split('-')) > 1 else ''
    if not main_id.isalnum():
        raise ValueError(f'Invalid accession id: "{accession_id}"')
    output_path = ensure_path(output_path)
    output_path.mkdir(exist_ok=True, parents=True)
    output_file = output_path / f'{accession_id}_features.txt'
    if output_file.exists() and skip_exist:
        return

    loger.info(f'Fetching {accession_id}')
    record = get_uniprot_record(accession_id)
    features = pd.DataFrame(extract_features(record))
    features.to_csv(output_file, sep='\t', index=False)
    loger.info(f'Writing to {output_file}')


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('accession_ids', type=str, nargs='+')
    parser.add_argument('--output_path', '-o', type=Path, required=True)
    parser.add_argument('--n_jobs', '-n', type=int, default=1)
    parser.add_argument('--log_level', '-l', type=str, default='INFO')
    parser.add_argument('--skip_exist', '-s', action='store_true')

    args = parser.parse_args()
    logging.basicConfig(level=logging._nameToLevel[args.log_level], 
                        format='%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s')
    acessions = map(lambda x: x.strip(), args.accession_ids)
    schedules = product(acessions, [args.output_path], [args.skip_exist])

    if args.n_jobs > 1:
        with Pool(args.n_jobs) as p:
            p.starmap(fetch_uniprot_features, schedules)

    else:
        list(starmap(fetch_uniprot_features, schedules))
