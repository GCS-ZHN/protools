#!/usr/bin/env python3
import pandas as pd
import logging

from Bio import ExPASy, SwissProt
from Bio.Seq import Seq
from pathlib import Path
from ..utils import max_retry, FilePath, ensure_path, catch_error
from multiprocessing import Pool
from itertools import starmap, product
from Bio.SeqFeature import ExactPosition


loger = logging.getLogger(__name__)


def extract_features(record: SwissProt.Record):
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


@catch_error(loger, err_types=(Exception,))
@max_retry(max=5)
def fetch_uniprot_features(accession_id: str, output_path: FilePath, skip_exist: bool = False):
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
    handle = ExPASy.get_sprot_raw(accession_id)
    record = SwissProt.read(handle)
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
