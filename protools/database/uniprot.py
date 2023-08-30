#!/usr/bin/env python3
import pandas as pd
import logging

from Bio import ExPASy, SwissProt
from Bio.Seq import Seq
from pathlib import Path
from ..utils import max_retry, FilePath, ensure_path
from multiprocessing import Pool
from itertools import starmap, product


loger = logging.getLogger(__name__)


def extract_features(record: SwissProt.Record):
    seq = Seq(record.sequence)

    for feature in record.features:
        feature_type = feature.type
        pos_start = feature.location.start
        pos_end = feature.location.end
        sub_seq = seq[pos_start: pos_end]
        qualifiers = feature.qualifiers
        res = {
            'type': feature_type,
            'subseq': str(sub_seq)
        }
        res.update(qualifiers)
        yield res


@max_retry(max=5)
def fetch_uniprot_features(accession_id: str, output_path: FilePath):
    loger.info(f'Fetching {accession_id}')
    handle = ExPASy.get_sprot_raw(accession_id)
    record = SwissProt.read(handle)
    output_path = ensure_path(output_path)
    output_path.mkdir(exist_ok=True, parents=True)
    features = pd.DataFrame(extract_features(record))
    output_file = output_path / f'{accession_id}_features.txt'
    features.to_csv(output_file, sep='\t', index=False)
    loger.info(f'Writing to {output_file}')


if __name__ == '__main__':
    from argparse import ArgumentParser
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s')
    parser = ArgumentParser()
    parser.add_argument('accession_ids', type=str, nargs='+')
    parser.add_argument('--output_path', '-o', type=Path, required=True)
    parser.add_argument('--n_jobs', '-n', type=int, default=1)

    args = parser.parse_args()
    acessions = map(lambda x: x.strip(), args.accession_ids)
    schedules = product(acessions, [args.output_path])

    if args.n_jobs > 1:
        with Pool(args.n_jobs) as p:
            p.starmap(fetch_uniprot_features, schedules)

    else:
        starmap(fetch_uniprot_features, schedules)
