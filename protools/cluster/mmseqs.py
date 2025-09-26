from protools.utils import CmdWrapperBase
from protools import seqio
from protools import typedef
from typing import Optional
import pandas as pd
import tempfile


class MMseqs2(CmdWrapperBase):
    """
    A Basic wrapper for MMseqs2 command line tool.
    """
    def __init__(self, cmd: str = 'mmseqs', num_workers: int = 1):
        super().__init__(
            cmd,
            num_workers=num_workers,
            underline2hyphen=True,
            install_help="Please refer to https://github.com/soedinglab/MMseqs2 for installation instructions.")
        

    def cluster(self, seqs: seqio.Fasta,
                output: Optional[typedef.FilePathType] = None,
                tmp_dir: Optional[typedef.FilePathType] = None,
                min_seq_id: float = 0.9, coverage: float = 0.8,
                coverage_mode: int = 0, threads: int = 1, **kwargs) -> pd.DataFrame:
        """
        Cluster sequences using MMseqs2.

        Parameters
        ----------
        seqs : seqio.Fasta
            Input sequences in Fasta format.
        output : Optional[typedef.FilePathType], optional
            Output file path to save clustering results, by default None.
        tmp_dir : Optional[typedef.FilePathType], optional
            Temporary directory for intermediate files, by default None.
        min_seq_id : float, optional
            Minimum sequence identity for clustering, by default 0.9.
        coverage : float, optional
            Minimum coverage for clustering, by default 0.8.
        coverage_mode : int, optional
            Coverage mode (0: coverage of query and target, 1: coverage of target, 2: coverage of query), by default 0.
        threads : int, optional
            Number of threads to use, by default 1.
        **kwargs
            Additional command line arguments for MMseqs2.

        Returns
        -------
        pd.DataFrame
            DataFrame containing clustering results.
        """
        with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
            db_path = f'{work_dir}/db'
            clu_path = f'{work_dir}/clu'
            tsv_path = f'{work_dir}/result.tsv'
            pipe_input = ''.join(seqs.to_fasta_str(two_line_mode=True))
            self.createdb('stdin', 
                          db_path,
                          pipe_input = pipe_input.encode('utf-8'))
            self('cluster', db_path, clu_path, f'{work_dir}/clu-tmp',
                 min_seq_id=min_seq_id,
                 c=coverage,
                 threads=threads,
                 cov_mode=coverage_mode,
                 **kwargs)
            idx2name = dict(enumerate(seqs.keys()))
            self.createtsv(db_path, clu_path, tsv_path,
                        threads=threads, first_seq_as_repr=1)
            df = pd.read_table(tsv_path, header=None,
                             names=['cluster_id', 'sequence_id'])
            df['cluster_id'] = df['cluster_id'].map(idx2name)
            df['sequence_id'] = df['sequence_id'].map(idx2name)
            for n, d in df.groupby('cluster_id'):
                assert n == d['cluster_id'].iloc[0], \
                f"cluster {d}'s first seq is not repr seq {n}"
            if output:
                df.to_csv(output, index=False)
            return df