---
name: protools-cluster
description: "Cluster or deduplicate sequences with CD-HIT and MMseqs2 through the protools library (`protools4py`): redundancy removal within one FASTA, 2D cross-dataset clustering, parsing .clstr cluster files into DataFrames, and MMseqs2 clustering with identity/coverage control returning representative-per-cluster tables. Trigger whenever the user wants to deduplicate FASTA files, cluster sequences, parse .clstr output, or use cd-hit / cd-hit-2d / mmseqs in the protools repo, even if they do not name the module (cluster/cdhit, cluster/mmseqs)."
---

# protools-cluster

Sequence clustering and redundancy removal in the `protools` package. Two
backends are wrapped: CD-HIT (`cluster/cdhit.py`) and MMseqs2
(`cluster/mmseqs.py`). Both run external command-line tools through
`utils.CmdWrapperBase`.

## How to work in this repo

- Run code with `uv run python`, tests with `uv run pytest test/cluster/`.
- External binaries are located via the `${CMD}_PATH` environment variable
  first, then `shutil.which`; a missing binary raises `CmdNotFoundError`
  whose message includes install instructions. Set e.g. `CD_HIT_PATH` when
  the binary lives outside `PATH`.
- Look up exact signatures with CodeGraph or the module source before coding.
- Cluster results are deterministic only for a fixed tool version and
  parameters — pin versions when a test depends on exact output.
- Real data over mocks: use FASTA fixtures from `data/` and the `Fasta`
  reader (`read_fasta`) rather than fabricating sequences.

## CD-HIT (`protools.cluster.cdhit`)

```python
from pathlib import Path
from protools.cluster.cdhit import CdHit, CdHit2D, unique_fasta, unique_fasta_2d, parse_cluster

# 1) Deduplicate a FASTA (also writes <output>.cluster.csv)
unique_fasta(Path("all.fasta"), Path("unique.fasta"), cutoff=0.9, num_threads=4)

# 2) 2D mode: keep source sequences that have no hit in target
unique_fasta_2d(Path("source.fasta"), Path("target.fasta"), Path("novel.fasta"))

# 3) Parse an existing .clstr file
df = parse_cluster(Path("out.clstr"))
# columns: cluster_id, sequence_id, sequence_length, representative, similarity
```

Lower-level wrappers:

- `CdHit(cmd="cd-hit", num_workers=1)` — call with keyword args mapped to
  CLI flags: `i`, `o`, `c` (cutoff, default 0.9), `n` (word length, default
  5), `M` (memory MB, default 800), `T` (threads, default 1), plus any extra
  `cd-hit` flags.
- `CdHit2D(cmd="cd-hit-2d")` — same idea for `-i2` cross-mode.
- The temporary-file plumbing in `unique_fasta*` renames ids to `0,1,2...`
  and maps them back afterwards, so **original ids are preserved in the
  cluster CSV** — don't re-map ids yourself.

## MMseqs2 (`protools.cluster.mmseqs`)

```python
from protools.cluster.mmseqs import MMseqs2
from protools.seqio import read_fasta

seqs = read_fasta("data/vdomain.fasta")
df = MMseqs2(num_workers=2).cluster(
    seqs,
    output=Path("clusters.csv"),     # optional; also returned as DataFrame
    min_seq_id=0.9,                  # --min-seq-id
    coverage=0.8,                    # --cov-mode 0 coverage on both
    coverage_mode=0,                 # 0: query+target, 1: target, 2: query
    threads=4,
)
```

- `cluster(seqs, output=None, tmp_dir=None, min_seq_id=0.9, coverage=0.8,
  coverage_mode=0, threads=1, **kwargs)` returns a DataFrame with columns
  `cluster_id`, `sequence_id`; the first row of each cluster is the
  representative (`first_seq_as_repr=1`).
- Input goes in through `mmseqs createdb` via stdin, so `seqs` must be a
  `Fasta` (or anything with `.to_fasta_str()` and `.keys()`).
- Intermediate files live in a temporary directory (`tmp_dir` if given);
  nothing is written to the working tree except the optional `output`.

## Pitfalls

- Both tools are external binaries; install them first
  (cdhit: https://github.com/weizhongli/cdhit,
  mmseqs: https://github.com/soedinglab/MMseqs2). The wrappers surface
  `CmdNotFoundError` with those install hints — don't suppress it.
- `unique_fasta` writes two artifacts: the deduplicated FASTA and a
  `<output>.cluster.csv` — keep both when a pipeline expects the cluster
  table.
- CD-HIT sequence-identity cutoffs and MMseqs2 `min_seq_id`/`coverage` are
  threshold choices that materially change results; surface them as
  explicit parameters rather than hard-coding new defaults.
- 2D mode semantics differ from 1D: it removes *source* sequences that hit
  *target* sequences; the output FASTA contains only surviving source
  sequences. Do not confuse it with plain deduplication.
