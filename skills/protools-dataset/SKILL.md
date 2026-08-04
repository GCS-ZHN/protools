---
name: protools-dataset
description: "Work with protein ML datasets in the protools library (`protools4py`)'s `protools.dataset` package: reproducible train/validation/test splits (`random_split` with explicit seeds, `grouped_split` that keeps sequence families/clusters intact with partial-random / random / sorted modes), plus any dataset utilities added to that package later. Trigger whenever the user wants train/val/test partitions, reproducible dataset splitting, clustered or grouped splitting, stratified-by-family partitioning, dataset deduplication or other dataset preparation in the protools repo, even if they do not name a specific module (e.g. dataset/split)."
---

# protools-dataset

Dataset preparation for ML workloads in the `protools` package. The
`protools.dataset` package currently contains one module, `split`
(`random_split` for plain shuffling, `grouped_split` for family-aware /
cluster-aware splits); more dataset utilities are expected to be added to
the package later. Keep this skill as the entry point for the whole package,
not just its current module.

## How to work in this repo

- Run code with `uv run python`, tests with `uv run pytest test/dataset/`.
- Look up exact signatures with CodeGraph or the module source before coding.
- Splits are the foundation of model training — always make them
  reproducible (explicit `seed`) and deterministic for a given seed/version.
- Real data over mocks: derive group sizes from real clustering output
  (e.g. `protools.cluster`) or fixtures in `data/`; don't invent protein
  families.
- The package is deliberately small; when new dataset functionality is
  added, put it under `protools/dataset/`, document it in the module
  docstring, and update `docs/architecture/overview.md` (the source of
  truth for module layout) in the same change.

## `random_split`

```python
from protools.dataset.split import random_split

train, val, test = random_split(
    list_of_sequences,        # any iterable; consumed once
    lengths=[800, 100, 100],  # sum must be <= len(data)
    seed=42,
)
```

- Returns a list of lists, one per requested length, sampled without
  replacement using `np.random.default_rng(seed)`.
- Every length must be a positive integer; a sum larger than the dataset
  raises `ValueError`.
- Order within each split is the shuffled order — don't assume the input
  order is preserved.

## `grouped_split`

Use this when items belong to groups (e.g. CD-HIT/MMseqs2 clusters, the same
protein family) and you want all members of a group to land in one split —
this prevents family leakage between train and test.

```python
from protools.dataset.split import grouped_split

# group_sizes: number of items in each group (cluster sizes)
split_group_ids = grouped_split(
    group_sizes=[120, 45, 210, 33, 90, 300],
    lengths=[400, 200, 198],      # must sum to sum(group_sizes)
    mode="partial_random",        # 'partial_random' | 'random' | 'sorted'
    seed=7,
    chunk_size=3,
)
# -> list of lists of group indices, e.g. [[2, 4], [0, 1, 3], [5]]
```

Semantics:

- The sum of `lengths` must **equal** the sum of `group_sizes`; the greedy
  assignment then minimizes the squared deviation of each split's final size
  from its target.
- `mode="partial_random"` (default): groups sorted by size descending, then
  shuffled within chunks of `chunk_size`, then greedily assigned — a
  balance of determinism and randomization.
- `mode="random"`: fully shuffled; `mode="sorted"`: strictly by descending
  size (deterministic).
- The result is a list of *group indices*, not items. Map back to your items
  with e.g. `[item for i in split_ids for item in groups[i]]` (or use the
  index list directly with your data structure).

## Typical pipeline

```python
seqs = read_fasta("data/vdomain.fasta")
cluster_df = MMseqs2().cluster(seqs, min_seq_id=0.9, coverage=0.8)
sizes = cluster_df.groupby("cluster_id").size().sort_index().tolist()
train_g, val_g, test_g = grouped_split(sizes, [0.8, 0.1, 0.1], seed=42)
```

## Pitfalls

- `random_split` does not know about groups; use `grouped_split` whenever
  leakage by sequence family matters for the task.
- Both functions require positive integer lengths; `grouped_split` also
  asserts the length totals match, so floating-point splits
  (e.g. `0.8 * total`) must be rounded explicitly.
- Seeding is per-call via `np.random.default_rng(seed)`; the same seed gives
  the same split for the same input order — sort inputs when order
  independence matters.
- Tests (`test/dataset/test_split.py`) exercise edge cases like oversized
  totals, zero/negative lengths, and `chunk_size` bounds; keep those green
  when extending the module.
