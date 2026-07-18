# protools — Architecture Overview

Single source of truth for project layout, module boundaries,
and key design notes. AGENTS.md and README.md link here so they
do not duplicate this content.

## Layout

```
protools/
  __init__.py          # re-exports: Fasta (seqio), Intervals (utils)
  seqio.py             # Fasta (OrderedDict[str, SeqRecord]) + read_fasta/save_fasta/temp_fasta
  seqanno.py seqconvert.py seqfixer.py seqfilter.py seqplot.py aa.py
  pdbio.py pdbanno.py pdbconvert.py pdbextract.py pdbfixer.py pdbplot.py
  dock.py              # HDock wrapper + DockQ scoring
  cluster/             # cdhit.py, mmseqs.py — CmdWrapperBase for external CLIs
  database/            # uniprot.py (ExPASy + Bio.SwissProt), patent.py
  dataset/split.py     # random_split, scaffold_split, etc.
  utils.py             # CmdWrapperBase, Intervals, ensure_path/ensure_fileio,
                       # auto_compression, require_package, progress_io (rich)
  typedef.py           # FilePathType / FilePathOrIOType / SeqLikeType
test/                  # mirrors protools/ structure; helper: test/tools.py
data/                  # committed test fixtures: 4zqk.cif, 4ozi.pdb,
                       # vdomain.fasta, ab1218*.a3m, test_protein.fasta, ...
```

## Module groupings (public surface)

| Group        | Modules |
|--------------|---------|
| Sequence     | `protools.seqio`, `seqanno`, `seqconvert`, `seqfixer`, `seqfilter`, `seqplot`, `aa` |
| Structure    | `protools.pdbio`, `pdbanno`, `pdbconvert`, `pdbextract`, `pdbfixer`, `pdbplot`, `dock` |
| Dataset      | `protools.dataset.split` |
| Cluster      | `protools.cluster.cdhit`, `protools.cluster.mmseqs` |
| Database     | `protools.database.uniprot`, `protools.database.patent` |

`protools/__init__.py` re-exports only `Fasta` (from `seqio`) and
`Intervals` (from `utils`); everything else is accessed via the
full module path.

## Key entrypoints & design notes

- **External CLIs** (`hdock`, `createpl`, `cd-hit`, `mmseqs`) are
  invoked via `protools.utils.CmdWrapperBase`. The lookup pattern
  is `${CMD_UPPER}_PATH` env var first, then `shutil.which`.
- **`Fasta`** (in `seqio.py`) subclasses `OrderedDict[str, SeqRecord]`
  and additionally accepts `int`, `slice`, and `Iterable[str]`
  indexing. Non-string keys are coerced to `str` with a
  `UserWarning`.
- **`Intervals`** (in `utils.py`) represents intervals with
  `zero_based` and `end_inclusive` flags. Operators: `~x`
  (invert), `x & y` (intersect), `x | y` (union), `x + n` /
  `x - n` (shift).
- **Sequence / structure I/O** (`seqio`, `pdbio`) wrap BioPython
  (`Bio.SeqIO`, `Bio.PDB`) and `gemmi` for mmCIF/PDB. Compression
  handling lives in `utils.auto_compression`.
- **Database access** (`database/uniprot.py`) uses `Bio.ExPASy` +
  `Bio.SwissProt` over the network — failures from 4xx/5xx or rate
  limits are infrastructure, not bugs.
- **Dock** (`dock.py`) wraps HDock for docking and DockQ for
  scoring.

## Fixtures

`data/` holds committed test fixtures (PDB/CIF/Fasta/A3M). Notable
files:

- `4zqk.cif`, `4zqk.cif.gz`, `4ozi.pdb`, `4ozi.pdb.gz`, `4i77.pdb`,
  `3inj.pdb` — real RCSB structures used by `test/test_pdbio.py`,
  `test/test_pdbfixer.py`, etc.
- `vdomain.fasta`, `vdomain_anno.csv` — antibody V-domain sequences
  + annotation.
- `ab1218*.a3m` — antibody sequence/structure alignments using IMGT
  / Kabat numbering (see Hard rules in `AGENTS.md`).
- `test_protein.fasta`, `test_dna.fasta`, `test_rna.fasta`,
  `test.fasta`, `test.fasta.gz` — generic sequence fixtures.

When adding a fixture, cite the source URL in the commit message
and webfetch-verify it before relying on it in tests.

## Update policy

This document is the source of truth for **layout, module
boundaries, and key entrypoints**. Update it whenever any of those
change (new module, new external CLI, new subpackage, change to
public surface). Do not duplicate this content in `AGENTS.md` or
`README.md` — link here instead.