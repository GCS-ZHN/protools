---
name: protools-dock
description: "Run protein-protein docking and scoring with the protools library (`protools4py`): HDock docking between receptor and ligand PDBs with optional active-site restraints (rsite/lsite), building docked complex PDBs via createpl, and DockQ scoring of models against native structures. Trigger whenever the user wants to dock two proteins, build a docked complex structure, score docking models with DockQ, or prepare HDock input/output in the protools repo, even if they do not name the module (dock.py, HDock, dockq_score)."
---

# protools-dock

Docking and scoring in the `protools` package, wrapped in
`protools.dock`:

- `HDock` — wraps the external `hdock` binary plus `createpl` to build
  complex PDBs.
- `dockq_score` — scores a model against a native structure with DockQ.

## How to work in this repo

- Run code with `uv run python`, tests with `uv run pytest test/test_dock.py`.
- External binaries are resolved via `${CMD}_PATH` env vars first
  (`HDOCK_PATH`, `CREATEPL_PATH`), then `shutil.which`; a missing binary
  raises `CmdNotFoundError`. HDock docs: http://hdock.phys.hust.edu.cn/
- Docking conventions are the user's domain: chain maps, interface
  residues, active-site restraints, and reference structures must come from
  the user or a cited source — never invent them.
- Look up exact signatures with CodeGraph or the module source before coding.
- Real data over mocks: dock with real receptor/ligand structures from RCSB
  fixtures or `pdbio.fetch`; cite the PDB IDs.

## Docking with HDock

```python
from protools.dock import HDock

hdock = HDock()
dockout = hdock.dock(
    ligand_pdb="ligand.pdb",
    receptor_pdb="receptor.pdb",
    output=Path("result/run1.dockout"),   # optional; defaults to temp file
    itscore=True,                          # ITScore scoring
    rsite=Path("rsite.txt"),               # optional receptor active-site restraints
    lsite=Path("lsite.txt"),               # optional ligand restraints
    angle=15,                              # sampling angle
)
```

Then build the complex structure(s):

```python
hdock.create_complex(
    dockout,
    pdb_name="result/run1_complex.pdb",
    model_num=100,      # -nmax
    rmsd=5.0,           # clustering RMSD
    complex=True,       # -complex flag
    models=True,        # -models flag
    chid=False,         # -chid flag
)
```

Details worth knowing:

- `dock()` copies inputs into a temp working directory, runs `hdock`, then
  rewrites the receptor/ligand paths inside the `.dockout` to be relative to
  the output — so keep the output file together with the input PDBs when
  downstream tools read the dockout.
- `create_complex()` re-locates those paths again (it copies the dockout
  and PDBs next to the output) and runs `createpl`; `pdb_name` is written
  into the output directory. `-complex` and `-models` change what is
  emitted; `-chid` preserves chain ids.
- If `output` is omitted, `dock()` returns a path inside a temporary
  directory owned by the `HDock` instance (cleaned up on destruction) —
  pass an explicit `output` when the result must outlive the object.

## Scoring with DockQ

```python
from protools.dock import dockq_score

scores = dockq_score(model, native, chain_map={"A": "A", "B": "C"})
```

- `model` / `native` are `Bio.PDB.Model` objects (parse with
  `pdbio.get_structure(...)[0]`).
- `chain_map` maps **native chain ids (keys) to model chain ids (values)** —
  the reverse of what you might guess; `dockq_score` validates that every
  key exists in `native` and every value exists in `model`. A wrong map
  silently produces wrong scores, so confirm the correspondence with the
  user when it is not obvious from the input files.
- The underlying `DockQ` package is function-based (`DockQ.run_on_all_native_interfaces`,
  `DockQ.load_PDB`), not a callable class — use `protools.dock.dockq_score`
  instead of guessing the package API.
- Returns a dict of DockQ metrics (see `DockQ.DockQ` output fields: DockQ
  score, fractions of native contacts, RMSD, etc.).

## Pitfalls

- `hdock`/`createpl` must be installed and executable; surface
  `CmdNotFoundError`'s install hint instead of silently falling back.
- Docking is stochastic-ish across versions: pin the HDock version when a
  test asserts specific scores or poses.
- `create_complex` writes into `dock_result.parent` — create the output
  directory before calling when the dockout lives elsewhere.
- Never hard-code PDB IDs, chain maps, or interface residues in new code or
  tests; they are dataset-specific facts that must be verified.
