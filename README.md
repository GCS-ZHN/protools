# protools

Personal protein processing tool collections. Python library covering
sequence / structure I/O, annotation, fix/convert, clustering,
AlphaFold & PDB databases, and ML dataset splits.

## Install

PyPI release is published on every GitHub release:

```bash
uv pip install protools4py
```

## Usage

For the full module list and layout, see
[`docs/architecture/overview.md`](docs/architecture/overview.md). A
quick index:

- **Sequence** — `protools.seqio`, `seqanno`, `seqconvert`,
  `seqfixer`, `seqfilter`, `seqplot`, `aa`
- **Structure** — `protools.pdbio`, `pdbanno`, `pdbconvert`,
  `pdbextract`, `pdbfixer`, `pdbplot`, `dock`
- **Dataset** — `protools.dataset.split`
- **Cluster** — `protools.cluster.cdhit`, `protools.cluster.mmseqs`
- **Database** — `protools.database.uniprot`,
  `protools.database.patent`

## Development

Build / install / test with `uv`:

```bash
uv sync --locked --all-extras --dev
uv run pytest
uv build   # produces dist/ for PyPI
```

## Skills

Agent skills for working with each protools module group ship in
[`skills/`](skills/) (module layout: see
[`docs/architecture/overview.md`](docs/architecture/overview.md)). Install
them into your Codex/Claude skills directory with:

```bash
npx skills add GCS-ZHN/protools
```

Skills cover sequence, structure, cluster, dataset, database, and docking
workflows; each ships with `evals/` test prompts.

Release pipeline (tag → PyPI) and agent workflows are documented
in [`AGENTS.md`](AGENTS.md).

## TODO

- [ ] Add per-module usage examples.
