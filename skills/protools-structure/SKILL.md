---
name: protools-structure
description: "Work with protein 3D structures (PDB/mmCIF) in the protools library (`protools4py`): parse and save structures, download from RCSB, convert PDB to FASTA/DataFrame, extract or delete residue ranges by selection string, split chains at linkers, translate/rotate coordinates, compute SASA, contacts and interactions, run TMalign alignment, renumber residues, and align design models with PyMOL. Trigger whenever the user handles .pdb/.cif files, structure I/O, residue extraction, chain splitting, SASA, TM-score alignment, residue renumbering, or structure visualization in the protools repo, even if they do not name the module (pdbio/pdbconvert/pdbextract/pdbanno/pdbfixer/pdbplot)."
---

# protools-structure

Helpers for structure work in the `protools` package: parsing and saving
PDB/mmCIF, extraction and editing of residue ranges and chains, geometric
operations, per-residue analysis, TMalign, residue renumbering, and PyMOL
based alignment.

## How to work in this repo

- Run code with `uv run python`, tests with `uv run pytest test/test_pdb*.py`.
- Look up exact signatures with CodeGraph (`codegraph explore "pdbio save_pdb"`)
  or the module source before writing code; do not guess parameter names.
- `docs/architecture/overview.md` is the source of truth for module layout.
- Real structures over mocks: use committed RCSB fixtures in `data/`
  (`4zqk.cif`, `4ozi.pdb`, `4i77.pdb`, `3inj.pdb`). New fixtures need a real
  source URL, webfetch-verified before use.
- PDB IDs, chain IDs, residue numbering, and antibody numbering schemes are
  the user's domain — never invent them. When a choice isn't pinned by code
  or a fixture, confirm with the user.
- Tests verify outputs via `test/tools.py::md5_equal` against hard-coded
  checksums; recompute the constant when a fixture changes.

## I/O (`pdbio`)

```python
from protools import pdbio

structure = pdbio.get_structure("data/4ozi.pdb")   # Bio.PDB Structure
```

- `get_structure(path, structure_id=None)` parses `.pdb` (PDBParser) or
  `.cif` (MMCIFParser), handles `.gz` transparently, and defaults the
  structure id to the file's md5.
- `save_pdb(output_path, *entities, remarks=None, seqres=None, modres=None)`
  accepts `Structure` | `Model` | `Chain` | `Residue` — all entities must be
  the same type, and `None` acts as a separator when building multi-chain /
  multi-model files. `remarks` is `{remark_number: [lines]}`, `seqres` is
  `{chain_id: one_letter_seq}`, `modres` is a DataFrame with columns
  `id_code, res_name, chain_id, sequence_number, insertion_code,
  standard_res_name, comment`.
- `fetch(pdb_id, target_dir, server="https://files.rcsb.org")` downloads
  PDB/mmCIF from RCSB and returns the local `Path`. Network failures
  (4xx/5xx, rate limits) are infrastructure issues, not code bugs.
- `pdb2fasta(fasta_path, *pdb_files, multimer_mode="joint"|"seperate",
  selected_chains=None, standard=True, unknown_aa="X", joint_sep=":")`
  writes a FASTA (`joint` joins chains with `:`).
- `pdb2seq(path)` returns a `Fasta` (SEQRES first, falls back to parsed
  atoms); `pdb2df(entity, *extra_attrs)` converts an entity to a DataFrame
  (see `RESIDUE_LEVEL_COLS`); `read_residue(pdb, mode="centroid")` extracts
  residue-level geometry.
- Missing/remark helpers: `read_remarks`, `read_missing_residues`,
  `read_missing_atoms`, `generate_missing_residues_remarks` /
  `generate_missing_atoms_remarks` (dicts keyed by remark id),
  `read_modified_residues` / `write_modified_residues`.
- `get_aa_residues(chain)` / `get_aa_sequence(chain, standard=True,
  unknown_aa="X")`; `is_aa` / `is_het` classify residues.
- `coord2chain(coord, seq, chain_id="A", atoms=None, init_residue_number=1,
  init_serial_number=1)` builds a `Chain` from an `(n_res, n_atoms, 3)`
  array (default atoms: `BACKBONE_ATOMS`).

## Editing geometry (`pdbconvert`)

All operations mutate the entity in place.

```python
from protools.pdbconvert import translate, rotate, rand_rotate, chain_split

translate(chain, x_translation=10.0)          # returns the translation vector
rotate(chain, [30, 0, 0], degrees=True, self_rotation=True)
rand_rotate(chain, seed=42)                   # returns the chosen angles
chains = chain_split(chain, "GGGGS")          # split at linker regex
```

- `rotate` angles are `(rx, ry, rz)`; `self_rotation=True` rotates around the
  entity's center of mass instead of the origin.
- `chain_split(chain, linker_pattern)` splits a chain at every match of the
  linker regex (`"GGGGS"`, `"G{3,4}S"`), re-attaching the pieces to the
  parent with fresh chain ids (A, B, C, ...).

## Extraction (`pdbextract`)

```python
from protools.pdbextract import extract, batch_extract, get_interface, find_seq

extract("in.pdb", "out.pdb", "H1A-4,F6,M8-10")        # keep these residues
extract("in.pdb", "out.pdb", "H1A-4", remain=False)   # delete them instead
```

- Selection strings start each range with a chain id, e.g. `"H1A-4,F6,M8-10"`
  (see `parse_resi`). `remain=True` keeps the selection, `False` removes it.
- `batch_extract(pdb_files, out_dir, resi_selection, ..., num_process=1)`
  parallelizes over files.
- `get_interface(model, threshold=8.0, dist_type="ca")` returns a DataFrame
  of inter-chain contacts; `find_seq(seq, model, findall=False)` locates a
  peptide sequence in a structure.

## Analysis (`pdbanno`)

- `distance(...)` — pairwise distance helper (see source; used by
  `get_interface` and `pdbextract`).
- `calc_sasa(...)` / `calc_sasa_from_pdbs(...)` — solvent-accessible surface
  area per residue (Bio.PDB `ShrakeRupley` under the hood); `neighbor_water_count(...)`
  counts water neighbors.
- `get_interactions(entity)` returns `Dict[str, pd.DataFrame]` of interaction
  types (e.g. hydrophobic, H-bonds) between chains.
- `TMalign(num_workers=1, byresi=0, clean=True, tmalign_path="TMalign")` is a
  `CmdWrapperBase` wrapper; calling it returns `(tm_score, rmsd)`:

```python
from protools.pdbanno import TMalign
tm, rmsd = TMalign()("query.pdb", "native.pdb", Path("out_dir"), prefix="aln")
```

It also offers `async_call` and `async_batch_align(pair_iterable, output_dir)`.
Requires the external `TMalign` binary (install help: zhanggroup.org/TM-align).

## Renumbering (`pdbfixer`)

```python
from protools.pdbfixer import renumber_residue
renumber_residue("in.pdb", chain_order=["A", "B"], start=1,
                 restart_by_chain=False, chain_id_map={"B": "L"})
```

Output defaults to `{stem}_renumbered{suffix}`; `chain_order` controls the
order (missing chains are appended), `restart_by_chain=True` restarts at
`start` per chain, `chain_id_map` renames chains. **Insertion codes are
removed** and remark 220 is written.

## Alignment & plotting (`pdbplot`)

```python
from protools.pdbplot import align_all_designs
align_all_designs(Path("designs"), Path("aligned.pdb"),
                  receptor_chain="A", ref_pdb=Path("ref.pdb"))
```

Aligns every `*.pdb` in a directory to a reference (default: first file)
using PyMOL. PyMOL is required and lazy-loaded via
`utils.require_package("pymol", "conda install -c schrodinger pymol")` — do
not silently swap in another aligner.

## Pitfalls

- `save_pdb` requires homogeneous entity types; use `None` separators to
  build multi-chain outputs from residues.
- PDB residue ids include insertion codes (hetatm flags); selection strings
  like `"H1A"` refer to insertion-coded residues. Keep `remain` semantics in
  mind — extraction mutates the structure's chains.
- `get_structure` defaults structure id to md5; pass an explicit id when
  downstream code relies on a stable name.
- Structure fixtures are real RCSB entries; verify any PDB ID against RCSB
  and cite the URL when adding fixtures.
