---
name: protools-sequence
description: "Work with protein, DNA, and RNA sequences in the protools library (`protools4py`): read/write FASTA and A3M, translate, reverse-translate and codon-optimize DNA, deduplicate and split FASTA files, filter sequences, number antibody/TCR sequences (IMGT/Kabat/Martin/Aho), annotate CDR regions and V/J genes, compute sequence identity and mutation lists, and draw amino-acid heatmaps. Trigger whenever the user works with FASTA/A3M files, sequence I/O, sequence conversion or optimization, antibody numbering, CDR annotation, mutation lists, or sequence visualization in the protools repo, even if they do not name the module (seqio/seqanno/seqconvert/seqfixer/seqfilter/seqplot/aa)."
---

# protools-sequence

Helpers for sequence work in the `protools` package: reading and writing
FASTA/A3M, converting and optimizing sequences, cleaning and filtering,
antibody/TCR numbering and CDR annotation, and simple visualization.

## How to work in this repo

- Run code with `uv run python`, tests with `uv run pytest test/test_seq*.py`.
- Before writing code that touches a function, look up its exact signature
  with CodeGraph (`codegraph explore "seqio Fasta read_fasta"`) or read the
  module source; do not guess parameter names or semantics.
- `docs/architecture/overview.md` is the source of truth for module layout;
  update it (and README/CHANGELOG) when a change touches the public surface.
- `.gz` files are handled transparently by `utils.auto_compression` inside the
  I/O functions — never decompress manually.
- Real data over mocks: prefer committed fixtures in `data/` (e.g.
  `vdomain.fasta`, `ab1218*.a3m`). New fixtures must cite a real source URL
  and be webfetch-verified before use.
- Tests often verify output via `test/tools.py::md5_equal` against a
  hard-coded checksum; if a fixture changes, recompute and update the
  constant.
- Antibody/TCR numbering is domain knowledge, not something to improvise.
  This library wraps `antpack`; always validate numbering choices against
  `numbering_seq` and the existing tests rather than inventing schemes.

## Core object: `protools.seqio.Fasta`

`Fasta` subclasses `OrderedDict[str, Bio.SeqRecord]`. Keys are always strings
(non-string keys are coerced with a `UserWarning`), values are coerced to
`SeqRecord` on assignment.

Indexing is richer than a plain dict:

```python
from protools.seqio import Fasta, read_fasta, save_fasta

fa = read_fasta("data/vdomain.fasta")
fa["seq1"]            # -> SeqRecord
fa[0]                 # -> first SeqRecord (int index)
fa[1:3]               # -> sub-Fasta (slice)
fa[["a", "c"]]        # -> sub-Fasta (iterable of keys)
```

Useful methods: `to_dict(*extra_attrs)`, `unique()` (dedupe by sequence),
`to_fasta(path, ...)`, `to_fasta_str(two_line_mode=True)`, `to_csv(path)`.

## I/O (`seqio`)

```python
fa = read_fasta("in.fasta")                 # also accepts file-like objects
save_fasta(fa.values(), "out.fasta", mkdir=True)
save_fasta(fa.values(), "out.fasta", two_line_mode=True)  # 2-line format
```

- `read_a3m(path, is_multimer=False)` — read an A3M alignment (query first).
- `read_seqres(path, auth=True)` — sequences from PDB/mmCIF SEQRES
  (`read_mmcif_seqres`/`read_pdb_seqres` are deprecated; prefer this).
- `df2fasta(df, fasta_path, *seq_cols, id_col=None, mode="seperate", sep="")`
  — DataFrame → FASTA. `mode="seperate"` makes one entry per (id, column);
  `mode="joint"` joins the selected columns (single column forces joint).
- `temp_fasta(path, id_prefix="")` — rewrite a FASTA with numeric ids
  (`0, 1, 2...`) and return `(temp_file, id_map)`; used by the cluster
  wrappers to avoid clashing ids.
- `create_complex_seq(seq_id, *seqs, linker=":")` and
  `cross_create(seqs1, seqs2, linker=":")` — build/cross complex sequences.

## Conversion (`seqconvert`)

```python
from protools.seqconvert import translate, reverse_translate, optimize_dna, clip_a3m

aa = translate(dna_fasta, codon_table="Standard", nucleotide_type="dna")
back = reverse_translate(aa_fasta, codon_table="Standard", nucleotide_type="dna")
opt = optimize_dna(dna_fasta, species="E. coli", avoid_patterns=["BsaI"])
```

- `translate` / `reverse_translate` accept `inplace=False` by default and
  return a **new** `Fasta` with ids suffixed `_translated` / `_reverse_translated`;
  pass `inplace=True` to mutate the input object.
- `optimize_dna` uses `dnachisel` (GC 30–70% over 50-aa windows,
  `EnforceTranslation`) and asserts the translated protein is unchanged.
- `clip_a3m(a3m, start=0, end=None)` clips an A3M alignment. `start` is
  inclusive, `end` exclusive, and coordinates refer to the **query**
  sequence; members are clipped to the matching aligned columns.
- `ComplexFasta` (multi-sequence-per-record files) and `AntibodyPSSM` are in
  the same module — read the source before using.

## Cleanup (`seqfixer`, `seqfilter`)

`seqfixer` writes next to the input file — note the output names:

```python
from protools.seqfixer import deduplicated_seq, split
deduplicated_seq(Path("in.fasta"))   # -> in_deduplicated.fasta (keeps first occurrence)
split(Path("in.fasta"), 4)           # -> in/seq_0.fasta ... in/seq_3.fasta (round-robin)
```

Because both helpers write beside the input, run them on a **copy** when the
input lives in a read-protected directory like `data/` (repo rule: never
write into `data/`), then move the generated file where you need it.

`seqfilter` is registry-based:

```python
from protools.seqfilter import filter_fasta, filter_register
ok = filter_fasta("in.fasta", "standard_aa")   # built-in filter

@filter_register("my_rule")
def my_rule(item):  # item is (key, SeqRecord)
    return len(item[1].seq) > 100
```

## Annotation (`seqanno`) — antibody/TCR numbering and CDRs

```python
from protools.seqanno import numbering_seq, anno_cdr, anno_vj_gene

numbering, pident, chain_type, err = numbering_seq(
    "EVQLVESGGG...", chain="H", scheme="imgt")
cdrs = anno_cdr("EVQLVESGGG...", chain="H", scheme="imgt")
# -> {"HCDR1": "...", "HCDR1_slice": slice(...), ..., "H_percent_identity": 87.5}
vj = anno_vj_gene("EVQLVESGGG...", chain_type="H", species="human")
# -> {"v_genes": ..., "j_genes": ..., "v_pident": ..., "j_pident": ..., "species": ...}
```

- `numbering_seq(seq, chain, is_tcr=False, scheme="imgt")` returns
  `(numbering_positions, percent_identity, chain_type, error_message)`.
  A non-empty `error_message` also emits a `SeqNumberingWarning`.
- Valid schemes: antibody → `imgt | martin | kabat | aho`; TCR → `imgt` only.
  `chothia` is accepted but redirected to `martin` with a warning.
- `anno_cdr` returns one entry per region (`HCDR1`, `HCDR2`, `HCDR3`,
  `HFR1`...) with the region sequence and a `_slice` key; region labels are
  derived from `assign_cdr_labels` (frameworks become `FRn`).
- `calc_seq_identity(s1, s2, mode="global"|"local", strategy="min_aligned_length",
  max_alignments_size=0)` — Biopython pairwise alignment, score / chosen
  length strategy; cap alignments for big sequences to avoid memory blowups.
- `get_mutations(s1, s2, comparison="type"|"property"|callable, chain_id="")`
  returns strings like `["A1W", "C2D"]`; positions are **1-indexed against
  the ungapped `s1`**. Insertions/deletions are expressed relative to `s1`
  (e.g. `"-2G"` = G inserted before position 2, `"C4-"` = deletion).

## Amino-acid utilities (`aa`)

```python
from protools.aa import validate_seq, aa_equal, RESIDUE2SIDECHAIN, AA2PROPERTIES
validate_seq("ACDEFGHIKLMNPQRSTVWY")
aa_equal("D", "E", comparsion="property")   # both negative -> True
```

Constants: `RESIDUE2SIDECHAIN` / `RESIDUE2SIDECHAIN_3LETTER` (side-chain
atom names per residue), `AA_PROPERTIES` / `AA2PROPERTIES`
(hydrophobic/positive/negative/polar).

## Visualization (`seqplot`)

```python
from protools.seqplot import aa_heatmap
ax = aa_heatmap(seqs, ref_seq="EVQLVESGGG...", title="VH alignment")
```

All sequences must share the same length; `ref_seq` is outlined with a red
box. `display_intervals` (an `Intervals` object) must be **one-based**
(`zero_based=False`). Returns a matplotlib `Axes`.

## Pitfalls

- `translate`/`optimize_dna` return copies unless `inplace=True`; ids get
  `_translated`/`_optimized` suffixes. Don't assume the input was modified.
- `Fasta` values are always `SeqRecord`; when you assign a `str` it is
  wrapped automatically.
- A3M clipping coordinates are query-based and column-aware; clip members by
  aligned position, not raw string index.
- Antibody numbering (IMGT/Kabat/Martin/Aho) must come from `antpack` via
  this module — never hard-code or invent numbering in new code or tests.
