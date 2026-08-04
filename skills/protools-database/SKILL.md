---
name: protools-database
description: "Fetch and parse protein data from UniProt and WIPO patent sequence listings with the protools library (`protools4py`): download SwissProt records by accession, extract feature tables and extracellular chains, fetch features to TSV, and parse WIPO ST.26-style sequence lists into structured records. Trigger whenever the user needs UniProt/SwissProt records or features, accession-based protein metadata, patent sequence lists, or sequence-list XML/text parsing in the protools repo, even if they do not name the module (database/uniprot, database/patent)."
---

# protools-database

External protein data access in the `protools` package:

- `protools.database.uniprot` — SwissProt records and feature extraction
  over the network (Bio.ExPASy + Bio.SwissProt).
- `protools.database.patent` — WIPO patent sequence-list parsing
  (ST.26-style tagged elements) into a structured object.

## How to work in this repo

- Run code with `uv run python`, tests with `uv run pytest test/database/`.
- Look up exact signatures with CodeGraph or the module source before coding.
- UniProt access is network I/O: 4xx/5xx and rate-limit failures are
  infrastructure problems, not code bugs. The module already decorates
  fetches with `@max_retry(max=5)` and `@catch_error`.
- Never fabricate accession IDs, feature coordinates, or patent data —
  these are real-world facts. When a value is not available from the source,
  say so rather than guessing.
- Real data over mocks: integration tests hit the live endpoints
  (`test/database/test_uniprot.py`); use real accessions in new tests.

## UniProt (`protools.database.uniprot`)

```python
from protools.database.uniprot import (
    get_uniprot_record, extract_features, extract_extracellular_chains,
    fetch_uniprot_features,
)

record = get_uniprot_record("P01308")           # Bio.SwissProt.Record
features = extract_features(record)             # generator of dicts
chains_df = extract_extracellular_chains(record)
fetch_uniprot_features("P01308", Path("out/"))  # writes out/P01308_features.txt
```

- `get_uniprot_record(accession_id)` — SwissProt record via ExPASy
  (`get_sprot_raw`); isoform ids like `P01308-2` are supported by the
  fetcher (split on `-`).
- `extract_features(record)` — yields one dict per feature with at least
  `type`, `pos_range`, `subseq`, plus all qualifiers merged in; non-exact
  positions produce an empty `subseq`.
- `extract_extracellular_chains(record)` — DataFrame of extracellular chain
  features (see source for the exact column set).
- `fetch_uniprot_features(accession_id, output_path, skip_exist=False)` —
  writes a TSV `{accession_id}_features.txt` into `output_path`
  (directory is created); `skip_exist=True` reuses an existing file. The
  module's `__main__` supports batch fetching with `n_jobs` multiprocessing.

## Patent sequence lists (`protools.database.patent`)

```python
from protools.database.patent import PatentSequenceList, SequenceListElement

psl = PatentSequenceList.from_file("sequence_list.txt")
for seq_id, seq_record in psl.items():   # OrderedDict-like
    ...
```

- `SequenceListElement` is an enum of WIPO ST.26 element tags — e.g.
  `110` applicant names, `120` invention title, `160` number of seq ids,
  `210` seq id, `211` seq length, `212` moltype, `213` organism, `220`
  feature, `221` feature name, `222` feature location, `400` sequence.
  Standard source: https://www.wipo.int/export/sites/www/standards/en/pdf/03-25-01.pdf
- `PatentSequenceList.from_file(file_path)` parses a tagged text/list file
  into an `OrderedDict` subclass keyed by sequence id; unknown tags map to
  `SequenceListElement.UNDEFINED`. Read the source before extending the
  parser — the format has mandatory vs optional elements
  (`SequenceListElementTuple(id, is_mandatory)`).
  Note: the parser prints each recognized element name to stdout as a side
  effect — don't treat that as an error, and don't duplicate it as your own
  logging.

Two parser quirks matter in practice: `<400>` blocks expect three-letter
amino-acid codes (`Ala Cys Asp ...`) for PRT sequences, and each
sequence-data line consumes the following line; a record's sequence is
materialized only when the next element tag arrives, so the final `<400>`
block needs a trailing tag (e.g. `<160>`) after it to get a non-empty
sequence.

## Pitfalls

- UniProt records are live data: accessions can be merged, split, or
  withdrawn. Re-verify against UniProt when a test or downstream consumer
  depends on a specific record's features.
- `extract_features` yields *all* feature types (chains, regions, sites,
  variants...); filter by `type` when a consumer expects only one kind.
- Patent sequence lists are WIPO-standardized; validate against the ST.26
  spec rather than assuming a generic XML parser applies.
- Keep network-dependent tests marked/isolated as integration tests per
  the repo conventions (rate limits are expected occasionally).
