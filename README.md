# protools
Personal protein processing tool collections.

## Installment

Necessary requirement is listed in [requirements.txt](./requirements.txt).
Optional requirement is listed below:
```bash
conda install -c schrodinger pymol
```
Also, pypi release will be updated every github release, just install it by `uv pip`:
```bash
uv pip install protools4py
```

## Usage

### Sequence Modules
- protools.seqio
- protools.seqanno
- protools.seqconvert
- protools.seqfixer
- protools.aa

### Structure Modules
- protools.pdbio
- protools.pdbanno
- protools.pdbconvert
- protools.pdbfixer
- protools.pdbplot
- protools.dock

### ML Dataset Modules
- protools.dataset.split

### Cluster Algorithm Modules
- protools.cluster.cdhit
- protools.cluster.mmseqs


### Protein Database Modules
- protools.database.uniprot
- protools.database.patent

### TODO

- [ ] Add documents for each module.

