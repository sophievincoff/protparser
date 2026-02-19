# protparser — a user-friendly package for parsing and comprehending RCSB PDB files

`protparser` is a Python toolkit for downloading and parsing RCSB PDB entries (mmCIF + FASTA). It produces structured, analysis-ready pandas DataFrames at three levels, and aggregates their information into a readable report.

- **Entry inventory** — one row per PDB structure; high-level info
- **Chain inventory** — one row per chain
- **Entity inventory** — one row per entity
- **Report** - one streamlined report combining chain, entity, and entry information.

It is designed for:

- Structural dataset curation
- Binary interaction detection
- PTM/noncanonical residue auditing
- Sequence validation and ambiguity detection
- Entity-level and structure-level summarization
- Automated reporting

This README documents the public API and explains **every major field**, with explicit references to the originating mmCIF fields.

---

# Preview of capabilities

Here, we preview some capabilities of `protparser` using RCSB entry 1AYC as an example. 

## Entry inventory

| pdb_id   | pdb_title                                                                                                | protein_chains   | protein_chains_auth   | protein_names                                                              | protein_uniprot_ids   | protein_species_names   | protein_species_taxids   | dna_names   | rna_names   | nonpolymer_names   | entry_has_protein_modification   | entry_has_protein_modification_raw   | n_chains_total   | n_protein_chains   | has_any_chain_ptm_evidence   | ptm_chains   | ptm_chains_auth   | has_any_chain_noncanonical_evidence   | noncanonical_chains   | is_binary_interaction   | binary_n_unique_protein_seqs   | binary_stoichiometry   | reduced_binary_stoichiometry   | dna_chains   | dna_chains_auth   | rna_chains   | rna_chains_auth   | water_chains   | water_chains_auth   | nonpolymer_other_chains   | nonpolymer_other_chains_auth   |
|----------|----------------------------------------------------------------------------------------------------------|------------------|-----------------------|----------------------------------------------------------------------------|-----------------------|-------------------------|--------------------------|-------------|-------------|--------------------|----------------------------------|--------------------------------------|------------------|--------------------|------------------------------|--------------|-------------------|---------------------------------------|-----------------------|-------------------------|--------------------------------|------------------------|--------------------------------|--------------|-------------------|--------------|-------------------|----------------|---------------------|---------------------------|--------------------------------|
| 1ayc     | CRYSTAL STRUCTURES OF PEPTIDE COMPLEXES OF THE AMINO-TERMINAL SH2 DOMAIN OF THE SYP TYROSINE PHOSPHATASE | A,B              | A,P                   | PEPTIDE PDGFR-740|PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN) | P05622|P35235         | Mus musculus            | 10090                    | nan         | nan         | water              | True                             | True                                 | 4                | 2                  | True                         | B            | P                 | False                                 | nan                   | True                    | 2                              | 1:1                    | 1:1                            | nan          | nan               | nan          | nan               | C,D            | C,D                 | nan                       | nan                            |


## Chain inventory

| pdb_id   | chain_id   | chain_id_auth   | entity_id   | entity_type   | entity_description                                       | entity_src_method                                           | entity_number_of_molecules   | entity_mutation   | entity_ec   | entity_fragment   | entity_details   | nonpoly_name   | chem_comp_name   | seq_can                                                                                               | seq_len   | is_valid_aa   | invalid_letters   | entity_poly_type   | has_ptm   | has_ptm_reason                             | ptm_detect_sources                   | ptm_resnames   | ptm_mf_chain_hit   | has_noncanonical_resname   | noncanonical_resnames   | ptm_categories             | ptm_types       | ptm_comp_ids   | ptm_n_features   | seq_resolved_masked                                                                                   | seq_resolved_only                                                                                    | n_resolved   | frac_resolved      | seq_can_annot                                                                                         | seq_resolved_masked_annot                                                                             | ptm_pos_can   | noncanon_pos_can   | ptm_pos_resolved_only   | noncanon_pos_resolved_only   | chain_id_auth_rest   | uniprot_ids   | mutation_pos_can   | rcsb_conflict_count   | rcsb_deletion_count   | rcsb_insertion_count   | rcsb_mutation_count   | rcsb_non_std_monomer_count   | rcsb_sample_sequence_length   | chain_id_auth_cif   |
|----------|------------|-----------------|-------------|---------------|----------------------------------------------------------|-------------------------------------------------------------|------------------------------|-------------------|-------------|-------------------|------------------|----------------|------------------|-------------------------------------------------------------------------------------------------------|-----------|---------------|-------------------|--------------------|-----------|--------------------------------------------|--------------------------------------|----------------|--------------------|----------------------------|-------------------------|----------------------------|-----------------|----------------|------------------|-------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|--------------|--------------------|-------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|---------------|--------------------|-------------------------|------------------------------|----------------------|---------------|--------------------|-----------------------|-----------------------|------------------------|-----------------------|------------------------------|-------------------------------|---------------------|
| 1ayc     | A          | A               | 1           | polymer       | PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN) | man (entity isolated from a genetically manipulated source) | 1                            | nan               | 3.1.3.48    | nan               | nan              | nan            | nan              | MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | 101       | True          | nan               | polypeptide(L)     | False     | nan                                        | nan                                  | nan            | nan                | False                      | nan                     | nan                        | nan             | nan            | 0                | -RRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | RRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | 100          | 0.9900990099009901 | MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | -RRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | nan           | nan                | nan                     | nan                          | A                    | P35235        | nan                | 0                     | 0                     | 0                      | 0                     | 0                            | 101                           | A                   |
| 1ayc     | B          | P               | 2           | polymer       | PEPTIDE PDGFR-740                                        | man (entity isolated from a genetically manipulated source) | 1                            | nan               | nan         | nan               | nan              | nan            | nan              | DGGYMDMSKGS                                                                                           | 11        | True          | nan               | polypeptide(L)     | True      | chain listed in _pdbx_modification_feature | pdbx_modification_feature,entry_flag | nan            | True               | False                      | nan                     | Named protein modification | Phosphorylation | PTR            | 1                | -GG-MDMS---                                                                                           | GGMDMS                                                                                               | 6            | 0.5454545454545454 | DGG<Y+PTR>MDMSKGS                                                                                     | -GG-MDMS---                                                                                           | 4             | nan                | nan                     | nan                          | P                    | P05622        | nan                | 0                     | 0                     | 0                      | 0                     | 1                            | 11                            | B                   |
| 1ayc     | C          | C               | 3           | water         | water                                                    | nat (entity isolated from a natural source)                 | 63                           | nan               | nan         | nan               | nan              | water          | WATER            | nan                                                                                                   | 0         | False         | no_sequence       | nan                | False     | nan                                        | nan                                  | nan            | nan                | False                      | nan                     | nan                        | nan             | nan            | 0                | nan                                                                                                   | nan                                                                                                  | 0            | 0.0                | nan                                                                                                   | nan                                                                                                   | nan           | nan                | nan                     | nan                          | nan                  | nan           | nan                | nan                   | nan                   | nan                    | nan                   | nan                          | nan                           | C                   |
| 1ayc     | D          | D               | 3           | water         | water                                                    | nat (entity isolated from a natural source)                 | 63                           | nan               | nan         | nan               | nan              | water          | WATER            | nan                                                                                                   | 0         | False         | no_sequence       | nan                | False     | nan                                        | nan                                  | nan            | nan                | False                      | nan                     | nan                        | nan             | nan            | 0                | nan                                                                                                   | nan                                                                                                  | 0            | 0.0                | nan                                                                                                   | nan                                                                                                   | nan           | nan                | nan                     | nan                          | nan                  | nan           | nan                | nan                   | nan                   | nan                    | nan                   | nan                          | nan                           | D                   |

## Entity inventory

| pdb_id   | entity_id   | chain_display   | chain_ids   | chain_id_auths   | entity_type   | entity_description                                       | entity_src_method                                           | entity_number_of_molecules   | entity_mutation   | entity_ec   | entity_fragment   | entity_details   | nonpoly_name   | chem_comp_name   | seq_can                                                                                               | seq_len   | is_valid_aa   | invalid_letters   | entity_poly_type   | has_ptm   | has_ptm_reason                             | ptm_detect_sources                   | ptm_resnames   | ptm_mf_chain_hit   | has_noncanonical_resname   | noncanonical_resnames   | ptm_categories             | ptm_types       | ptm_comp_ids   | ptm_n_features   | seq_resolved_masked                                                                                   | seq_resolved_only                                                                                    | n_resolved   | frac_resolved      | seq_can_annot                                                                                         | seq_resolved_masked_annot                                                                             | ptm_pos_can   | noncanon_pos_can   | ptm_pos_resolved_only   | noncanon_pos_resolved_only   | chain_id_auth_rest   | uniprot_ids   | mutation_pos_can   | rcsb_conflict_count   | rcsb_deletion_count   | rcsb_insertion_count   | rcsb_mutation_count   | rcsb_non_std_monomer_count   | rcsb_sample_sequence_length   | chain_id_auth_cif   | entity_varying_columns   |
|----------|-------------|-----------------|-------------|------------------|---------------|----------------------------------------------------------|-------------------------------------------------------------|------------------------------|-------------------|-------------|-------------------|------------------|----------------|------------------|-------------------------------------------------------------------------------------------------------|-----------|---------------|-------------------|--------------------|-----------|--------------------------------------------|--------------------------------------|----------------|--------------------|----------------------------|-------------------------|----------------------------|-----------------|----------------|------------------|-------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|--------------|--------------------|-------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|---------------|--------------------|-------------------------|------------------------------|----------------------|---------------|--------------------|-----------------------|-----------------------|------------------------|-----------------------|------------------------------|-------------------------------|---------------------|--------------------------|
| 1ayc     | 1           | A               | A           | A                | polymer       | PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN) | man (entity isolated from a genetically manipulated source) | 1                            | nan               | 3.1.3.48    | nan               | nan              | nan            | nan              | MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | 101       | True          | nan               | polypeptide(L)     | False     | nan                                        | nan                                  | nan            | nan                | False                      | nan                     | nan                        | nan             | nan            | 0                | -RRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | RRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | 100          | 0.9900990099009901 | MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | -RRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN | nan           | nan                | nan                     | nan                          | A                    | P35235        | nan                | 0                     | 0                     | 0                      | 0                     | 0                            | 101                           | A                   | nan                      |
| 1ayc     | 2           | B[auth: P]      | B           | P                | polymer       | PEPTIDE PDGFR-740                                        | man (entity isolated from a genetically manipulated source) | 1                            | nan               | nan         | nan               | nan              | nan            | nan              | DGGYMDMSKGS                                                                                           | 11        | True          | nan               | polypeptide(L)     | True      | chain listed in _pdbx_modification_feature | pdbx_modification_feature,entry_flag | nan            | True               | False                      | nan                     | Named protein modification | Phosphorylation | PTR            | 1                | -GG-MDMS---                                                                                           | GGMDMS                                                                                               | 6            | 0.5454545454545454 | DGG<Y+PTR>MDMSKGS                                                                                     | -GG-MDMS---                                                                                           | 4             | nan                | nan                     | nan                          | P                    | P05622        | nan                | 0                     | 0                     | 0                      | 0                     | 1                            | 11                            | B                   | nan                      |
| 1ayc     | 3           | C,D             | C|D         | C|D              | water         | water                                                    | nat (entity isolated from a natural source)                 | 63                           | nan               | nan         | nan               | nan              | water          | WATER            | nan                                                                                                   | 0         | False         | no_sequence       | nan                | False     | nan                                        | nan                                  | nan            | nan                | False                      | nan                     | nan                        | nan             | nan            | 0                | nan                                                                                                   | nan                                                                                                  | 0            | 0.0                | nan                                                                                                   | nan                                                                                                   | nan           | nan                | nan                     | nan                          | nan                  | nan           | nan                | nan                   | nan                   | nan                    | nan                   | nan                          | nan                           | C|D                 | chain_id_auth_cif        |

## FASTA inventory

| pdb_id   | record_id   | chains_pdb   | chains_auth   | pdb_name                                                 | species      |   species_id | sequence                                                                                              |   seq_len |
|:---------|:------------|:-------------|:--------------|:---------------------------------------------------------|:-------------|-------------:|:------------------------------------------------------------------------------------------------------|----------:|
| 1ayc     | 1AYC_1      | A            | A             | PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN) | Mus musculus |        10090 | MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN |       101 |
| 1ayc     | 1AYC_2      | B            | P             | PEPTIDE PDGFR-740                                        | Mus musculus |        10090 | DGGYMDMSKGS                                                                                           |        11 |


## Report

```
1AYC
- Entry flag has_protein_modification: True
- Chains: 4 (protein chains: 2)
- Binary interaction (protein seqs): True (n_unique=2, stoich=1:1)
- Entities:
  - entity 1 | polymer | PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN)
    - Chains: A
    - Varying columns: none
  - entity 2 | polymer | PEPTIDE PDGFR-740
    - Chains: B[auth: P]
    - Varying columns: none
  - entity 3 | water | water
    - Chains: C,D
    - Varying columns: chain_id_auth_cif
- Chains:
  - Chain A | entity 1 | polypeptide(L)
    - Lengths: canonical=101, resolved=100 (frac=0.990) | UniProt=P35235
    - Invalids: none
    - Noncanonicals: none detected
    - PTMs: none detected
  - Chain B[auth: P] | entity 2 | polypeptide(L)
    - Lengths: canonical=11, resolved=6 (frac=0.545) | UniProt=P05622
    - Invalids: none
    - Noncanonicals: none detected
    - PTMs: pos=4; comp_ids=PTR; types=Phosphorylation; categories=Named protein modification
  - Chain C | entity 3 | non-polymer
  - Chain D | entity 3 | non-polymer

```

# Installation

Requirements:

- Python 3.10+
- `gemmi`
- `pandas`
- `numpy`
- `requests`

```
pip install gemmi pandas numpy requests
````

If developing locally:

```
pip install -e .
```

---

# Core Class: `RCSBStructure`

```python
from protparser.rcsb import RCSBStructure

r = RCSBStructure("2LTV", cache_dir="rcsb_cache")

chain_df  = r.build_inventory(include_rcsb_rest=True)
entity_df = r.build_entity_inventory(include_rcsb_rest=True)
entry_df  = r.build_entry_inventory()

print(r.generate_report())
```

### Constructor

```python
RCSBStructure(
    pdb_id,
    cache_dir=None,
    cif_path=None,
    download=True,
    valid_extra_tokens=None
)
```

* `pdb_id`: PDB identifier (case-insensitive)
* `cache_dir`: where mmCIF/FASTA files are cached
* `cif_path`: explicit mmCIF path (skips download)
* `download`: whether to auto-download mmCIF
* `valid_extra_tokens`: extend allowed canonical AA alphabet

---

# Typical Workflows

## Analyze a Single PDB

```python
from pathlib import Path
from protparser.rcsb import RCSBStructure

out = Path("outputs/2LTV")
out.mkdir(parents=True, exist_ok=True)

r = RCSBStructure("2LTV", cache_dir="rcsb_cache")

chain_df  = r.build_inventory(include_rcsb_rest=True)
entity_df = r.build_entity_inventory(include_rcsb_rest=True)
entry_df  = r.build_entry_inventory()

chain_df.to_csv(out / "chain_inventory.csv", index=False)
entity_df.to_csv(out / "entity_inventory.csv", index=False)
entry_df.to_csv(out / "entry_inventory.csv", index=False)

r.generate_report(output="file", out_path=out / "report.txt")
```

---

## Batch Analysis (Loop)

This workflow is most useful if you have a small list of PDBs, and cannot take advantage of multiprocessing on your computer system. 

```python
import pandas as pd
from protparser.rcsb import RCSBStructure

pdb_ids = ["2LTV", "4OO8", "1A2B"]

all_entries = []

for pid in pdb_ids:
    r = RCSBStructure(pid, cache_dir="rcsb_cache")
    r.build_inventory(include_rcsb_rest=True)
    entry_df = r.build_entry_inventory()
    all_entries.append(entry_df)

pd.concat(all_entries).to_csv("entry_inventory_all.csv", index=False)
```

---

## Batch Analysis (Parallel)

For analyzing many PDBs at once, use the functions in `parallel_inventory.py`. This builds inventories **in parallel**, using multiprocessing, and returns a single results object with combined DataFrames + per-PDB inventories.

### 1) Build inventories in parallel

```python
from protparser.parallel_inventory import build_inventories

pdb_ids = ["2LTV", "4OO8", "1A2B"]

res = build_inventories(
    pdb_ids,
    cache_dir="downloads",
    include_rcsb_rest=True,
    parallel=True,
    executor="thread",   # "thread" or "process" (if supported)
)
```

### 2) Access the combined DataFrames

`build_inventories(...)` returns an object containing **group-level** (combined) inventories across all PDBs:

```python
chain_all  = res.chain_inventory_df
entity_all = res.entity_inventory_df
entry_all  = res.entry_inventory_df

# Example: save combined inventories
chain_all.to_csv("chain_inventory_all.csv", index=False)
entity_all.to_csv("entity_inventory_all.csv", index=False)
entry_all.to_csv("entry_inventory_all.csv", index=False)
```

---

## Writing inventory reports

Use `write_inventory_reports(res, ...)` to produce nice, human-readable inventory summaries. You can run it in:

* **Individual PDB mode** (one report per PDB)
* **Group mode** (one report summarizing the whole batch)

And you can either:

* **Return a string** (for printing / logging / notebooks)
* **Save to file(s)**

### A) Group report (single report for the whole batch)

**Return as a string:**

```python
from protparser.parallel_inventory import write_inventory_reports

txt = write_inventory_reports(
    res,
    mode="group",          # group summary across all pdb_ids
    output="string",       # return a string instead of writing files
)
print(txt)
```

**Save to a file:**

```python
write_inventory_reports(
    res,
    mode="group",
    output="file",
    out_path="inventory_group_report.txt",
)
```

---

### B) Individual PDB reports (one report per PDB)

**Return strings (e.g., a dict of `{pdb_id: report_text}`):**

```python
reports = write_inventory_reports(
    res,
    mode="individual",
    output="string",
)

print(reports["4OO8"])   # view a single PDB’s report
```

**Save one file per PDB (e.g., to a folder):**

```python
write_inventory_reports(
    res,
    mode="individual",
    output="file",
    out_dir="inventory_reports/",
    filename_template="{pdb_id}_inventory.txt",
)
```

---

# Important Conceptual Distinction

## Invalid Amino Acids vs Non-Canonical Residues

These are **different biological concepts** in this codebase.

### 1. Invalid / Ambiguous Canonical Tokens

Derived from canonical sequence
(`_entity_poly.pdbx_seq_one_letter_code_can`)

Examples:

* `X`, `B`, `J`, `Z`, `*`

Columns:

* `is_valid_aa`
* `invalid_letters`

These represent **ambiguous or unknown sequence tokens**, not structural chemical modifications.

---

### 2. Non-Canonical Residue Names

Derived from structure model
(`_atom_site.label_comp_id`, ATOM records)

Examples:

* `MSE`, `SEP`, `TPO`, `PTR`

Columns:

* `has_noncanonical_resname`
* `noncanonical_resnames`
* `noncanon_pos_can`

These represent **actual non-standard monomers present in the structural model.**

---

# Public API Overview: RCSBStructure Object

## `build_inventory()`

Returns **chain inventory** (`inventory_df`).

Source mmCIF categories:

* `_struct_asym`
* `_entity`
* `_entity_poly`
* `_atom_site`
* `_pdbx_modification_feature`
* `_pdbx_entry_details`
* `_pdbx_entity_nonpoly`
* `_chem_comp`

REST enrichment: uses the public REST API of the RCSB PDB to collect additional information. You may toggle off the REST API query by calling `build_inventory(include_rcsb_rest=False)`. The URLS used are shown below. 

* Entry-level information: `https://data.rcsb.org/rest/v1/core/entry/{PDB}`
* Polymer entity information: `https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB}/{polymer_entity}`
* FASTA: `https://www.rcsb.org/fasta/entry/{PDB}` 

---

## `build_entity_inventory()`

Aggregates chain inventory rows per entity (`entity_id`).

---

## `build_entry_inventory()`

Produces one row per PDB entry.

Includes structure-level rollups and metadata.

---

## FASTA Utilities

| Method                    | Description                                                    |
| ------------------------- | -------------------------------------------------------------- |
| `load_fasta()`            | Download/load FASTA (`https://www.rcsb.org/fasta/entry/{PDB}`) |
| `write_fasta(path)`       | Save FASTA                                                     |
| `build_fasta_inventory()` | Parse FASTA into DataFrame                                     |

For example, the RCSB FASTA for 1AYC: 

```
>1AYC_1|Chain A|PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN)|Mus musculus (10090)
MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN
>1AYC_2|Chain B[auth P]|PEPTIDE PDGFR-740|Mus musculus (10090)
DGGYMDMSKGS
```

and the corresponding DataFrame from `build_fasta_inventory()`:

|    | pdb_id   | record_id   | chains_pdb   | chains_auth   | pdb_name                                                 | species      |   species_id | sequence                                                                                              |   seq_len |
|---:|:---------|:------------|:-------------|:--------------|:---------------------------------------------------------|:-------------|-------------:|:------------------------------------------------------------------------------------------------------|----------:|
|  0 | 1ayc     | 1AYC_1      | A            | A             | PROTEIN-TYROSINE PHOSPHATASE SYP (N-TERMINAL SH2 DOMAIN) | Mus musculus |        10090 | MRRWFHPNITGVEAENLLLTRGVDGSFLARPSKSNPGDFTLSVRRNGAVTHIKIQNTGDYYDLYGGEKFATLAELVQYYMEHHGQLKEKNGDVIELKYPLN |       101 |
|  1 | 1ayc     | 1AYC_2      | B            | P             | PEPTIDE PDGFR-740                                        | Mus musculus |        10090 | DGGYMDMSKGS                                                                                           |        11 |

---

# Chain Inventory (One Row per Chain)

### Core Identity

| Column          | Meaning          | mmCIF Source                       |
| --------------- | ---------------- | ---------------------------------- |
| `pdb_id`        | Lowercase PDB ID | derived                            |
| `chain_id`      | Label chain ID   | (`_struct_asym.id`)                |
| `chain_id_auth` | Author chain ID  | (`_struct_asym.pdbx_auth_asym_id`) |
| `entity_id`     | Entity reference | (`_struct_asym.entity_id`)         |

---

### Entity Metadata

| Column                       | mmCIF Source                         |
| ---------------------------- | ------------------------------------ |
| `entity_type`                | (`_entity.type`)                     |
| `entity_description`         | (`_entity.pdbx_description`)         |
| `entity_src_method`          | (`_entity.src_method`)               |
| `entity_number_of_molecules` | (`_entity.pdbx_number_of_molecules`) |
| `entity_mutation`            | (`_entity.pdbx_mutation`)            |
| `entity_ec`                  | (`_entity.pdbx_ec`)                  |
| `entity_fragment`            | (`_entity.pdbx_fragment`)            |
| `entity_details`             | (`_entity.details`)                  |
| `nonpoly_name`               | (`_pdbx_entity_nonpoly.name`)        |
| `chem_comp_name`             | (`_chem_comp.name`)                  |

---

### Polymer Type

| Column             | Source                |
| ------------------ | --------------------- |
| `entity_poly_type` | (`_entity_poly.type`) |

Common values:

* `polypeptide(L)`
* `polydeoxyribonucleotide`
* `polyribonucleotide`

---

### Canonical Sequence

| Column    | Source                                        |
| --------- | --------------------------------------------- |
| `seq_can` | (`_entity_poly.pdbx_seq_one_letter_code_can`) |
| `seq_len` | derived                                       |

---

### Sequence Validation

| Column            | Meaning                             |
| ----------------- | ----------------------------------- |
| `is_valid_aa`     | All tokens in allowed alphabet      |
| `invalid_letters` | Sorted list of offending characters |

⚠ Only relevant for `polypeptide(L)`.

---

### Resolved Structure Coverage

Derived from:

* (`_atom_site.group_PDB == "ATOM"`)
* (`_atom_site.label_asym_id`)
* (`_atom_site.label_seq_id`)

Columns:

* `seq_resolved_masked`
* `seq_resolved_only`
* `n_resolved`
* `frac_resolved`

---

### PTM / Modification Evidence

Sources:

* (`_pdbx_modification_feature`)
* (`_pdbx_entry_details.has_protein_modification`)
* residue-name heuristics from (`_atom_site.label_comp_id`)

Columns:

| Column               | Meaning                                      |
| -------------------- | -------------------------------------------- |
| `has_ptm`            | Chain-local PTM evidence                     |
| `has_ptm_reason`     | Human-readable reason                        |
| `ptm_detect_sources` | detection mechanisms                         |
| `ptm_resnames`       | heuristic residue names                      |
| `ptm_categories`     | (`_pdbx_modification_feature.category`)      |
| `ptm_types`          | (`_pdbx_modification_feature.type`)          |
| `ptm_comp_ids`       | (`_pdbx_modification_feature.label_comp_id`) |
| `ptm_n_features`     | Count of feature rows                        |

---

### Non-Canonical Residues (Structure-Level)

| Column                       | Source                                    |
| ---------------------------- | ----------------------------------------- |
| `has_noncanonical_resname`   | derived from (`_atom_site.label_comp_id`) |
| `noncanonical_resnames`      | derived                                   |
| `noncanon_pos_can`           | derived                                   |
| `noncanon_pos_resolved_only` | derived                                   |

---

### REST Enrichment (Optional)

Source:
`https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB}/{entity}`

Includes:

* `uniprot_ids`
* `uniprot_alignments_json`
* `pfam_json`
* `disordered_regions_json`
* `disordered_binding_sites_json`
* mutation positions
* conflict/deletion/insertion counts

---

# Entity Inventory (One Row per Entity)

Aggregates chain rows per (`_struct_asym.entity_id`).

Key columns:

* `entity_id`
* `chain_ids`
* `chain_id_auths`
* `chain_display`
* `entity_varying_columns`

Any column with inconsistent values across chains is pipe-joined.

---

# Entry Inventory (One Row per PDB)

### Core Structure Metadata

| Column                           | mmCIF Source                                     |
| -------------------------------- | ------------------------------------------------ |
| `pdb_id`                         | derived                                          |
| `pdb_title`                      | (`_struct.title`)                                |
| `entry_has_protein_modification` | (`_pdbx_entry_details.has_protein_modification`) |

---

### Chain Counts

* `n_chains_total`
* `n_protein_chains`

---

### Chain Category Lists

Derived from `entity_poly_type` and `entity_description`:

* `protein_chains`
* `dna_chains`
* `rna_chains`
* `water_chains`
* `nonpolymer_other_chains`

(Auth variants also included.)

---

### Names

Derived from (`_entity.pdbx_description`) or nonpoly fields:

* `protein_names`
* `dna_names`
* `rna_names`
* `nonpolymer_names`

---

### Species Metadata (Proteins)

From:

* `_entity_src_nat`
* `_entity_src_gen`
* `_entity_src_syn`

Columns:

* `protein_species_names`
* `protein_species_taxids`

---

### UniProt Rollup

* `protein_uniprot_ids`

(From REST polymer_entity alignment.)

---

### PTM Rollups

* `has_any_chain_ptm_evidence`
* `ptm_chains`
* `ptm_chains_auth`

---

### Noncanonical Rollups

* `has_any_chain_noncanonical_evidence`
* `noncanonical_chains`

---

### Invalid Amino Acid Rollups

(Only `polypeptide(L)` chains)

* `has_any_chain_invalid_aa`
* `invalid_aa_chains`
* `invalid_aa_chains_auth`

---

### Binary Interaction Summary

Derived from canonical protein sequences (`seq_can`):

* `binary_n_unique_protein_seqs`
* `is_binary_interaction`
* `binary_stoichiometry`

Definition:

> Binary interaction = exactly two unique protein sequences, regardless of chain multiplicity.

---

# Reporting

`generate_report()` produces a structured summary including:

* Entry-level modification flag
* Chain counts
* Binary interaction summary
* Entity summary
* Per-chain sequence, PTM, invalid, and noncanonical status

---

# Design Philosophy

* Explicit mmCIF field provenance
* Clear separation of entry/entity/chain levels
* Strict separation of:

  * canonical sequence validation
  * structural residue-name analysis
  * entry-level metadata flags
* Robust to partial mmCIF schemas
* CSV-friendly

---

# Future Directions

Potential extensions:

* Symmetry-aware assembly stoichiometry
* Resolution and experimental method extraction (`_refine`, `_exptl`)
* Ligand chemical property enrichment
* Inter-chain contact summaries

---

# Author

Sophia Vincoff

Bioengineering PhD Candidate in Dr. Pranam Chatterjee's lab

University of Pennsylvania, School of Engineering and Applied Science

---
