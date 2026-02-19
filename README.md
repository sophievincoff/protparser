# protparser — a user-friendly package for parsing and comprehending RCSB PDB files

`protparser` is a Python toolkit for downloading and parsing RCSB PDB entries (mmCIF + FASTA) and producing structured, analysis-ready pandas DataFrames at three levels:

- **Chain inventory** — one row per chain (mmCIF's `label_asym_id`)
- **Entity inventory** — one row per entity (mmCIF's `entity_id`)
- **Entry inventory** — one row per PDB structure

It is designed for:

- Structural dataset curation
- Binary interaction detection
- PTM/noncanonical residue auditing
- Sequence validation and ambiguity detection
- Entity-level and structure-level summarization
- Automated reporting

This README documents the public API and explains **every major field**, with explicit references to the originating mmCIF fields (shown in parentheses).

---

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

## Batch Analysis

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

# Important Conceptual Distinction

## Invalid Amino Acids vs Non-Canonical Residues

These are **different biological concepts** in this codebase.

### 1️⃣ Invalid / Ambiguous Canonical Tokens

Derived from canonical sequence
(`_entity_poly.pdbx_seq_one_letter_code_can`)

Examples:

* `X`, `B`, `J`, `Z`, `*`

Columns:

* `is_valid_aa`
* `invalid_letters`

These represent **ambiguous or unknown sequence tokens**, not structural chemical modifications.

---

### 2️⃣ Non-Canonical Residue Names

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

# Public API Overview

## `build_inventory(include_rcsb_rest=True)`

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

Optional REST enrichment:

* `https://data.rcsb.org/rest/v1/core/polymer_entity/...`

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
