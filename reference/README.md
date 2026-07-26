# Reference: Stereochemistry Featurization in RDKit

This folder contains the empirical analysis and design rationale for stereochemistry featurization decisions in this project. The key conclusion is that chemprop's atom chirality encoding (`@`/`@@` → CW/CCW) does not map to absolute stereochemistry (R/S), and a CIP-based encoding should be added.

## Quick Start

**Read this first**: [`stereochemistry_quick_reference.md`](stereochemistry_quick_reference.md) — skimmable summary of all RDKit stereo functions, what they return, and what action to take for each.

## Reports

| File | Contents |
|---|---|
| [`stereochemistry_quick_reference.md`](stereochemistry_quick_reference.md) | One-page summary: all RDKit functions, enum tables, action items |
| [`atom_chirality_analysis_report.md`](atom_chirality_analysis_report.md) | Full atom chirality analysis: @/@@ vs R/S, why they don't map, dataset results |
| [`bond_stereo_analysis_report.md`](bond_stereo_analysis_report.md) | Full bond stereo analysis: CIS/TRANS vs E/Z, why current encoding is adequate |

## Key Findings

1. **Atom chirality (`GetChiralTag`)**: CW/CCW does NOT map to R/S (~50% mismatch even with canonical SMILES). Must use `rdCIPLabeler.AssignCIPLabels(mol)` to get absolute R/S labels.
2. **Bond stereochemistry (`GetStereo`)**: CIS/TRANS maps to E/Z 99.93% of the time. Current chemprop encoding is adequate — no change needed.
3. **Implementation**: Added optional `cip_codes` parameter to chemprop's `MultiHotAtomFeaturizer` (backward compatible, opt-in).

## Data

Analysis was run across 19 MoleculeNet andChEMBL benchmark datasets (1.23M molecules total) from chemprop v1 https://github.com/kspieks/chemprop/blob/barrier_prediction/data.tar.gz. See the full reports for dataset-level breakdowns.

## Reproducibility

Scripts and outputs are in subdirectories:

```
scripts/
├── analyze_atom_chirality_datasets.py   # Atom chirality analysis (all 19 datasets)
└── analyze_bond_stereo_datasets.py      # Bond stereo analysis (all 19 datasets)

results/
├── atom_chirality_analysis_results.json # Structured output (counts, examples, mappings)
├── bond_stereo_analysis_results.json
├── atom_chirality_analysis.log          # Full printed output
└── bond_stereo_analysis.log
```

To re-run:
```bash
conda activate chemprop
python scripts/analyze_atom_chirality_datasets.py 2>&1 | tee results/atom_chirality_analysis.log
python scripts/analyze_bond_stereo_datasets.py 2>&1 | tee results/bond_stereo_analysis.log
```

**Environment**: Python 3.11, RDKit 2026.03.3, chemprop (pip install -e .)

**Data source**: https://github.com/kspieks/chemprop/blob/barrier_prediction/data.tar.gz (19 CSV files, 1.23M molecules)
