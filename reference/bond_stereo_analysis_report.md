# Bond Stereochemistry Encoding in RDKit and Chemprop: CIS/TRANS vs E/Z

## Motivation

Bond stereochemistry (E/Z, cis/trans) describes the spatial arrangement of substituents around a double bond and is distinct from atom chirality (R/S). In molecular property prediction, E/Z isomers can have vastly different biological activity (e.g., cis-platin vs trans-platin, tamoxifen E vs Z). Chemprop encodes bond stereo via `int(b.GetStereo())` as a one-hot bond-level feature. This report documents all possible values, their relationships to CIP E/Z labels, and empirical prevalence across 19 benchmark datasets (1.23M molecules) taken from chemprop v1.

## Methods

**Tools**: RDKit 2026.03.3, `rdCIPLabeler.AssignCIPLabels(mol)` for bond-level CIP (E/Z).

**Process**:
1. Read SMILES from each CSV dataset https://github.com/kspieks/chemprop/blob/barrier_prediction/data.tar.gz
2. Filter: skip molecules without `@`, `/`, or `\` in raw SMILES (catches both atom and bond stereo)
3. Parse with `Chem.MolFromSmiles()`, canonicalize, re-parse from canonical SMILES
4. Call `rdCIPLabeler.AssignCIPLabels(mol)` to assign bond CIP codes
5. For each bond: record `b.GetStereo()`, `int(b.GetStereo())`, `b.GetBondDir()`, and `bond.GetProp("_CIPCode")`

**Filter validation**: The combined filter (`@` or `/` or `\`) skipped 92.6% of molecules (1,141,618 / 1,232,875). Molecules with only `/`/`\` (bond stereo but no atom stereo) are correctly included — e.g., muv has 8,109 such molecules. Total runtime: 44 seconds.

## RDKit BondStereo Reference

### All 8 BondStereo Enum Values

| int | Enum Name | Description |
|---|---|---|
| 0 | `STEREONONE` | No bond stereochemistry |
| 1 | `STEREOANY` | Unspecified/either (wiggly bond in MOL files) |
| 2 | `STEREOZ` | Z isomer (older RDKit assignment path, from MOL/SDF files) |
| 3 | `STEREOE` | E isomer (older RDKit assignment path, from MOL/SDF files) |
| 4 | `STEREOCIS` | Cis isomer (canonical RDKit assignment from SMILES) |
| 5 | `STEREOTRANS` | Trans isomer (canonical RDKit assignment from SMILES) |
| 6 | `STEREOATROPCW` | Atropisomerism clockwise (restricted single bond rotation) |
| 7 | `STEREOATROPCCW` | Atropisomerism counter-clockwise |

### Python Dictionary Mapping

```python
BOND_STEREO_INT_TO_NAME = {
    0: "STEREONONE",
    1: "STEREOANY",
    2: "STEREOZ",
    3: "STEREOE",
    4: "STEREOCIS",
    5: "STEREOTRANS",
    6: "STEREOATROPCW",
    7: "STEREOATROPCCW",
}
```

### Why Both Z/E (2/3) and CIS/TRANS (4/5)?

RDKit has two code paths for assigning bond stereochemistry:
- **From SMILES** (`/` and `\`): assigns `STEREOCIS` (4) or `STEREOTRANS` (5)
- **From MOL/SDF files**: assigns `STEREOZ` (2) or `STEREOE` (3)

These are semantically equivalent — STEREOCIS ≈ STEREOZ, STEREOTRANS ≈ STEREOE — but have different integer values. Since chemprop reads SMILES as input, you will almost exclusively see values 4 and 5.

### Observed Values in Data

Across **1.23M molecules** from 19 datasets:

| Value | Name | Count |
|---|---|---|
| 5 | STEREOTRANS | 43,097 |
| 4 | STEREOCIS | 18,908 |
| 2 | STEREOZ | 3 |
| 1 | STEREOANY | 0 |
| 3 | STEREOE | 0 |
| 6 | STEREOATROPCW | 0 |
| 7 | STEREOATROPCCW | 0 |

Data-driven mapping confirmed:
```
STEREOCIS   -> int values observed: [4]
STEREOTRANS -> int values observed: [5]
STEREOZ     -> int values observed: [2]
```

## RDKit BondDir Reference

### All 7 BondDir Enum Values

| int | Enum Name | Description |
|---|---|---|
| 0 | `NONE` | No directional marker |
| 1 | `BEGINWEDGE` | Wedge bond (atom chirality, not stereo bonds) |
| 2 | `BEGINDASH` | Dashed wedge bond (atom chirality) |
| 3 | `ENDDOWNRIGHT` | `\` in SMILES (pointing down-right) |
| 4 | `ENDUPRIGHT` | `/` in SMILES (pointing up-right) |
| 5 | `EITHERDOUBLE` | Crossed double bond (unspecified stereo) |
| 6 | `UNKNOWN` | Unknown direction |

### Relationship Between BondDir and BondStereo

These are **different levels of abstraction**:
- **BondDir** is the low-level SMILES encoding: `/` and `\` markers placed on the **single bonds adjacent** to a double bond
- **BondStereo** is the interpreted result: CIS or TRANS assigned to the **double bond itself**

Each stereo double bond typically requires **two** BondDir markers (one on each side). In conjugated systems (`C=C-C=C`), a single bond's BondDir can serve both adjacent double bonds, so the ratio drops below 2:1.

**Empirical ratio**: BondDir total (121,468) / BondStereo total (62,008) = **1.96×** (close to 2.0, reduced by conjugated systems).

### Observed BondDir Values

| Value | Name | Count |
|---|---|---|
| 4 | ENDUPRIGHT | 98,298 |
| 3 | ENDDOWNRIGHT | 23,170 |

Values 1, 2, 5, 6 were not observed (BEGINWEDGE/BEGINDASH are for atom chirality wedge drawings, not bond stereo in SMILES context).

## Bond CIP Code Reference

### How to Assign

After calling `rdCIPLabeler.AssignCIPLabels(mol)`, stereo bonds get a `_CIPCode` property:

```python
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

mol = Chem.MolFromSmiles(smiles)
rdCIPLabeler.AssignCIPLabels(mol)

for bond in mol.GetBonds():
    if bond.HasProp("_CIPCode"):
        cip = bond.GetProp("_CIPCode")  # "E", "Z", or "z"
```

### Possible Values

| CIPCode | Count | Meaning |
|---|---|---|
| `"E"` | 43,107 | Entgegen — highest priority groups on opposite sides |
| `"Z"` | 18,897 | Zusammen — highest priority groups on same side |
| `"z"` | 1 | Pseudo-asymmetric bond (lowercase, analogous to atom r/s) |
| *(absent)* | — | Not a stereo bond |

Like atom CIP codes, bond CIP codes are stored as **string properties**, not an integer enum.

## Key Finding: CIS/TRANS Maps to Z/E ~99.93% of the Time

Unlike atom chirality where CW/CCW mapped to R/S only ~50% of the time, bond stereo is much better behaved:

| BondStereo | Maps to E | Maps to Z | Maps to z | Agreement |
|---|---|---|---|---|
| STEREOCIS | 29 | 18,345 | 1 | 99.8% → Z |
| STEREOTRANS | 42,064 | 14 | 0 | 99.97% → E |

**Why the discrepancy exists (rare cases)**:

CIS/TRANS is defined relative to `GetStereoAtoms()` — the specific pair of reference atoms RDKit uses. E/Z is defined by CIP priority rules. When the reference atoms are **not** the highest-priority substituents on each side of the double bond, CIS/TRANS can disagree with E/Z.

Example mismatch: `CSc1ccc(/C=C2/C(C)=C(CC(=O)O)c3cc(F)ccc32)cc1` — STEREOTRANS but CIP says Z, because the StereoAtoms chosen by RDKit are not the highest-priority groups.

**This is the same conceptual issue as atom chirality (CW/CCW vs R/S)**, but empirically it's negligible for bonds (0.07% mismatch) vs catastrophic for atoms (~50% mismatch). The reason: for canonical SMILES, RDKit's reference atom selection happens to align with CIP priorities in the vast majority of cases.

## Dataset Results

| Dataset | Molecules | Atom Stereo | Bond Stereo | % Any Stereo |
|---|---|---|---|---|
| clintox | 1,478 | 616 | 150 | 47.0% |
| bbbp | 2,039 | 637 | 148 | 35.2% |
| lipo | 4,200 | 1,127 | 80 | 28.2% |
| bace | 1,513 | 383 | 55 | 26.4% |
| tox21 | 7,831 | 1,321 | 465 | 21.5% |
| toxcast | 8,576 | 1,461 | 550 | 21.6% |
| sider | 1,427 | 288 | 30 | 21.4% |
| pcba | 437,929 | 30,817 | 47,534 | 17.4% |
| freesolv | 642 | 49 | 11 | 9.3% |
| muv | 93,087 | 0 | 8,109 | 8.7% |
| delaney | 1,128 | 0 | 6 | 0.5% |

**Datasets with NO stereochemistry at all**: chembl (456K), hiv (41K), pdbbind_core/full/refined, qm7, qm8, qm9.

### Notable Findings

- **Bond stereo is more prevalent than atom stereo**: 57,138 molecules have bond stereo vs 36,699 with atom chirality
- **muv** has 8,109 molecules with bond stereo but zero atom chirality — previously reported as "no stereocenters" in the atom-only analysis
- **delaney** has 6 molecules with bond stereo — also previously missed
- Breakdown: 34,108 molecules with atom stereo only, 54,547 with bond stereo only, 2,591 with both

## Chemprop's Current Bond Stereo Encoding

From `chemprop/featurizers/bond.py`:

```python
self.stereo = stereos or range(6)  # [0, 1, 2, 3, 4, 5]
# ...
stereo_bit, _ = self.one_hot_index(int(b.GetStereo()), self.stereo)
x[i + stereo_bit] = 1
```

Chemprop encodes `int(b.GetStereo())` as a one-hot vector over 6 values + 1 unknown pad, for a total of 7 bits:

| Slot | Value | Activated in practice? |
|---|---|---|
| 0 | STEREONONE | Yes (vast majority of bonds) |
| 1 | STEREOANY | No (never observed) |
| 2 | STEREOZ | Rarely (3 bonds across 1.23M molecules) |
| 3 | STEREOE | No (never observed) |
| 4 | STEREOCIS | Yes (18,908 bonds) |
| 5 | STEREOTRANS | Yes (43,097 bonds) |
| 6 | Unknown pad | For values ≥6 (atropisomerism, if encountered) |

### Is `range(6)` Sensible?

**Yes, for typical organic/drug datasets.** Here's why:

1. **Values 6–7 (atropisomerism)** describe restricted rotation around single bonds in biaryl systems. Not observed in any of our 1.23M molecules. Extremely rare outside specialized medicinal chemistry. If encountered, they fall into the "unknown" pad slot — which is acceptable graceful degradation.

2. **The Z/E (2/3) vs CIS/TRANS (4/5) redundancy** is harmless. Since chemprop reads SMILES input, only values 4/5 activate. Slots 2/3 are dead features the model learns to ignore. This wastes 2 one-hot dimensions but has no accuracy impact.

3. **STEREOANY (1)** represents intentionally unspecified stereo (wiggly bonds in MOL files). Never seen in curated SMILES datasets.

**If you wanted to optimize**, you could collapse to 3 meaningful states for SMILES input: `{none, cis, trans}`. But the current approach is compatible with multiple input formats (MOL, SDF, SMILES) and handles edge cases gracefully.

**When to expand to `range(8)`**: Only if working with datasets containing atropisomeric compounds (axially chiral biaryls, increasingly relevant in modern drug design). For standard MoleculeNet benchmarks, this is unnecessary.

## Comparison: Atom Stereo vs Bond Stereo Encoding Issues

| | Atom (CW/CCW vs R/S) | Bond (CIS/TRANS vs E/Z) |
|---|---|---|
| Conceptual issue | Same: encoding is relative to reference atoms, not CIP priorities | Same |
| Empirical mismatch | ~50% (CW doesn't reliably predict R or S) | 0.07% (CIS reliably predicts Z) |
| Why different? | Canonical atom ordering rarely aligns with CIP priority ordering | Canonical SMILES picks reference atoms that usually ARE highest priority |
| Need CIP labels? | **Yes** — essential for meaningful encoding | **No** — current encoding is 99.93% equivalent to E/Z |
| Recommendation | Replace with CIP R/S encoding | Current CIS/TRANS encoding is adequate |

## Reproducibility

| File | Purpose |
|---|---|
| `analyze_bond_stereo_datasets.py` | Full bond + atom stereo analysis across all 19 datasets |
| `bond_stereo_analysis.log` | Complete printed output from the full run |
| `bond_stereo_analysis_results.json` | Structured results: counts, examples, smiles_to_datasets dict |
| `bond_stereo_report_workspec.md` | Work spec for this report |

**Environment**: conda env `chemprop`, Python 3.11, RDKit 2026.03.3
