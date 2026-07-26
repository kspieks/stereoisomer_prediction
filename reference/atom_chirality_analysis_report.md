# Chirality Encoding in RDKit and Chemprop: @/@@ vs R/S

## Motivation

Chemprop encodes atom chirality as a one-hot feature over `int(a.GetChiralTag())`, which represents the local parity of a stereocenter (`@`/`@@` in SMILES). This captures clockwise (CW) vs counter-clockwise (CCW) spatial arrangement relative to the SMILES atom ordering. However, the chemically meaningful descriptor of absolute stereochemistry is the CIP assignment (R/S), which is determined by substituent priority rules. This report investigates whether @/@@ can serve as a proxy for R/S, documents all possible values for both representations, and summarizes empirical findings across 19 benchmark datasets (1.23M molecules) taken from chemprop v1.

## Methods

**Tools**: RDKit 2026.03.3 (via `pip install -e .` of chemprop), `rdCIPLabeler.AssignCIPLabels(mol)` for CIP assignment.

**Process**:
1. Read SMILES from each CSV dataset https://github.com/kspieks/chemprop/blob/barrier_prediction/data.tar.gz
2. Filter: skip molecules without `@` in the raw SMILES string
3. Parse with `Chem.MolFromSmiles()`, canonicalize, re-parse from canonical SMILES
4. Call `rdCIPLabeler.AssignCIPLabels(mol)` to assign CIP codes
5. For each atom: record `a.GetChiralTag()`, `int(a.GetChiralTag())`, and `atom.GetProp("_CIPCode")`

**Optimization**: Filtering by `@` presence skipped 97% of molecules (1,196,171 / 1,232,875) without RDKit parsing. This is a safe optimization because `@` or `@@` is the only way to encode a stereocenter in SMILES. Total runtime: 25 seconds.

## Datasets

19 CSV files from chemprop v1 datasets (standard MoleculeNet + ChEMBL benchmarks):

| Dataset | Molecules | % with Chirality |
|---|---|---|
| clintox | 1,478 | 41.7% |
| bbbp | 2,039 | 31.2% |
| lipo | 4,200 | 26.8% |
| bace | 1,513 | 25.3% |
| sider | 1,427 | 20.2% |
| toxcast | 8,576 | 17.0% |
| tox21 | 7,831 | 16.9% |
| freesolv | 642 | 7.6% |
| pcba | 437,929 | 7.0% |

**Datasets with NO stereocenters**: chembl (456K), delaney (1.1K), hiv (41K), muv (93K), pdbbind_core (168), pdbbind_full (9.9K), pdbbind_refined (3K), qm7 (6.8K), qm8 (21.8K), qm9 (134K).

Why these datasets lack stereocenters:
- **qm7/qm8/qm9**: Small molecules for quantum chemistry — mostly lack tetrahedral stereocenters
- **delaney/freesolv**: Small simple organic molecules (solubility datasets)
- **chembl/hiv/muv/pdbbind**: Stereochemistry was stripped during dataset creation (the raw molecules certainly have stereocenters, but the SMILES in these files don't encode them)

## RDKit ChiralTag Reference

### All 9 ChiralType Enum Values

| int | Enum Name | Description |
|---|---|---|
| 0 | `CHI_UNSPECIFIED` | No chirality specified (achiral atom or unassigned) |
| 1 | `CHI_TETRAHEDRAL_CW` | Clockwise tetrahedral — corresponds to `@@` in canonical SMILES |
| 2 | `CHI_TETRAHEDRAL_CCW` | Counter-clockwise tetrahedral — corresponds to `@` in canonical SMILES |
| 3 | `CHI_OTHER` | Other/unrecognized chirality |
| 4 | `CHI_TETRAHEDRAL` | Generic tetrahedral (direction unspecified) |
| 5 | `CHI_ALLENE` | Axial chirality (allenes, C=C=C) |
| 6 | `CHI_SQUAREPLANAR` | Square planar geometry |
| 7 | `CHI_TRIGONALBIPYRAMIDAL` | Trigonal bipyramidal geometry |
| 8 | `CHI_OCTAHEDRAL` | Octahedral geometry |

### Python Dictionary Mapping

```python
CHIRAL_TAG_INT_TO_NAME = {
    0: "CHI_UNSPECIFIED",
    1: "CHI_TETRAHEDRAL_CW",
    2: "CHI_TETRAHEDRAL_CCW",
    3: "CHI_OTHER",
    4: "CHI_TETRAHEDRAL",
    5: "CHI_ALLENE",
    6: "CHI_SQUAREPLANAR",
    7: "CHI_TRIGONALBIPYRAMIDAL",
    8: "CHI_OCTAHEDRAL",
}
```

### Observed Values in Data

Across **1.23M molecules** from 19 datasets, only values **1 (CW)** and **2 (CCW)** were observed for chiral atoms. Values 3–8 were **never** encountered. This was confirmed by the data-driven mapping:

```
CHI_TETRAHEDRAL_CW   -> int values observed: [1]  (57,160 atoms)
CHI_TETRAHEDRAL_CCW  -> int values observed: [2]  (59,805 atoms)
```

### Is It Sensible for Chemprop to Use Only 4 Tags?

Chemprop's default `chiral_tags=list(range(4))` encodes `[0, 1, 2, 3]` — unspecified, CW, CCW, and other. This is a reasonable choice:

- **Values 4–8 are exotic geometries** (allene, square planar, trigonal bipyramidal, octahedral) that occur almost exclusively in organometallic or coordination chemistry — not in typical drug-like or organic molecules found in these datasets.
- **Empirically confirmed**: zero occurrences of values 3–8 across 1.23M molecules.
- **For drug discovery / organic chemistry tasks**, 4 slots is sufficient. For inorganic or materials datasets, you'd want to expand to cover values 4–8.

### Example SMILES

| ChiralTag | Example (canonical SMILES) |
|---|---|
| CHI_TETRAHEDRAL_CW (1) | `CC[C@]1(O)CC[C@H]2[C@@H]3CCC4=CCCC[C@@H]4[C@H]3CC[C@@]21C` |
| CHI_TETRAHEDRAL_CCW (2) | `O=C(O)[C@H](O)c1ccccc1` |
| CHI_OTHER (3) | *(not observed in data)* |
| CHI_ALLENE (5) | *(not observed — would require allene stereochemistry)* |
| CHI_SQUAREPLANAR (6) | *(not observed — would require square planar metal complexes)* |

## CIP Code Reference

### What is CIP?

The **Cahn-Ingold-Prelog (CIP)** priority rules assign absolute stereochemistry labels to stereocenters based on substituent atomic numbers. Unlike @/@@, CIP labels are:
- **Molecule-intrinsic**: depend only on the molecular structure, not on how SMILES was written
- **Absolute**: R is always R regardless of atom ordering or canonicalization

### Possible Values

| CIPCode | Meaning | Count in Data |
|---|---|---|
| `"R"` | Rectus (right-handed) | 51,718 |
| `"S"` | Sinister (left-handed) | 64,987 |
| `"r"` | Pseudo-asymmetric center, r | 194 |
| `"s"` | Pseudo-asymmetric center, s | 66 |
| *(absent)* | Not a stereocenter | — |

### Python Dictionary Mapping

RDKit stores CIP codes as **string properties** (`atom.GetProp("_CIPCode")`), not as an integer enum. There is no built-in integer mapping. A proposed integer encoding for featurization:

```python
CIP_CODE_TO_INT = {
    None: 0,   # not a stereocenter (no _CIPCode property)
    "R": 1,
    "S": 2,
    "r": 3,    # pseudo-asymmetric
    "s": 4,    # pseudo-asymmetric
}
```

### What Are Pseudo-Asymmetric Centers (r/s)?

There are **four** possible `_CIPCode` values, not just two. Lowercase `r` and `s` denote [pseudo-asymmetric centers](https://en.wikipedia.org/wiki/Pseudoasymmetric_center) — stereocenters where two substituents differ only in their own chirality (e.g., an R-configured branch and an S-configured branch). The labels `r` and `s` (lowercase) are used instead of R/S. They are rare (260 atoms vs 116,705 R/S atoms in our data, ~0.2%) but chemically meaningful — they distinguish diastereomers.

### How to Assign CIP Labels

```python
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

mol = Chem.MolFromSmiles(smiles)
rdCIPLabeler.AssignCIPLabels(mol)

for atom in mol.GetAtoms():
    if atom.HasProp("_CIPCode"):
        cip = atom.GetProp("_CIPCode")  # "R", "S", "r", or "s"
```

## Key Finding: @/@@ Does NOT Map to R/S

### The Problem

CW/CCW (`@`/`@@`) encodes the **spatial arrangement of neighbors relative to the order they appear in the SMILES string**. R/S is determined by **CIP priority rules** which rank substituents by atomic number, mass, etc. These are fundamentally different: one depends on atom ordering, the other depends on substituent chemistry.

### Canonicalizing SMILES Does NOT Fix the Problem

Even with **canonical SMILES**, CW/CCW does not cleanly map to R/S:

| Canonical SMILES | ChiralTag | CIPCode |
|---|---|---|
| `C[C@H](N)C(=O)O` (L-alanine) | CCW | S |
| `N[C@@H](CS)C(=O)O` (L-cysteine) | CW | **R** |
| `N[C@@H](CC(=O)O)C(=O)O` (L-aspartate) | CW | **S** |
| `O[C@@H](Cl)Br` | CW | **S** |
| `F[C@H](Cl)Br` | CCW | **R** |

**CW maps to both R and S** depending on the molecule. Same for CCW.

Canonicalization fixes the *atom ordering* problem (so the same molecule always gets the same CW/CCW assignment), but it doesn't fix the fundamental disconnect: whether CW maps to R or S depends on whether the canonical atom ordering happens to agree with CIP priority ordering for that specific molecule. That relationship changes from molecule to molecule based on what the substituents are.

**Concrete example**: L-cysteine `N[C@@H](CS)C(=O)O` → canonical CW → **R** (because sulfur changes the CIP priorities). L-aspartate `N[C@@H](CC(=O)O)C(=O)O` → canonical CW → **S**. Same `@@` tag, same canonical ordering pattern, different R/S. The CIP priority ranking of substituents is what determines R vs S, and that's molecular-structure-dependent.

### Why This Happens

Whether the canonical atom ordering happens to agree with CIP priority ordering is molecule-dependent. For L-cysteine, the presence of sulfur (higher priority than oxygen) flips the CIP ranking relative to the canonical atom order, so the same `@@` (CW) pattern maps to R instead of S.

### Bottom Line

There is no fixed mapping from `@`/`@@` to R/S. If you want to encode absolute stereochemistry (R/S), you must explicitly compute CIP labels via `rdCIPLabeler.AssignCIPLabels(mol)` and use `atom.GetProp("_CIPCode")`. There is no shortcut through `@`/`@@`.

### Sanity Check

Total chiral atoms from ChiralTag (CW + CCW) = 57,160 + 59,805 = **116,965**
Total CIP-labeled atoms (R + S + r + s) = 51,718 + 64,987 + 194 + 66 = **116,965**

Every atom with a ChiralTag gets a CIP label, and vice versa. The two representations are alternative views of the same stereocenters.

## Investigating Molecules with @ But No Chirality

Two molecules in freesolv had `@` in their raw SMILES but **no chirality after parsing**:

```
Raw:       CCC[N@@](CC1CC1)c2c(cc(cc2[N+](=O)[O-])C(F)(F)F)[N+](=O)[O-]
Canonical: CCCN(CC1CC1)c1c([N+](=O)[O-])cc(C(F)(F)F)cc1[N+](=O)[O-]

Raw:       CCCC[N@](CC)c1c(cc(cc1[N+](=O)[O-])C(F)(F)F)[N+](=O)[O-]
Canonical: CCCCN(CC)c1c([N+](=O)[O-])cc(C(F)(F)F)cc1[N+](=O)[O-]
```

**Explanation**: Both have `@` on a trivalent nitrogen (`N@@`, `N@`). RDKit correctly recognizes that **sp2 nitrogens with lone pairs are not stereocenters** — the nitrogen is planar (conjugated into the aromatic ring), so the `@` annotation is chemically invalid. RDKit silently removes it during parsing. This is expected and correct behavior. The `@` filter is still safe; these are rare edge cases (2/1.23M) where the input SMILES has incorrect stereochemistry annotations.

## E/Z (Bond) Stereochemistry: Should It Be Encoded?

This analysis focused on **atom** chirality (tetrahedral stereocenters). A separate type of stereochemistry exists for **double bonds**: E/Z (cis/trans) isomerism, encoded in SMILES as `/` and `\`.

**Should E/Z be encoded as a feature?**

Yes, there is evidence it would be useful:

1. **Chemical relevance**: E/Z isomers can have vastly different biological activity (e.g., cis-platin vs trans-platin, tamoxifen E vs Z). For drug discovery datasets, ignoring bond stereochemistry loses meaningful information.
2. **Literature support**: The extended atomic featurization study (Wojtuch et al., 2023) and the D-MPNN paper (Yang et al., 2019) both include bond stereo features. The D-MPNN bond representation includes 7 stereo features and the paper shows that "rich bond representation can substitute for additional atomic features."
3. **Chemprop already encodes it**: In chemprop, bond stereochemistry is handled via `bond.GetStereo()` as a bond-level feature (separate from atom features). The possible values are `STEREONONE`, `STEREOANY`, `STEREOZ`, `STEREOE`, `STEREOCIS`, `STEREOTRANS`.

The same concern applies: chemprop encodes the RDKit stereo enum directly rather than the absolute E/Z assignment. However, unlike R/S, the E/Z assignments from `GetStereo()` are already absolute (E and Z are defined by CIP priorities on the bond substituents), so the current encoding is more reliable than the @/@@ → R/S situation.

## Related RDKit Properties: `_ChiralityPossible` and `_CIPRank`

### `_ChiralityPossible`

An atom with `ChiralTag = CHI_UNSPECIFIED` (int=0) could be one of three things:
- **A) Truly achiral** — symmetric substituents (e.g., central C in `CC(C)C`)
- **B) Potential stereocenter with unspecified chirality** — 4 different groups but no `@` in SMILES (e.g., `CC(O)(F)Cl`)
- **C) Non-tetrahedral atom** — sp2 carbon, aromatic atom, etc.

**Can `_ChiralityPossible` distinguish these?** No. This property is only set on atoms that *already have* a ChiralTag (i.e., `@` is present in SMILES) and is an internal validation flag confirming the tagged atom is a valid stereocenter. It is **not** set on unspecified-but-possible stereocenters.

**The right tool is `FindPotentialStereo(mol)`**:

```python
from rdkit.Chem import FindPotentialStereo

mol = Chem.MolFromSmiles("OC(F)C(O)F")  # Two potential centers, no @ marks
stereo_info = FindPotentialStereo(mol)
for si in stereo_info:
    print(f"Atom {si.centeredOn}: specified={si.specified}")
# Atom 1: specified=Unspecified
# Atom 3: specified=Unspecified
```

| Input | FindPotentialStereo result |
|---|---|
| `CC(C)C` (symmetric) | No potential stereocenters |
| `CC(O)(F)Cl` (unspecified) | Atom 1: `specified=Unspecified` |
| `C[C@H](O)F` (specified) | Atom 1: `specified=Specified` |
| `CC=CC` (sp2) | No potential stereocenters |

**Is this relevant for chemprop featurization?** It's overkill for most use cases:
- If datasets have **complete stereo annotations** (every real stereocenter has `@`), then `CHI_UNSPECIFIED` reliably means "achiral." No ambiguity.
- If datasets have **stripped stereo** (chembl, hiv, muv — zero `@`), the entire molecule is missing stereo. Flagging individual atoms doesn't help.
- It would only matter in a mixed scenario where a molecule has *some* `@` marks but is *missing* others — a data quality issue, not a featurization question.

### `_CIPRank`

Every atom gets a `_CIPRank` integer after stereo assignment. This is the **internal ranking** used by RDKit to compute R/S:

```python
mol = Chem.MolFromSmiles("C[C@H](O)CC")
Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
for atom in mol.GetAtoms():
    rank = atom.GetIntProp("_CIPRank")
    # Atom 0 (C): rank=1, Atom 1 (C*): rank=3, Atom 2 (O): rank=4, ...
```

**Not useful as a feature** for three reasons:
1. It's a **relative** ranking within a molecule — not an absolute atomic property
2. Its scale changes molecule to molecule (rank 3 in one molecule means something different from rank 3 in another)
3. The useful output of the ranking is already captured in `_CIPCode` (R/S)

`_CIPRank` is an intermediate computation while `_CIPCode` is the final answer.

## Implications for Chemprop Featurization

| | Current (ChiralTag) | Proposed (CIPCode) |
|---|---|---|
| **Encoding** | One-hot over [0, 1, 2, 3] | One-hot over [0, 1, 2, 3, 4] |
| **Values** | unspecified, CW, CCW, other | unspecified, R, S, r, s |
| **Meaning** | Local parity (depends on atom ordering) | Absolute stereochemistry |
| **Preprocessing** | None (from SMILES parse) | Requires `rdCIPLabeler.AssignCIPLabels(mol)` |
| **Trade-off** | Fast, but not chemically absolute | Slightly slower, chemically meaningful |

## Reproducibility

| File | Purpose |
|---|---|
| `explore_chirality.py` | Initial exploration: ChiralTag enum values, CIP codes, @/@@→R/S counter-examples |
| `analyze_chirality_datasets.py` | Full analysis across all 19 datasets with timing and @ filter optimization |
| `chirality_analysis.log` | Complete printed output from the full run |
| `chirality_analysis_results.json` | Structured results: counts, examples, smiles_to_datasets dict |
| `chirality_report_workspec.md` | Work spec for this report |

**Environment**: conda env `chemprop`, Python 3.11, RDKit 2026.03.3
