# Stereochemistry Featurization: Quick Reference & Action Items

## Preprocessing Functions

### `rdCIPLabeler.AssignCIPLabels(mol)`

The **newer, more correct** CIP implementation. Use this one.

```python
from rdkit.Chem import rdCIPLabeler

mol = Chem.MolFromSmiles(smiles)
rdCIPLabeler.AssignCIPLabels(mol)
# Now atom.GetProp("_CIPCode") returns R/S/r/s
# And bond.GetProp("_CIPCode") returns E/Z/e/z
```

- Ported from the Centres library, aims to cover all 6 CIP rules
- Sets `_CIPCode` on **both** atoms (R/S/r/s) and bonds (E/Z/e/z)
- Must be called explicitly — not automatic
- Supersedes the legacy function below

### `Chem.AssignStereochemistry(mol, cleanIt=True, force=True)`

**Legacy** CIP implementation. You generally don't need to call this.

- Only covers basic CIP rules (1–3 of 6)
- Sets `_CIPCode` on atoms only (NOT bonds)
- Called **automatically** by `Chem.MolFromSmiles()` — so `_CIPCode` is often already present after parsing
- `rdCIPLabeler.AssignCIPLabels()` overwrites its results with more correct assignments

**Action item**: Use `rdCIPLabeler.AssignCIPLabels(mol)` for preprocessing. Don't call `AssignStereochemistry` separately.

### `FindPotentialStereo(mol)`

Identifies atoms/bonds that **could** be stereocenters, regardless of whether `@` is in the SMILES.

```python
from rdkit.Chem import FindPotentialStereo

stereo_info = FindPotentialStereo(mol)
for si in stereo_info:
    # si.centeredOn = atom/bond index
    # si.specified = Specified | Unspecified
```

- Returns `specified=Specified` for centers with `@` in SMILES
- Returns `specified=Unspecified` for centers that *could* be chiral but have no `@`
- Returns nothing for truly achiral atoms (symmetric substituents)

**Action item**: Not needed for featurization. Only useful for data quality auditing (finding molecules with missing stereo annotations).

---

## Atom Stereochemistry

### `a.GetChiralTag()` / `int(a.GetChiralTag())`

Encodes local parity (`@`/`@@` in SMILES). Does **NOT** map to R/S.

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

**Data finding**: Only values 0, 1, 2 observed across 1.23M molecules. Values 3–8 never seen.

**Action item**: Chemprop's default `chiral_tags=list(range(4))` ([source](https://github.com/chemprop/chemprop/blob/main/chemprop/featurizers/atom.py#L167)) is sufficient for drug-like/organic molecules. Values 4–8 are only relevant for organometallic or coordination chemistry datasets.

**Critical caveat**: CW/CCW does NOT reliably map to R/S (~50% mismatch even with canonical SMILES). Canonicalization fixes atom ordering but not the fundamental disconnect — whether CW maps to R or S depends on whether the canonical atom ordering agrees with CIP priority ordering, which is molecule-specific.

---

### `atom.GetProp("_CIPCode")`

Encodes absolute stereochemistry (R/S). Requires `rdCIPLabeler.AssignCIPLabels(mol)` before featurization.

| CIPCode | Meaning |
|---|---|
| `"R"` | Rectus (right-handed) stereocenter |
| `"S"` | Sinister (left-handed) stereocenter |
| `"r"` | Pseudo-asymmetric center, r — [wiki](https://en.wikipedia.org/wiki/Pseudoasymmetric_center) |
| `"s"` | Pseudo-asymmetric center, s |
| *(absent)* | Not a stereocenter |

```python
# No integer enum exists in RDKit — it's a string property.
# Proposed encoding for featurization:
CIP_CODE_TO_INT = {
    None: 0,   # not a stereocenter
    "R": 1,
    "S": 2,
    "r": 3,
    "s": 4,
}
```

**Action item**: ✅ **Implemented.** Added `cip_codes` parameter to `MultiHotAtomFeaturizer`. Usage:

```python
from rdkit.Chem import rdCIPLabeler

mol = Chem.MolFromSmiles(smiles)
rdCIPLabeler.AssignCIPLabels(mol)  # MUST be called before featurization

featurizer = MultiHotAtomFeaturizer.v2(cip_codes=[None, "R", "S", "r", "s"])
```

This adds 6 bits (5 values + 1 unknown pad). Backward compatible — defaults to `None` (disabled).

---

### `_ChiralityPossible` and `_CIPRank`

**Not necessary for featurization.**

- **`_ChiralityPossible`**: Only set on atoms that already have a ChiralTag (`@` in SMILES). Cannot discover unspecified-but-possible stereocenters. Use `FindPotentialStereo(mol)` if you need that (but it's overkill for featurization).

- **`_CIPRank`**: Internal per-molecule relative ranking used to compute R/S. Not an absolute atomic property — rank 3 in one molecule means something different from rank 3 in another. The useful output is already in `_CIPCode`.

**Action item**: No changes needed. Ignore both for featurization.

---

## Bond Stereochemistry

### `bond.GetStereo()` / `int(bond.GetStereo())`

Encodes cis/trans stereochemistry on double bonds.

| int | Enum Name | Description |
|---|---|---|
| 0 | `STEREONONE` | No bond stereochemistry |
| 1 | `STEREOANY` | Unspecified/either (wiggly bond) |
| 2 | `STEREOZ` | Z isomer (older RDKit code path, from MOL/SDF files) |
| 3 | `STEREOE` | E isomer (older RDKit code path, from MOL/SDF files) |
| 4 | `STEREOCIS` | Cis isomer (canonical RDKit assignment from SMILES) |
| 5 | `STEREOTRANS` | Trans isomer (canonical RDKit assignment from SMILES) |
| 6 | `STEREOATROPCW` | Atropisomerism clockwise |
| 7 | `STEREOATROPCCW` | Atropisomerism counter-clockwise |

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

**Data finding**: Only values 0, 2, 4, 5 observed (primarily 4 and 5). Values 2/3 are redundant with 4/5 (same info from different code paths). Values 6/7 never seen.

**Action item**: Chemprop's default `stereos=range(6)` ([source](https://github.com/chemprop/chemprop/blob/main/chemprop/featurizers/bond.py)) is sufficient. No change needed. Includes slots for values 2/3 (harmless dead bits when input is SMILES) and values 6/7 fall to the unknown pad.

**Key difference from atom chirality**: Unlike CW/CCW→R/S (50% mismatch), CIS/TRANS maps to E/Z **99.93%** of the time. The current encoding is already nearly equivalent to absolute E/Z.

---

### `bond.GetBondDir()` / `int(bond.GetBondDir())`

Low-level SMILES encoding: `/` and `\` markers on **single bonds adjacent** to a double bond.

| int | Enum Name | Description |
|---|---|---|
| 0 | `NONE` | No directional marker |
| 1 | `BEGINWEDGE` | Wedge bond (for atom chirality drawing, not bond stereo) |
| 2 | `BEGINDASH` | Dashed wedge bond (for atom chirality drawing) |
| 3 | `ENDDOWNRIGHT` | `\` in SMILES |
| 4 | `ENDUPRIGHT` | `/` in SMILES |
| 5 | `EITHERDOUBLE` | Crossed double bond (unspecified stereo) |
| 6 | `UNKNOWN` | Unknown direction |

**Action item**: **Do NOT add.** `GetBondDir()` is the raw SMILES-level encoding. `GetStereo()` is the interpreted result at the correct abstraction level. Chemprop already uses `GetStereo()` which is the right choice.

---

### `bond.GetProp("_CIPCode")` (E/Z)

Absolute E/Z labels on bonds after `rdCIPLabeler.AssignCIPLabels(mol)`. Possible values: `"E"`, `"Z"`, `"z"` (pseudo, extremely rare).

**Action item**: **Do NOT add.** The existing `GetStereo()` encoding (CIS/TRANS) already agrees with CIP E/Z 99.93% of the time. Adding bond CIP codes would require the same `AssignCIPLabels` preprocessing cost for only 0.07% improvement. Not worth it.

---

## Summary of Action Items

| Feature | Action | Rationale |
|---|---|---|
| `a.GetChiralTag()` range | No change | `range(4)` covers all observed values |
| `atom._CIPCode` (R/S/r/s) | ✅ **Added** | Essential — CW/CCW has ~50% mismatch with R/S |
| `_ChiralityPossible` | No change | Not useful for featurization |
| `_CIPRank` | No change | Internal relative ranking, not a feature |
| `b.GetStereo()` range | No change | `range(6)` covers all observed values |
| `bond.GetBondDir()` | Do NOT add | Lower abstraction than `GetStereo()` |
| `bond._CIPCode` (E/Z) | Do NOT add | `GetStereo()` already 99.93% equivalent |
