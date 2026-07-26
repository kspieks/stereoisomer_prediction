# Stereoisomer Prediction

Exploring how stereochemistry featurization affects molecular property prediction, with a focus on comparing chiral parity tags (`@`/`@@`) vs CIP absolute configuration (R/S) in chemprop and Random Forest models.

## Dataset

CMRT (Chiral Molecular Retention Time) dataset from [Baimacheva et al. 2025](https://link.springer.com/article/10.1186/s13321-025-01080-7): 1,929 enantiomeric pairs with chiral HPLC elution order labels. See `cmrt/README.md` for details.

## Environment Setup

Currently tested with Python v3.11 and Chemprop v2.2.2

### Option A: Install from requirements.txt

```bash
conda create --name chemprop python=3.11
conda activate chemprop
pip install -r requirements.txt
```

### Option B: Install chemprop from source (for development)

```bash
git clone https://github.com/chemprop/chemprop.git
cd chemprop
conda create --name chemprop python=3.11
conda activate chemprop
pip install -e .
```

Then install additional packages:
```bash
pip install scikit-learn pandas jupyter
```

> **Note:** Installing chemprop pulls in PyTorch, RDKit, Lightning, and NumPy automatically. The `pip install -e .` editable install is useful if you want to modify chemprop's featurization code directly.

## Project Structure
- `cmrt/` contains CMRT dataset and experiments comparing random forest with different chemprop featurizations.
- `reference/`contains background analysis on RDKit stereochemistry


## Key Findings

- Chemprop's default atom featurization encodes the SMILES parity tag (`@`/`@@` → CW/CCW), which does **not** reliably map to absolute stereochemistry (R/S).
- Across 1.23M molecules, CW maps to R ~50% and S ~50% of the time.
- Encoding CIP R/S labels directly gives the model access to chemically meaningful stereochemistry information. See `reference/stereochemistry_quick_reference.md` additional analysis.
