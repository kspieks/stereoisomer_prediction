"""
Analyze bond stereochemistry across chemprop v1 datasets.

For each molecule:
- Count BondStereo values (GetStereo()) for all bonds
- Count BondDir values (GetBondDir()) for all bonds
- After AssignCIPLabels, count bond-level _CIPCode (E/Z)
- Store up to 5 example molecules per category
- Track BondStereo -> int mapping for verification
- Track which dataset each canonical SMILES comes from

Optimization: skip molecules without @, /, or \ in raw SMILES
(these are the only characters that encode stereochemistry in SMILES).

Usage:
    conda activate chemprop
    python analyze_bond_stereo_datasets.py
    python analyze_bond_stereo_datasets.py --data-dir /path/to/data --output ../results/bond_stereo_analysis_results.json
"""

import argparse
import csv
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.DisableLog("rdApp.*")

# ============================================================
# Configuration
# ============================================================
# Download data from here: https://github.com/kspieks/chemprop/blob/barrier_prediction/data.tar.gz
# then unzip and update _DEFAULT_DATA_DIR

# Defaults are relative to this script's location (reference/scripts/)
_SCRIPT_DIR = Path(__file__).resolve().parent
_DEFAULT_DATA_DIR = _SCRIPT_DIR / "../../data"
_DEFAULT_OUTPUT = _SCRIPT_DIR / "../results/bond_stereo_analysis_results.json"

# (now controlled via --single-file CLI arg)
MAX_EXAMPLES = 5

# Reference mappings (from RDKit enums)
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

BOND_DIR_INT_TO_NAME = {
    0: "NONE",
    1: "BEGINWEDGE",
    2: "BEGINDASH",
    3: "ENDDOWNRIGHT",
    4: "ENDUPRIGHT",
    5: "EITHERDOUBLE",
    6: "UNKNOWN",
}


def get_smiles_column(header):
    """Determine which column contains SMILES."""
    header_lower = [h.lower().strip() for h in header]
    if "smiles" in header_lower:
        return header_lower.index("smiles")
    elif "mol" in header_lower:
        return header_lower.index("mol")
    else:
        return 0


def has_stereo_chars(smi):
    """Check if SMILES contains any stereochemistry indicators (@ for atom, / \\ for bond)."""
    return "@" in smi or "/" in smi or "\\" in smi


def add_example(examples_dict, key, canon_smi):
    """Add canon_smi to examples_dict[key] if under MAX_EXAMPLES and not duplicate."""
    if len(examples_dict[key]) < MAX_EXAMPLES:
        if canon_smi not in examples_dict[key]:
            examples_dict[key].append(canon_smi)


def process_datasets(data_dir: Path = _DEFAULT_DATA_DIR, output_path: Path = _DEFAULT_OUTPUT, single_file: str | None = None):
    # ============================================================
    # Data structures
    # ============================================================
    # Bond stereo (GetStereo) — only for bonds with stereo != STEREONONE
    bond_stereo_counts = Counter()          # str(stereo) -> count
    bond_stereo_int_counts = Counter()      # int(stereo) -> count

    # Bond direction (GetBondDir) — only for bonds with dir != NONE
    bond_dir_counts = Counter()             # str(dir) -> count
    bond_dir_int_counts = Counter()         # int(dir) -> count

    # Bond CIP codes (after AssignCIPLabels)
    bond_cip_counts = Counter()             # E/Z -> count

    # Atom chirality (included for completeness)
    chiral_tag_counts = Counter()
    chiral_tag_int_counts = Counter()
    atom_cip_counts = Counter()

    # Examples
    bond_stereo_examples = defaultdict(list)
    bond_dir_examples = defaultdict(list)
    bond_cip_examples = defaultdict(list)
    chiral_tag_examples = defaultdict(list)
    atom_cip_examples = defaultdict(list)

    # Mapping verification
    bond_stereo_to_int_mapping = defaultdict(set)
    bond_dir_to_int_mapping = defaultdict(set)

    # Dataset tracking (only for molecules with any stereo)
    smiles_to_datasets = defaultdict(list)

    # Per-dataset stats
    dataset_stats = {}

    # Global stats
    total_molecules = 0
    total_stereo_molecules = 0
    total_skipped = 0
    failed_parses = 0
    mols_with_bond_stereo = 0
    mols_with_atom_stereo = 0
    mols_with_both = 0

    # ============================================================
    # Determine files
    # ============================================================
    if single_file:
        csv_files = [data_dir / single_file]
    else:
        csv_files = sorted(data_dir.glob("*.csv"))

    print(f"Processing {len(csv_files)} file(s)...")
    print()

    overall_start = time.time()

    # ============================================================
    # Main loop
    # ============================================================
    for csv_path in csv_files:
        dataset_name = csv_path.stem
        file_start = time.time()
        file_mol_count = 0
        file_stereo_count = 0
        file_fail_count = 0
        file_skipped = 0
        file_bond_stereo_count = 0
        file_atom_stereo_count = 0

        with open(csv_path, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            smi_col = get_smiles_column(header)

            for row in reader:
                if not row or smi_col >= len(row):
                    continue

                raw_smi = row[smi_col].strip()
                if not raw_smi:
                    continue

                file_mol_count += 1
                total_molecules += 1

                # Filter: must have @, /, or \ to have any stereochemistry
                if not has_stereo_chars(raw_smi):
                    file_skipped += 1
                    total_skipped += 1
                    continue

                # Parse and canonicalize
                mol = Chem.MolFromSmiles(raw_smi)
                if mol is None:
                    file_fail_count += 1
                    failed_parses += 1
                    continue

                canon_smi = Chem.MolToSmiles(mol)

                # Re-parse from canonical SMILES
                mol2 = Chem.MolFromSmiles(canon_smi)
                rdCIPLabeler.AssignCIPLabels(mol2)

                has_atom_chiral = False
                has_bond_stereo = False

                # --- Atom chirality ---
                for atom in mol2.GetAtoms():
                    tag = atom.GetChiralTag()
                    tag_str = str(tag)
                    tag_int = int(tag)

                    if tag_int != 0:
                        has_atom_chiral = True
                        chiral_tag_counts[tag_str] += 1
                        chiral_tag_int_counts[tag_int] += 1
                        add_example(chiral_tag_examples, tag_str, canon_smi)

                    if atom.HasProp("_CIPCode"):
                        cip = atom.GetProp("_CIPCode")   # R, S, r, or s
                        atom_cip_counts[cip] += 1
                        add_example(atom_cip_examples, cip, canon_smi)

                # --- Bond stereochemistry ---
                for bond in mol2.GetBonds():
                    # GetStereo (E/Z/cis/trans on double bonds)
                    stereo = bond.GetStereo()
                    stereo_str = str(stereo)
                    stereo_int = int(stereo)

                    if stereo_int != 0:  # Not STEREONONE
                        has_bond_stereo = True
                        bond_stereo_counts[stereo_str] += 1
                        bond_stereo_int_counts[stereo_int] += 1
                        bond_stereo_to_int_mapping[stereo_str].add(stereo_int)
                        add_example(bond_stereo_examples, stereo_str, canon_smi)

                    # GetBondDir (/ and \ directions on single bonds adjacent to double bond)
                    bdir = bond.GetBondDir()
                    bdir_str = str(bdir)
                    bdir_int = int(bdir)

                    if bdir_int != 0:  # Not NONE
                        bond_dir_counts[bdir_str] += 1
                        bond_dir_int_counts[bdir_int] += 1
                        bond_dir_to_int_mapping[bdir_str].add(bdir_int)
                        add_example(bond_dir_examples, bdir_str, canon_smi)

                    # Bond-level CIP code (after AssignCIPLabels)
                    if bond.HasProp("_CIPCode"):
                        bcip = bond.GetProp("_CIPCode")   # E, Z, e, or z
                        bond_cip_counts[bcip] += 1
                        add_example(bond_cip_examples, bcip, canon_smi)

                # Track stats
                if has_atom_chiral or has_bond_stereo:
                    file_stereo_count += 1
                    total_stereo_molecules += 1
                    smiles_to_datasets[canon_smi].append(dataset_name)

                if has_atom_chiral:
                    file_atom_stereo_count += 1
                    mols_with_atom_stereo += 1
                if has_bond_stereo:
                    file_bond_stereo_count += 1
                    mols_with_bond_stereo += 1
                if has_atom_chiral and has_bond_stereo:
                    mols_with_both += 1

        file_elapsed = time.time() - file_start
        avg_ms = (file_elapsed / file_mol_count * 1000) if file_mol_count > 0 else 0
        processed_count = file_mol_count - file_skipped

        dataset_stats[dataset_name] = {
            "total_molecules": file_mol_count,
            "skipped_no_stereo_chars": file_skipped,
            "processed_with_rdkit": processed_count,
            "failed_parses": file_fail_count,
            "molecules_with_any_stereo": file_stereo_count,
            "molecules_with_atom_stereo": file_atom_stereo_count,
            "molecules_with_bond_stereo": file_bond_stereo_count,
            "pct_any_stereo": (file_stereo_count / file_mol_count * 100) if file_mol_count > 0 else 0,
            "time_seconds": round(file_elapsed, 2),
            "avg_ms_per_molecule": round(avg_ms, 3),
        }

        print(
            f"  {csv_path.name:<25s} | "
            f"{file_mol_count:>7d} mols | "
            f"{file_skipped:>7d} skip | "
            f"{processed_count:>7d} parsed | "
            f"atom:{file_atom_stereo_count:>5d} bond:{file_bond_stereo_count:>5d} | "
            f"{file_elapsed:.2f}s"
        )

    overall_elapsed = time.time() - overall_start

    # ============================================================
    # Report (printed to stdout; tee to .log file for persistence)
    # ============================================================
    print(f"\n{'=' * 90}")
    print("RESULTS")
    print(f"{'=' * 90}")

    print(f"\nTotal molecules across all datasets: {total_molecules}")
    print(f"Skipped (no @, /, or \\ in SMILES):  {total_skipped}")
    print(f"Processed with RDKit:                {total_molecules - total_skipped}")
    print(f"Failed parses:                       {failed_parses}")
    print(f"Molecules with any stereochemistry:  {total_stereo_molecules}")
    print(f"  - Atom chirality only:             {mols_with_atom_stereo - mols_with_both}")
    print(f"  - Bond stereo only:                {mols_with_bond_stereo - mols_with_both}")
    print(f"  - Both:                            {mols_with_both}")
    print(f"Unique stereo canonical SMILES:      {len(smiles_to_datasets)}")
    print(f"\nTotal time: {overall_elapsed:.2f}s")

    # --- Bond Stereo (GetStereo) ---
    print(f"\n{'─' * 90}")
    print("bond.GetStereo() value counts (bonds with stereo != STEREONONE):")
    print(f"{'─' * 90}")
    for s, count in bond_stereo_counts.most_common():
        print(f"  {s:<25s}: {count:>10d}")

    print(f"\n{'─' * 90}")
    print("int(bond.GetStereo()) value counts:")
    print(f"{'─' * 90}")
    for i, count in bond_stereo_int_counts.most_common():
        print(f"  {i}: {count:>10d}")

    # --- Bond Dir ---
    print(f"\n{'─' * 90}")
    print("bond.GetBondDir() value counts (bonds with dir != NONE):")
    print(f"{'─' * 90}")
    for s, count in bond_dir_counts.most_common():
        print(f"  {s:<25s}: {count:>10d}")

    print(f"\n{'─' * 90}")
    print("int(bond.GetBondDir()) value counts:")
    print(f"{'─' * 90}")
    for i, count in bond_dir_int_counts.most_common():
        print(f"  {i}: {count:>10d}")

    # --- Bond CIP ---
    print(f"\n{'─' * 90}")
    print("bond.GetProp('_CIPCode') value counts (after AssignCIPLabels):")
    print(f"{'─' * 90}")
    for cip, count in bond_cip_counts.most_common():
        print(f"  {cip}: {count:>10d}")

    # --- Atom chirality (summary) ---
    print(f"\n{'─' * 90}")
    print("Atom chirality summary (same as atom analysis):")
    print(f"{'─' * 90}")
    print("  a.GetChiralTag():")
    for s, count in chiral_tag_counts.most_common():
        print(f"    {s:<30s}: {count:>10d}")
    print("  atom.GetProp('_CIPCode'):")
    for cip, count in atom_cip_counts.most_common():
        print(f"    {cip}: {count:>10d}")

    # --- Mapping verification ---
    print(f"\n{'─' * 90}")
    print("BondStereo str -> int mapping:")
    print(f"{'─' * 90}")
    for s, ints in sorted(bond_stereo_to_int_mapping.items(), key=lambda x: min(x[1])):
        print(f"  {s:<25s} -> {sorted(ints)}")

    print(f"\n{'─' * 90}")
    print("BondDir str -> int mapping:")
    print(f"{'─' * 90}")
    for s, ints in sorted(bond_dir_to_int_mapping.items(), key=lambda x: min(x[1])):
        print(f"  {s:<25s} -> {sorted(ints)}")

    # --- Examples ---
    print(f"\n{'─' * 90}")
    print(f"Example molecules per BondStereo (up to {MAX_EXAMPLES}):")
    print(f"{'─' * 90}")
    for s in sorted(bond_stereo_examples.keys()):
        print(f"\n  {s}:")
        for smi in bond_stereo_examples[s]:
            print(f"    {smi}")

    print(f"\n{'─' * 90}")
    print(f"Example molecules per BondDir (up to {MAX_EXAMPLES}):")
    print(f"{'─' * 90}")
    for s in sorted(bond_dir_examples.keys()):
        print(f"\n  {s}:")
        for smi in bond_dir_examples[s]:
            print(f"    {smi}")

    print(f"\n{'─' * 90}")
    print(f"Example molecules per bond CIPCode (up to {MAX_EXAMPLES}):")
    print(f"{'─' * 90}")
    for cip in sorted(bond_cip_examples.keys()):
        print(f"\n  {cip}:")
        for smi in bond_cip_examples[cip]:
            print(f"    {smi}")

    # --- Per-dataset summary ---
    print(f"\n{'─' * 90}")
    print("Per-dataset summary:")
    print(f"{'─' * 90}")
    print(f"  {'Dataset':<20s} | {'Total':>7s} | {'Skip':>7s} | {'Parsed':>7s} | {'AtomSt':>6s} | {'BondSt':>6s} | {'%Stereo':>7s} | {'Time':>5s}")
    print(f"  {'─' * 20}-+-{'─' * 7}-+-{'─' * 7}-+-{'─' * 7}-+-{'─' * 6}-+-{'─' * 6}-+-{'─' * 7}-+-{'─' * 5}")
    for name in sorted(dataset_stats.keys()):
        s = dataset_stats[name]
        print(
            f"  {name:<20s} | {s['total_molecules']:>7d} | {s['skipped_no_stereo_chars']:>7d} | "
            f"{s['processed_with_rdkit']:>7d} | {s['molecules_with_atom_stereo']:>6d} | "
            f"{s['molecules_with_bond_stereo']:>6d} | {s['pct_any_stereo']:>6.1f}% | {s['time_seconds']:>4.1f}s"
        )

    # Datasets with no stereo at all
    no_stereo_datasets = [name for name, s in dataset_stats.items() if s["molecules_with_any_stereo"] == 0]
    if no_stereo_datasets:
        print(f"\n  Datasets with NO stereochemistry at all:")
        for name in no_stereo_datasets:
            print(f"    - {name} ({dataset_stats[name]['total_molecules']} molecules)")

    # ============================================================
    # Save structured results to JSON for later inspection
    # ============================================================
    results = {
        "config": {
            "data_dir": str(data_dir),
            "files_processed": [p.name for p in csv_files],
            "single_file_mode": single_file,
            "filter": "@ or / or \\ in raw SMILES",
        },
        "timing": {
            "total_seconds": round(overall_elapsed, 2),
        },
        "stats": {
            "total_molecules": total_molecules,
            "skipped_no_stereo_chars": total_skipped,
            "processed_with_rdkit": total_molecules - total_skipped,
            "failed_parses": failed_parses,
            "molecules_with_any_stereo": total_stereo_molecules,
            "molecules_with_atom_stereo": mols_with_atom_stereo,
            "molecules_with_bond_stereo": mols_with_bond_stereo,
            "molecules_with_both": mols_with_both,
            "unique_stereo_canonical_smiles": len(smiles_to_datasets),
        },
        "datasets_with_no_stereo": no_stereo_datasets,
        "dataset_stats": dataset_stats,
        "bond_stereo_int_to_name": BOND_STEREO_INT_TO_NAME,
        "bond_dir_int_to_name": BOND_DIR_INT_TO_NAME,
        "bond_stereo_counts": dict(bond_stereo_counts.most_common()),
        "bond_stereo_int_counts": {str(k): v for k, v in bond_stereo_int_counts.most_common()},
        "bond_dir_counts": dict(bond_dir_counts.most_common()),
        "bond_dir_int_counts": {str(k): v for k, v in bond_dir_int_counts.most_common()},
        "bond_cip_counts": dict(bond_cip_counts.most_common()),
        "chiral_tag_counts": dict(chiral_tag_counts.most_common()),
        "chiral_tag_int_counts": {str(k): v for k, v in chiral_tag_int_counts.most_common()},
        "atom_cip_counts": dict(atom_cip_counts.most_common()),
        "bond_stereo_to_int_mapping": {k: sorted(v) for k, v in bond_stereo_to_int_mapping.items()},
        "bond_dir_to_int_mapping": {k: sorted(v) for k, v in bond_dir_to_int_mapping.items()},
        "bond_stereo_examples": dict(bond_stereo_examples),
        "bond_dir_examples": dict(bond_dir_examples),
        "bond_cip_examples": dict(bond_cip_examples),
        "chiral_tag_examples": dict(chiral_tag_examples),
        "atom_cip_examples": dict(atom_cip_examples),
        "smiles_to_datasets": dict(smiles_to_datasets),
    }

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n\nResults saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze bond stereochemistry across datasets.")
    parser.add_argument("--data-dir", type=Path, default=_DEFAULT_DATA_DIR,
                        help="Directory containing CSV files (default: relative to script)")
    parser.add_argument("--output", type=Path, default=_DEFAULT_OUTPUT,
                        help="Output JSON path (default: ../results/bond_stereo_analysis_results.json)")
    parser.add_argument("--single-file", type=str, default=None,
                        help="Process only this CSV file (e.g., 'tox21.csv'). Default: all CSVs.")
    args = parser.parse_args()

    process_datasets(data_dir=args.data_dir, output_path=args.output, single_file=args.single_file)
