"""
Analyze atom chirality across chemprop v1 datasets.

For each molecule (canonical SMILES):
- Count ChiralTag values and int(ChiralTag) values across all atoms
- Count CIPCode (R/S/r/s) values across all atoms
- Store up to 5 example molecules per ChiralTag and CIPCode
- Store the ChiralTag -> int mapping for every observed case
- Track which dataset each canonical SMILES comes from

Optimization: only fully process molecules whose SMILES contains '@'
(non-@ molecules cannot have stereocenters).

Usage:
    conda activate chemprop
    python analyze_atom_chirality_datasets.py
    python analyze_atom_chirality_datasets.py --data-dir /path/to/data --output ../results/atom_chirality_analysis_results.json
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

# Suppress RDKit warnings (e.g., valence errors for Al)
RDLogger.DisableLog("rdApp.*")

# ============================================================
# Configuration
# ============================================================
# Download data from here: https://github.com/kspieks/chemprop/blob/barrier_prediction/data.tar.gz
# then unzip and update _DEFAULT_DATA_DIR

# Defaults are relative to this script's location (reference/scripts/)
_SCRIPT_DIR = Path(__file__).resolve().parent
_DEFAULT_DATA_DIR = _SCRIPT_DIR / "../../data"
_DEFAULT_OUTPUT = _SCRIPT_DIR / "../results/atom_chirality_analysis_results.json"

# Set to a filename to test on one file, or None to process all files
# (now controlled via --single-file CLI arg)

MAX_EXAMPLES = 5  # Max example molecules stored per category

# Reference mapping (from RDKit ChiralType enum)
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


def get_smiles_column(header):
    """Determine which column contains SMILES."""
    header_lower = [h.lower().strip() for h in header]
    if "smiles" in header_lower:
        return header_lower.index("smiles")
    elif "mol" in header_lower:
        return header_lower.index("mol")
    else:
        return 0


def add_example(examples_dict, key, canon_smi):
    """Add canon_smi to examples_dict[key] if under MAX_EXAMPLES and not duplicate."""
    if len(examples_dict[key]) < MAX_EXAMPLES:
        if canon_smi not in examples_dict[key]:
            examples_dict[key].append(canon_smi)


def process_datasets(data_dir: Path = _DEFAULT_DATA_DIR, output_path: Path = _DEFAULT_OUTPUT, single_file: str | None = None):
    # ============================================================
    # Data structures
    # ============================================================
    # Value counts (only for chiral atoms, i.e. tag != 0)
    chiral_tag_counts = Counter()        # str(tag) -> count
    chiral_tag_int_counts = Counter()    # int(tag) -> count
    cip_code_counts = Counter()          # CIPCode string -> count

    # Examples: up to MAX_EXAMPLES canonical SMILES per category
    chiral_tag_examples = defaultdict(list)  # str(tag) -> [smiles, ...]
    cip_code_examples = defaultdict(list)    # CIPCode -> [smiles, ...]

    # Mapping verification: str(tag) -> set of int values observed
    tag_to_int_mapping = defaultdict(set)

    # Dataset tracking: canonical_smiles -> list of dataset names
    # (only for molecules with chirality)
    smiles_to_datasets = defaultdict(list)

    # Per-dataset stats
    dataset_stats = {}

    # Global stats
    total_molecules = 0
    total_chiral_molecules = 0
    total_skipped_no_at = 0
    failed_parses = 0

    # ============================================================
    # Determine which files to process
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
        file_chiral_count = 0
        file_fail_count = 0
        file_skipped_no_at = 0

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

                # Optimization: skip molecules without @ (no stereocenters possible)
                if "@" not in raw_smi:
                    file_skipped_no_at += 1
                    total_skipped_no_at += 1
                    continue

                # Parse and canonicalize
                mol = Chem.MolFromSmiles(raw_smi)
                if mol is None:
                    file_fail_count += 1
                    failed_parses += 1
                    continue

                canon_smi = Chem.MolToSmiles(mol)

                # Re-parse from canonical SMILES for consistent analysis
                mol2 = Chem.MolFromSmiles(canon_smi)
                rdCIPLabeler.AssignCIPLabels(mol2)

                has_chirality = False
                for atom in mol2.GetAtoms():
                    tag = atom.GetChiralTag()
                    tag_str = str(tag)
                    tag_int = int(tag)

                    if tag_int != 0:  # Has some chirality
                        has_chirality = True

                        # Count
                        chiral_tag_counts[tag_str] += 1
                        chiral_tag_int_counts[tag_int] += 1

                        # Mapping verification
                        tag_to_int_mapping[tag_str].add(tag_int)

                        # Examples
                        add_example(chiral_tag_examples, tag_str, canon_smi)

                    # CIPCode is assigned only to atoms with a valid ChiralTag
                    # (empirically confirmed: ChiralTag count == CIPCode count)
                    if atom.HasProp("_CIPCode"):
                        cip = atom.GetProp("_CIPCode")   # R, S, r, or s
                        cip_code_counts[cip] += 1
                        add_example(cip_code_examples, cip, canon_smi)

                if has_chirality:
                    file_chiral_count += 1
                    total_chiral_molecules += 1
                    smiles_to_datasets[canon_smi].append(dataset_name)

        file_elapsed = time.time() - file_start
        avg_ms = (file_elapsed / file_mol_count * 1000) if file_mol_count > 0 else 0
        processed_count = file_mol_count - file_skipped_no_at

        dataset_stats[dataset_name] = {
            "total_molecules": file_mol_count,
            "skipped_no_at": file_skipped_no_at,
            "processed_with_rdkit": processed_count,
            "failed_parses": file_fail_count,
            "molecules_with_chirality": file_chiral_count,
            "pct_chiral": (file_chiral_count / file_mol_count * 100) if file_mol_count > 0 else 0,
            "time_seconds": round(file_elapsed, 2),
            "avg_ms_per_molecule": round(avg_ms, 3),
        }

        print(
            f"  {csv_path.name:<25s} | "
            f"{file_mol_count:>7d} mols | "
            f"{file_skipped_no_at:>7d} skipped (no @) | "
            f"{processed_count:>7d} parsed | "
            f"{file_chiral_count:>6d} chiral ({dataset_stats[dataset_name]['pct_chiral']:.1f}%) | "
            f"{file_elapsed:.2f}s ({avg_ms:.3f} ms/mol)"
        )

    overall_elapsed = time.time() - overall_start

    # ============================================================
    # Report Results (printed to stdout; tee to .log file for persistence)
    # ============================================================
    print(f"\n{'=' * 80}")
    print("RESULTS")
    print(f"{'=' * 80}")

    print(f"\nTotal molecules across all datasets: {total_molecules}")
    print(f"Skipped (no @ in SMILES):            {total_skipped_no_at}")
    print(f"Processed with RDKit:                {total_molecules - total_skipped_no_at}")
    print(f"Failed parses:                       {failed_parses}")
    print(f"Molecules with chirality:            {total_chiral_molecules}")
    print(f"Unique chiral canonical SMILES:      {len(smiles_to_datasets)}")
    print(f"\nTotal time: {overall_elapsed:.2f}s")
    print(f"Avg time per molecule (all): {overall_elapsed / total_molecules * 1000:.3f} ms/mol")
    processed_total = total_molecules - total_skipped_no_at
    if processed_total > 0:
        print(f"Avg time per RDKit-processed molecule: {overall_elapsed / processed_total * 1000:.3f} ms/mol")

    # Datasets with no stereocenters
    no_stereo_datasets = [name for name, stats in dataset_stats.items() if stats["molecules_with_chirality"] == 0]
    if no_stereo_datasets:
        print(f"\nDatasets with NO stereocenters at all:")
        for name in no_stereo_datasets:
            print(f"  - {name} ({dataset_stats[name]['total_molecules']} molecules)")

    # ChiralTag counts
    print(f"\n{'─' * 80}")
    print("a.GetChiralTag() value counts (atoms with chirality):")
    print(f"{'─' * 80}")
    for tag_str, count in chiral_tag_counts.most_common():
        print(f"  {tag_str:<35s}: {count:>10d}")

    # int(ChiralTag) counts
    print(f"\n{'─' * 80}")
    print("int(a.GetChiralTag()) value counts (atoms with chirality):")
    print(f"{'─' * 80}")
    for tag_int, count in chiral_tag_int_counts.most_common():
        print(f"  {tag_int}: {count:>10d}")

    # CIP code counts
    print(f"\n{'─' * 80}")
    print("atom.GetProp('_CIPCode') value counts:")
    print(f"{'─' * 80}")
    for cip, count in cip_code_counts.most_common():
        print(f"  {cip}: {count:>10d}")

    # Mapping verification
    print(f"\n{'─' * 80}")
    print("ChiralTag str -> int mapping (data-driven verification):")
    print(f"{'─' * 80}")
    for tag_str, int_vals in sorted(tag_to_int_mapping.items(), key=lambda x: min(x[1])):
        print(f"  {tag_str:<35s} -> int values observed: {sorted(int_vals)}")

    # Examples
    print(f"\n{'─' * 80}")
    print(f"Example molecules per ChiralTag (up to {MAX_EXAMPLES}):")
    print(f"{'─' * 80}")
    for tag_str in sorted(chiral_tag_examples.keys()):
        print(f"\n  {tag_str}:")
        for smi in chiral_tag_examples[tag_str]:
            print(f"    {smi}")

    print(f"\n{'─' * 80}")
    print(f"Example molecules per CIPCode (up to {MAX_EXAMPLES}):")
    print(f"{'─' * 80}")
    for cip in sorted(cip_code_examples.keys()):
        print(f"\n  {cip}:")
        for smi in cip_code_examples[cip]:
            print(f"    {smi}")

    # Per-dataset summary table
    print(f"\n{'─' * 80}")
    print("Per-dataset summary:")
    print(f"{'─' * 80}")
    print(f"  {'Dataset':<20s} | {'Total':>7s} | {'No @':>7s} | {'Parsed':>7s} | {'Chiral':>7s} | {'%Chiral':>7s} | {'Time':>6s}")
    print(f"  {'─' * 20}-+-{'─' * 7}-+-{'─' * 7}-+-{'─' * 7}-+-{'─' * 7}-+-{'─' * 7}-+-{'─' * 6}")
    for name in sorted(dataset_stats.keys()):
        s = dataset_stats[name]
        print(
            f"  {name:<20s} | {s['total_molecules']:>7d} | {s['skipped_no_at']:>7d} | "
            f"{s['processed_with_rdkit']:>7d} | {s['molecules_with_chirality']:>7d} | "
            f"{s['pct_chiral']:>6.1f}% | {s['time_seconds']:>5.1f}s"
        )

    # ============================================================
    # Save structured results to JSON for later inspection
    # ============================================================
    results = {
        "config": {
            "data_dir": str(data_dir),
            "files_processed": [p.name for p in csv_files],
            "single_file_mode": single_file,
            "at_filter_optimization": True,
        },
        "timing": {
            "total_seconds": round(overall_elapsed, 2),
            "avg_ms_per_molecule_all": round(overall_elapsed / total_molecules * 1000, 3) if total_molecules > 0 else 0,
            "avg_ms_per_rdkit_processed": round(overall_elapsed / processed_total * 1000, 3) if processed_total > 0 else 0,
        },
        "stats": {
            "total_molecules": total_molecules,
            "skipped_no_at": total_skipped_no_at,
            "processed_with_rdkit": processed_total,
            "failed_parses": failed_parses,
            "molecules_with_chirality": total_chiral_molecules,
            "unique_chiral_canonical_smiles": len(smiles_to_datasets),
        },
        "datasets_with_no_stereocenters": no_stereo_datasets,
        "dataset_stats": dataset_stats,
        "chiral_tag_counts": dict(chiral_tag_counts.most_common()),
        "chiral_tag_int_counts": {str(k): v for k, v in chiral_tag_int_counts.most_common()},
        "cip_code_counts": dict(cip_code_counts.most_common()),
        "tag_to_int_mapping": {k: sorted(v) for k, v in tag_to_int_mapping.items()},
        "chiral_tag_int_to_name": CHIRAL_TAG_INT_TO_NAME,
        "chiral_tag_examples": dict(chiral_tag_examples),
        "cip_code_examples": dict(cip_code_examples),
        "smiles_to_datasets": dict(smiles_to_datasets),
    }

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n\nResults saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze atom chirality across datasets.")
    parser.add_argument("--data-dir", type=Path, default=_DEFAULT_DATA_DIR,
                        help="Directory containing CSV files (default: relative to script)")
    parser.add_argument("--output", type=Path, default=_DEFAULT_OUTPUT,
                        help="Output JSON path (default: ../results/atom_chirality_analysis_results.json)")
    parser.add_argument("--single-file", type=str, default=None,
                        help="Process only this CSV file (e.g., 'tox21.csv'). Default: all CSVs.")
    args = parser.parse_args()

    process_datasets(data_dir=args.data_dir, output_path=args.output, single_file=args.single_file)
