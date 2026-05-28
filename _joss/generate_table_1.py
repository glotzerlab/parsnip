r"""Measure the accuracy of various CIF parsers using a consensus-based system.

To reproduce the figure from the paper, download the files from `cifs.txt`. The
following command will download the required files to a new parsnip/cifs directory:

```bash
> mkdir -p cifs && sed 's|^|http://www.crystallography.net/cod/|' \
  _joss/table_1_cifs.txt | wget -i - -P ./cifs

> python ./_joss/generate_table_1.py
```
"""

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import gemmi
from ase.io.cif import read_cif
from pymatgen.io.cif import CifParser
from tqdm import tqdm

import parsnip

COD_FILES = list(Path("./cifs").glob("*.cif"))

# Filter out ASE setting warnings, which do not affect the results.
warnings.filterwarnings(action="ignore", category=UserWarning, message="crystal system")
warnings.filterwarnings(
    action="ignore", category=UserWarning, message="scaled_positions"
)
warnings.filterwarnings(action="ignore", category=UserWarning)
warnings.filterwarnings(action="ignore", category=parsnip._errors.ParseWarning)


@dataclass
class Parser:
    """A parser with its function and tracking lists."""

    name: str
    func: Callable
    mismatches: list = field(default_factory=list)
    errors: list = field(default_factory=list)


def count_atoms_ase(cif_path: Path | str) -> int:
    """Count atoms using ASE's read_cif with full unit cell."""
    with open(cif_path) as f:
        atoms = read_cif(f, index=None, reader="ase")
    return len(atoms)


def count_atoms_gemmi(cif_path: Path | str) -> int:
    """Count atoms using gemmi's read_small_structure with full unit cell."""
    # st = gemmi.read_small_structure(str(cif_path))
    st = gemmi.make_small_structure_from_block(
        gemmi.cif.read_file(str(cif_path), check_level=0).sole_block()
    )
    return len(st.get_all_unit_cell_sites())


def count_atoms_parsnip_rational_3(cif_path: Path | str) -> int:
    """Count atoms using parsnip's build_unit_cell."""
    cif = parsnip.CifFile(str(cif_path))
    unit_cell = cif.build_unit_cell(n_decimal_places=3, parse_mode="rational")
    return len(unit_cell)


def count_atoms_pymatgen(cif_path: Path | str) -> int:
    """Count atoms using pymatgen's CifParser."""
    parser = CifParser(str(cif_path), check_cif=False)
    structure = parser.parse_structures(primitive=False)[0]
    return len(structure)


def get_cell_formula_units_z(cif_path: str) -> int | None:
    """Extract the _cell_formula_units_Z value from CIF.

    Returns Z if present and not 1, otherwise None.
    """
    try:
        z_line = parsnip.CifFile(cif_path)["_cell_formula_units_Z"]
        if z_line:
            match = re.search(r"(\d+)", z_line)
            if match:
                z = int(match.group(1))
                if z != 1:
                    return z
    except Exception:
        return None
    return None


# def fetch_cod_files(cod_dir: Path, limit: int = 1000) -> list[Path]:
#     """Fetch CIF files from the COD database.
#
#     NOTE: This was written before I realized there was an rsync endpoint. Using that
#     is much preferred to SQL. See https://wiki.crystallography.net/howtoobtaincod/

#     Args:
#         cod_dir: Directory to store CIF files
#         limit: Maximum number of files to fetch

#     Returns
#     -------
#         List of paths to downloaded CIF files
#     """
#     import traceback
#     from urllib.request import urlopen
#     import mysql.connector

#     cod_dir.mkdir(exist_ok=True)

#     # Connect to COD database
#     cnx = mysql.connector.connect(
#         user="cod_reader", host="www.crystallography.net", database="cod"
#     )
#     cursor = cnx.cursor()

#     query = f"SELECT file, formula FROM data WHERE 1=1 LIMIT {limit}"
#     cursor.execute(query)

#     downloaded = []
#     base_url = "http://www.crystallography.net/cod/"

#     for file_id, _ in tqdm(cursor, desc="Fetching COD files"):
#         url = f"{base_url}{file_id}.cif"
#         output_path = cod_dir / f"{file_id}.cif"

#         if output_path.exists():
#             downloaded.append(output_path)
#             continue

#         try:
#             with urlopen(url) as response:
#                 content = response.read().decode("utf-8")
#                 with open(output_path, "w") as f:
#                     f.write(content)
#             downloaded.append(output_path)
#         except Exception as e:
#             print(f"Failed to download {file_id}: {e}")

#     cursor.close()
#     cnx.close()
#     return downloaded


def print_summary_table(parsers: list[Parser], num_files: int):
    """Print a summary table for a list of parsers."""
    print(f"Checked {num_files} CIF files\n")

    for parser in parsers:
        if parser.mismatches:
            print(f"{parser.name.upper()} MISMATCHES ({len(parser.mismatches)}):")
            print("-" * 60)
            for name, expected, actual, is_hr in parser.mismatches[:10]:
                hr_note = " (hR, 3x applied)" if is_hr else ""
                print(f"{name}: expected {expected}, got {actual}{hr_note}")
            if len(parser.mismatches) > 10:
                print(f"... and {len(parser.mismatches) - 10} more")
            print()

        if parser.errors:
            print(f"{parser.name.upper()} ERRORS ({len(parser.errors)}):")
            print("-" * 60)
            for name, msg in parser.errors[:10]:  # Show first 10 errors per parser
                print(f"{name}: {msg}")
            if len(parser.errors) > 10:
                print(f"... and {len(parser.errors) - 10} more")
            print()

    # Summary table
    print("=" * 95)
    print(
        f"{'Library':<22} | {'Correct Crystals':<18} {'Incorrect Crystals':<18}"
        f" {'Failed to Parse':<12} {'Correct %':>15}"
    )
    print("-" * 95)

    total_correct = 0
    total_incorrect = 0
    total_failed = 0

    for parser in parsers:
        correct_crystals = num_files - len(parser.errors) - len(parser.mismatches)
        incorrect_crystals = len(parser.mismatches)
        failed = len(parser.errors)
        pct = (correct_crystals / num_files * 100) if num_files > 0 else 0

        print(
            f"{parser.name:<22} | {correct_crystals:<18} {incorrect_crystals:<18} "
            f"{failed:<12} {pct:>14.1f}%"
        )

        total_correct += correct_crystals
        total_incorrect += incorrect_crystals
        total_failed += failed

    total_pct = (
        (total_correct / (num_files * len(parsers)) * 100) if num_files > 0 else 0
    )

    print("-" * 95)
    print(
        f"{'Total':<22} | {total_correct:<12} {total_incorrect:<12} "
        f"{total_failed:<12} {total_pct:>14.1f}%"
    )
    print("=" * 95)


if __name__ == "__main__":
    # Optionally fetch more COD files
    # _ = fetch_cod_files(cod_dir, limit=150_000)
    # cod_files = list(cod_dir.glob("*.cif"))

    cod_parsers = [
        Parser("parsnip (rational, 3)", count_atoms_parsnip_rational_3),
        Parser("ASE", count_atoms_ase),
        Parser("gemmi", count_atoms_gemmi),
        Parser("pymatgen", count_atoms_pymatgen),
    ]

    # Process COD files, using a consensus between parsers as the correct output.
    # Note that, in cases where structures may have multiple formula units, we allow for
    # the replicated (_cell_formula_units_Z>=1) or minimal structure (both are valid).
    cod_file_count = 0
    for cif_file in tqdm(sorted(COD_FILES)[::-1], desc="Processing COD files"):
        cod_file_count += 1
        counts = {}
        for parser in cod_parsers:
            try:
                counts[parser.name] = parser.func(cif_file)
            except Exception as e:
                parser.errors.append((cif_file.name, f"{parser.name} read error: {e}"))

        if not counts:
            continue

        z = get_cell_formula_units_z(cif_file)

        # Use the most common value between parsers as consensus
        most_common = max(set(counts.values()), key=list(counts.values()).count)

        for parser in cod_parsers:
            if parser.name not in counts:
                continue

            if z is not None and z > 1:
                # Either consensus OR consensus/Z could be correct
                consensus_div_z = most_common // z
                is_valid = (
                    counts[parser.name] == most_common
                    or counts[parser.name] == consensus_div_z
                )
                if not is_valid:
                    parser.mismatches.append(
                        (
                            cif_file.name,
                            f"{most_common} or {consensus_div_z}",
                            counts[parser.name],
                            False,
                        )
                    )
            elif counts[parser.name] != most_common:
                parser.mismatches.append(
                    (
                        cif_file.name,
                        most_common,
                        counts[parser.name],
                        False,
                    )
                )

    print("\n" + "=" * 70)
    print("COD FILES (ref = consensus)")
    print("=" * 70)
    print_summary_table(cod_parsers, len(COD_FILES))
