# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

r"""Measure the space group symmetry accuracy of various CIF parsers.

For each parser, reconstruct the crystal structure from the COD subset and use
spglib to determine the space group. This compares against the reference space group
number stored in the CIF metadata and categorizes results as correct, too high
(spurious symmetry), or too low (lost symmetry).

Requires ``spglib`` and ``matplotlib`` in addition to the dependencies in
``generate_table_1.py``.

To reproduce the figure from the paper, first download the CIF files as described
in ``generate_table_1.py``, then run:

.. code-block:: bash

    python ./_joss/generate_table_2.py
"""

import re
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np
import spglib
from ase.data import atomic_numbers
from ase.io.cif import read_cif
from pymatgen.io.cif import CifParser
from tqdm import tqdm

import parsnip
from parsnip.patterns import _matrix_from_lengths_and_angles

COD_FILES = list(Path("./cifs").glob("*.cif"))
SYMPREC = 1e-3

# Filter out ASE setting warnings, which do not affect the results.
warnings.filterwarnings(action="ignore", category=UserWarning)
warnings.filterwarnings(action="ignore", category=parsnip._errors.ParseWarning)


@dataclass
class SpaceGroupParser:
    """Track space group accuracy for a single parser."""

    name: str
    func: Callable
    correct: int = 0
    too_high: int = 0
    too_low: int = 0
    errors: list = field(default_factory=list)


def get_reference_spacegroup(cif_path: str) -> int | None:
    """Extract the International Tables space group number from CIF metadata.

    Uses parsnip wildcard queries to match across naming conventions.
    """
    try:
        cif = parsnip.CifFile(cif_path)
        for pattern in ("_space_group_*", "_symmetry_Int*"):
            val = cif[pattern]
            if val is None:
                continue
            for v in val if isinstance(val, list) else [val]:
                match = re.match(r"^\d+$", str(v).strip())
                if match:
                    return int(match.group(0))
    except Exception as e:
        print(e)
    return None


def _extract_element(symbol: str) -> str:
    """Extract an element symbol from a type symbol or site label.

    Handles charge states like 'Fe3+' and numeric suffixes like 'Cu1'.
    """
    m = re.match(r"^([A-Z][a-z]?)", symbol.strip())
    return m.group(1) if m else ""


def get_cell_parsnip(cif_path):
    """Build an spglib cell tuple from parsnip."""
    cif = parsnip.CifFile(str(cif_path))
    lattice = cif.lattice_vectors.T

    # Use _atom_site_type_symbol when available (most reliable element identification),
    # fall back to regex extraction from _atom_site_label.
    type_sym = cif.get_from_loops("_atom_site?type_symbol")
    if type_sym is not None and len(type_sym) > 0:
        aux_data, positions = cif.build_unit_cell(
            n_decimal_places=3,
            parse_mode="rational",
            additional_columns=["_atom_site_type_symbol"],
        )
    else:
        aux_data, positions = cif.build_unit_cell(
            n_decimal_places=3,
            parse_mode="rational",
            additional_columns=["_atom_site_label"],
        )
    numbers = np.array(
        [atomic_numbers.get(_extract_element(str(s)), 0) for s in aux_data.flatten()]
    )
    return _validate_cell(lattice, positions, numbers)


def get_cell_ase(cif_path):
    """Build an spglib cell tuple from ASE."""
    with open(cif_path) as f:
        atoms = read_cif(f, index=None, reader="ase")
    return _validate_cell(
        atoms.cell.array, atoms.get_scaled_positions(), atoms.get_atomic_numbers()
    )


def _validate_cell(lattice, positions, numbers):
    """Return None if the cell has NaN/inf or zero-length lattice vectors."""
    if (
        np.any(~np.isfinite(lattice))
        or np.any(~np.isfinite(positions))
        or np.any(~np.isfinite(numbers))
        or np.linalg.det(lattice) == 0
        or len(positions) == 0
    ):
        return None
    return (
        np.asarray(lattice, dtype=float),
        np.asarray(positions, dtype=float),
        np.asarray(numbers, dtype=int),
    )


def get_cell_gemmi(cif_path):
    """Build an spglib cell tuple from gemmi."""
    block = gemmi.cif.read_file(str(cif_path), check_level=0).sole_block()
    st = gemmi.make_small_structure_from_block(block)
    cell = st.cell
    lattice = _matrix_from_lengths_and_angles(
        cell.a,
        cell.b,
        cell.c,
        np.radians(cell.alpha),
        np.radians(cell.beta),
        np.radians(cell.gamma),
    )
    sites = st.get_all_unit_cell_sites()
    positions = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in sites])
    numbers = np.array([s.element.atomic_number for s in sites])
    return _validate_cell(lattice, positions, numbers)


def get_cell_pymatgen(cif_path):
    """Build an spglib cell tuple from pymatgen."""
    parser = CifParser(str(cif_path), check_cif=False)
    structure = parser.parse_structures(primitive=False)[0]
    numbers = []
    for site in structure:
        try:
            numbers.append(site.specie.Z)
        except AttributeError:
            numbers.append(next(iter(site.species.keys())).Z)
    return _validate_cell(
        structure.lattice.matrix, structure.frac_coords, np.array(numbers)
    )


if __name__ == "__main__":
    parsers = [
        SpaceGroupParser("parsnip (rational, 3)", get_cell_parsnip),
        SpaceGroupParser("ASE", get_cell_ase),
        SpaceGroupParser("gemmi", get_cell_gemmi),
        SpaceGroupParser("pymatgen", get_cell_pymatgen),
    ]

    skipped = 0
    for cif_file in tqdm(sorted(COD_FILES), desc="Processing COD files"):
        ref_sg = get_reference_spacegroup(str(cif_file))
        if ref_sg is None:
            skipped += 1
            continue

        for parser in parsers:
            try:
                cell = parser.func(cif_file)
                if cell is None:
                    parser.errors.append((cif_file.name, "invalid_cell"))
                    continue
                dataset = spglib.get_symmetry_dataset(cell, symprec=SYMPREC)
                if dataset is None:
                    parser.errors.append((cif_file.name, "spglib_none"))
                    continue

                sg = dataset.number
                if sg == ref_sg:
                    parser.correct += 1
                elif sg > ref_sg:
                    parser.too_high += 1
                else:
                    parser.too_low += 1
            except Exception as e:
                parser.errors.append((cif_file.name, str(e)))

    total = len(COD_FILES) - skipped
    print(
        f"\nChecked {total} CIF files ({skipped} skipped, no space group in metadata)"
    )
    print(f"symprec = {SYMPREC}\n")
    print(
        f"{'Library':<25} {'Correct':>8} {'Too High':>10} {'Too Low':>9} {'Errors':>8}"
    )
    print("-" * 65)
    for p in parsers:
        print(
            f"{p.name:<25} {p.correct:>8} {p.too_high:>10} "
            f"{p.too_low:>9} {len(p.errors):>8}"
        )

    # Print error breakdown per parser
    for p in parsers:
        if p.errors:
            invalid = sum(1 for _, reason in p.errors if reason == "invalid_cell")
            spglib_fail = sum(1 for _, reason in p.errors if reason == "spglib_none")
            exceptions = [
                (f, r) for f, r in p.errors if r not in ("invalid_cell", "spglib_none")
            ]
            print(
                f"\n{p.name} errors ({len(p.errors)}): "
                f"{invalid} invalid_cell, {spglib_fail} spglib_none, "
                f"{len(exceptions)} exceptions"
            )
            for f, r in exceptions[:5]:
                print(f"  {f}: {r[:100]}")
            if len(exceptions) > 5:
                print(f"  ... and {len(exceptions) - 5} more exceptions")
