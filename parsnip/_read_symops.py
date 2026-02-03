# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

import json
import re
from fractions import Fraction

import requests
import tqdm
from bs4 import BeautifulSoup

# Standard centering vectors from International Tables
CENTERING_VECTORS = {
    "P": [(Fraction(0), Fraction(0), Fraction(0))],  # Primitive
    "I": [
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction("1/2"), Fraction("1/2"), Fraction("1/2")),
    ],  # Body-centered
    "F": [
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction("1/2"), Fraction("1/2")),
        (Fraction("1/2"), Fraction(0), Fraction("1/2")),
        (Fraction("1/2"), Fraction("1/2"), Fraction(0)),
    ],  # Face-centered
    "A": [
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction("1/2"), Fraction("1/2")),
    ],  # A-face centered
    "B": [
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction("1/2"), Fraction(0), Fraction("1/2")),
    ],  # B-face centered
    "C": [
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction("1/2"), Fraction("1/2"), Fraction(0)),
    ],  # C-face centered
    "R": [
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction("2/3"), Fraction("1/3"), Fraction("1/3")),
        (Fraction("1/3"), Fraction("2/3"), Fraction("2/3")),
    ],  # Rhombohedral (hex setting)
}


def add_shift(op_comp, shift):
    """Add a shift (Fraction) to a symmetry component string (e.g. 'x+1/2', '1/2-y')."""
    if shift == 0:
        return op_comp

    # Parse existing component into (variable_part, constant_part)
    # [+-]?            Optional sign
    # (?:              Start non-capturing group for value types
    #   \d+(?:/\d+)?   Fraction or integer
    #   | [xyz]        Variable
    # )
    terms = re.findall(r"[+-]?(?:(?:\d+(?:/\d+)?)|[xyz])", op_comp.strip())

    var_str, current_frac = "", Fraction(0)

    # Accumulate the fractional parts and store the variable
    for term in terms:
        if any(c in term for c in "xyz"):
            var_str = term
        else:
            current_frac += Fraction(term)

    # Shift based on our centering, then wrap into [0, 1)
    final_frac = (current_frac + shift) % Fraction(1)

    # Reconstruct the string
    # If const is 0, just return variable. Otherwise, return shift ± xyz.
    if final_frac == 0:
        return var_str.lstrip("+")

    frac_str = str(final_frac)

    if var_str.startswith("-"):
        return f"{frac_str}{var_str}"

    # Variable is positive (e.g. "x" or "+x")
    clean_var = var_str.lstrip("+")
    if final_frac > 0:
        return f"{frac_str}+{clean_var}"
    # const_str already has sign (e.g. "-1/4")
    return f"{frac_str}{var_str}"


def expand_symmetry(scraped_ops, lattice_type):
    """Expand symmetry operations based on lattice type."""
    vectors = CENTERING_VECTORS.get(lattice_type, CENTERING_VECTORS["P"])
    full_set = []

    for vec in vectors:
        vx, vy, vz = vec
        # If we are the identity vector, we can just add everything at once
        if vx == 0 and vy == 0 and vz == 0:
            full_set.extend(scraped_ops)
            continue

        for op in scraped_ops:
            # "x,y,z" -> ["x", "y", "z"]
            comps = op.split(",")

            new_comps = [
                add_shift(comps[0], vx),
                add_shift(comps[1], vy),
                add_shift(comps[2], vz),
            ]
            full_set.append(",".join(new_comps))

    return full_set


if __name__ == "__main__":

    def build_sg_dict():
        """Extract symops from the hypertext crystallography book."""
        data = {}
        for i in tqdm.tqdm(range(1, 231)):
            # c-unique centering where applicable
            url = f"http://img.chem.ucl.ac.uk/sgp/LARGE/{i:03d}az3.htm"
            response = requests.get(url, timeout=10)
            soup = BeautifulSoup(response.content, "html.parser")

            # Extract raw symops
            text = soup.find("font", {"face": "COURIER"}).get_text()
            symops = [
                line.replace(" ", "") for line in text.splitlines() if "," in line
            ]

            # Extract Hermann-Mauguin notation
            header_font = soup.find("font", {"size": "+2"})
            if header_font:
                header_text = header_font.get_text(" ", strip=True)
                # Normalize whitespace
                header_text = " ".join(header_text.split()).removeprefix("Space Group")
                hm_symbol = header_text.strip()
            else:
                raise ValueError("Did not find HM symbol!")

            # Determine lattice type from the first character of the HM symbol
            lattice_type = hm_symbol[0]

            expanded_symops = expand_symmetry(symops, lattice_type)
            print(i, expanded_symops)

            # { "1": { "xyz": [...], "HM": "P 1" } }
            data[str(i)] = {"xyz": expanded_symops, "HM": hm_symbol}

        return data

    OUTPUT_FN = "parsnip/symops.json"

    data = build_sg_dict()

    with open(OUTPUT_FN, "w") as f:
        json.dump(data, f, indent="  ")
    with open(OUTPUT_FN) as f:
        data = json.load(f)
    print(len(data["230"]["xyz"]))
    print(len(data["5"]["xyz"]))
