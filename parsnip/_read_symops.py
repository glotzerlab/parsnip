# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Parse and process the OpenBabel space groups database.

Similar to *gemmi*, we use this as a complete reference for the 530+ standard space
groups and settings. This data is processed into json, which is used when CIF files do
not store their own symmetry operations. Note that we cross-reference against [Table 6
from this site](https://cci.lbl.gov/sginfo/hall_symbols.html) to determine the standard
setting

https://github.com/openbabel/openbabel/blob/889c350feb179b43aa43985799910149d4eaa2bc/data/space-groups.txt
"""

import re


def build_sg_dict(openbabel_db_path, crossref_db_path):
    r"""Process the OpenBabel space group dict into a JSON.

    As of writing this, the OB database has 541 settings+sgs. There are 530 standard
    settings and infinitely many possible options, so this is a happy medium.

    This database stores entries in the following format:
    ```
    <InternationalTablesNumber>
    <HallSymbol>
    <ExtendedHermannMauguinSymbol>
    (<Symop>\n)+
    \n
    ```

    We convert this to a JSON with the following format:
    ```
    "<HallSymbol>" : {
        "IT": <InternationalTablesNumber>,
        "HM": <ExtendedHermannMauguinSymbol>,
        "symops": [ <Symop>(, <Symop>)* ],
        "is_default_setting": (true)|(false)
    }
    ```
    Note that the cif key `_space_group_name_H-M_alt` does support both regular and
    extended HM notation, but `_symmetry_space_group_name_H-M` only supports the
    standard variant.
    """
    spacegroup_data = {}
    default_settings_it_to_hall = {}
    hall_to_it_setting = {}
    with open(openbabel_db_path) as f:
        sg_chunks = re.split(r"\n\n+", f.read().rstrip())
    with open(crossref_db_path) as f:
        ref_chunks = [re.findall(r"\S+(?: \S+)*", row) for row in f][
            2:
        ]  # skip 2 header rows
        for it_with_setting, _, hall in ref_chunks:
            # We rely on the fact that the default centering is first in the IT.
            it = it_with_setting.split(":")[0]
            if not default_settings_it_to_hall.get(it):
                default_settings_it_to_hall[it] = hall
            hall_to_it_setting[hall] = it_with_setting

    for entry in sg_chunks:
        parsed_entry = entry.split("\n")
        it, hall, ehm = parsed_entry[:3]
        # if "," in ehm: # Short name != full name
        hm_short, hm_full = ehm.split(",") if "," in ehm else (ehm, ehm)
        symops = parsed_entry[3:]
        spacegroup_data[hall] = {
            "table_number": it,
            "tables_number_with_setting": hall_to_it_setting.get(hall),
            "hermann_mauguin_short": hm_short,
            "hermann_mauguin_full": hm_full,
            "is_default_setting": default_settings_it_to_hall[it] == hall,
            "symops": symops,
        }

    return spacegroup_data


if __name__ == "__main__":
    import json

    OUTPUT_FN = "parsnip/symops.json"
    INPUT_FN = "space-groups.txt"

    data = build_sg_dict(INPUT_FN, "temp.txt")

    # with open(OUTPUT_FN, "w") as f:
    #     json.dump(data, f, indent="  ")

    with open(OUTPUT_FN) as f:
        json.load(f)

    # Properties to test
    # [ ] only one default setting per sg
    # [ ] all sgs represented
    # [ ] short symbols are correct? may be hard
    # [ ] number of symops is consistent
    # [x] Queries reconstruct the original inputs
