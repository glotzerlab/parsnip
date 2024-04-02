import numpy as np
from conftest import cif_files_mark
from gemmi import cif

from parsnip.parse import read_table


@cif_files_mark
def test_read_symop(cif_data):
    parsnip_data = read_table(filename=cif_data.filename, keys=cif_data.symop_keys)

    gemmi_data = np.array(
        cif.read_file(cif_data.filename).sole_block().find(cif_data.symop_keys)
    )

    # We replace ", " strings with "," to ensure data is collected properly
    # We have to apply this same transformation to the gemmi data to check correctness.
    if "CCDC_1446529_Pm-3m.cif" in cif_data.filename:
        gemmi_data = np.array(
            [[item.replace(", ", ",") for item in row] for row in gemmi_data]
        )

    np.testing.assert_array_equal(parsnip_data, gemmi_data)
