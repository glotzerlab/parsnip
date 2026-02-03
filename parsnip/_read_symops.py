# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.


if __name__ == "__main__":
    import json

    import requests
    import tqdm
    from bs4 import BeautifulSoup

    def build_sg_dict():
        """Extract symops from the hypertext crystallography book."""
        data = {}
        for i in tqdm.tqdm(range(1, 231)):
            url = f"http://img.chem.ucl.ac.uk/sgp/LARGE/{i:03d}az3.htm"
            response = requests.get(url, timeout=10)
            soup = BeautifulSoup(response.content, "html.parser")

            # Extract symops
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

            # { "1": { "xyz": [...], "HM": "P 1" } }
            data[str(i)] = {"xyz": symops, "HM": hm_symbol}

        return data

    OUTPUT_FN = "parsnip/symops.json"

    data = build_sg_dict()

    with open(OUTPUT_FN, "w") as f:
        json.dump(data, f, indent="  ")
