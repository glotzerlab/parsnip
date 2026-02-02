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

        # Parse text inside a specific FONT tag. Ideally we'd take something more
        # clearly structured, but this does uniquely extract the data we want
        soup = BeautifulSoup(response.content, "html.parser")
        text = soup.find("font", {"face": "COURIER"}).get_text()

        # data formatted "x, y, z x...,y...,z..." -> [["x", "y", "z"], ...]
        data[str(i)] = [
            line.replace(" ", "") for line in text.splitlines() if "," in line
        ]
    return data


if __name__ == "__main__":
    OUTPUT_FN = "parsnip/symops.json"

    data = build_sg_dict()

    with open(OUTPUT_FN, "w") as f:
        json.dump(data, f, indent="  ")
