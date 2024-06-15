import os, requests
import numpy as np
from bs4 import BeautifulSoup
from astropy.table import Table


def get_latest_vserion():
    """Update BLACKCAT"""
    request = requests.get("https://www.astro.puc.cl/BlackCAT/transients.php")
    body = BeautifulSoup(request.text, "html.parser")

    # Get all rows of the table
    rows = body.find_all("tbody")

    # Get the columns
    # columns = rows[0].find_all("td")
    # colnames = [i.attrs["class"][0] for i in columns]
    tab = Table(
        names=[
            "Id",
            "Name",
            "RA",
            "Dec",
            "GL",
            "GB",
            "Mag burst",
            "Mag qui",
            "dis",
            "flux_xray",
            "period",
        ],
        dtype=["U50"] * 11,
    )

    for r in rows:
        row_content = r.find_all("td")
        id = int(row_content[0].find_all("a")[0].text)
        name = row_content[1].find_all("a")[0].text
        name = name.replace(" ", "").replace("\n", ",").replace("=", ",")
        # print([i.contents for i in row_content[2:]])
        values = [i.text.replace("\n", "").replace(" ", "") for i in row_content[2:]]

        row = np.concatenate(([id, name], values))
        tab.add_row(row)

    tab["Id"] = tab["Id"].astype(int)

    if not os.path.isdir("data/"):
        os.mkdir("data/")
    tab.write("data/blackcat.txt", format="ascii", delimiter="\t", overwrite=True)
