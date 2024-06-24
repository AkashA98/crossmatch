from astroquery.simbad import Simbad
from astropy import units as u
import numpy as np

Simbad.TIMEOUT = 100000
Simbad.ROW_LIMIT = -1


class SimbadClass:
    """Simbad class to cross match sources with simbad"""

    def __init__(self, coords, radius=5) -> None:
        """Simbad class to handle all of Simbad queries

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): cross-match radius. Defaults
                to 5 arcsec.
            plot (bool, Optional): Plot WISE color plor? Defaults to False.
        """
        self.coords = coords
        self.rad = radius
        return None

    def _do_corssmatch(self):
        coords = self.coords
        rad = self.rad

        s = Simbad.query_region(coordinates=coords, radius=rad * u.arcsec)

        self.simbad = s
        has_match = np.zeros(len(coords)).astype(bool)
        match_name = np.array([""] * len(coords)).astype("U48")

        if not s is None:
            has_match[s["SCRIPT_NUMBER_ID"] - 1] = True
            match_name[s["SCRIPT_NUMBER_ID"] - 1] = s["MAIN_ID"]

        self.simbad_match = has_match
        self.simbad_match_name = match_name
