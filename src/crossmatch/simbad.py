from astroquery.simbad import Simbad
from astropy import units as u
import numpy as np


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

        if not s is None:
            has_match = np.zeros(len(coords)).astype(bool)
            has_match[s["SCRIPT_NUMBER_ID"] - 1] = True

            match_name = np.array([""] * len(coords)).astype("U48")
            match_name[s["SCRIPT_NUMBER_ID"] - 1] = s["MAIN_ID"]

            self.simbad = s
            self.simbad_match = has_match
            self.simbad_match_name = match_name
        else:
            self.simbad = None
            self.simbad_match = None
            self.simbad_match_name = None
