from astroquery.ipac.ned import Ned
from astropy import units as u
import numpy as np

Ned.TIMEOUT = 100000
Ned.clear_cache()


class NEDClass:
    """NED class to cross match sources with simbad"""

    def __init__(self, coords, radius=5) -> None:
        """NED class to handle all of crossmatch queries

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

        s = Ned.query_region(coordinates=coords, radius=rad * u.arcsec)

        if len(s) > 0:
            has_match = np.zeros(len(coords)).astype(bool)
            has_match[s["No."] - 1] = True

            match_name = np.array([""] * len(coords)).astype("U48")
            match_name[s["No."] - 1] = s["Object Name"]

            self.ned = s
            self.ned_match = has_match
            self.ned_match_name = match_name
        else:
            self.ned = None
            self.ned_match = None
            self.ned_match_name = None
