from astropy.table import Table, vstack
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

        # Ned queries do not work for multiple objects, so serialize it
        res = Table()
        ned_match = np.zeros(len(coords)).astype(bool)
        ned_match_name = np.array([None] * len(coords))
        for i, c in enumerate(coords):
            s = Ned.query_region(coordinates=c, radius=rad * u.arcsec)
            if len(s) > 0:
                s.sort("Separation")
                s.replace_column("No.", np.zeros(len(s)) + i)
                ned_match[i] = True
                ned_match_name[i] = s["Object Name"][0]
                res = vstack([res, s])
        self.ned = res
        self.ned_match = ned_match
        self.ned_match_name = ned_match_name
        return None
