import numpy as np
from astropy import units as u
from astroquery.vizier import Vizier
from astropy import units as u, constants as c

Vizier.TIMEOUT = 1000
Vizier.ROW_LIMIT = -1
from .utils import catalogs, dict_value_to_key


class AGNClass:
    """AGN class to cross match with known AGN catalogs"""

    def __init__(self, coords, radius=5) -> None:
        """Funtion to cross match given coordinates with variuos
            AGN catalogs

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): The cross match radius. Defaults
                to 5 arcsec.
        """
        self.coords = coords
        self.rad = radius

    def _do_crossmatch(self):
        """Main function to do the cross match"""
        cat = catalogs.AGN
        coords = self.coords
        rad = self.rad

        # Get the Vizier DOI entries
        dois = list(cat.values())

        # Do the search
        agn = Vizier(columns=["all"]).query_region(
            coordinates=coords, radius=rad * u.arcsec, catalog=dois
        )

        # Parse the output
        match_dict = {}
        match_masks = []

        res_cat = list(agn.keys())
        if len(res_cat) > 0:
            for i, n in enumerate(res_cat):
                try:
                    name = dict_value_to_key(cat, n)[0]
                    ind_cat = agn[i]
                    match_dict[name] = ind_cat
                    ind_mask = np.zeros(len(coords), dtype="bool")
                    ind_mask[ind_cat["_q"] // 5] = True
                    match_masks.append(ind_mask)
                except KeyError:
                    continue

        match_masks = np.array(match_masks)
        self.match_res = match_dict
        self.match_ind_mask = match_masks
        self.match_mask = match_masks.sum(0).astype(bool)
