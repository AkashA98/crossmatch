from astropy.time import Time
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c

Vizier.TIMEOUT = 1000
Vizier.ROW_LIMIT = -1
import numpy as np

from .utils import correct_for_proper_motion
from .utils import catalogs


class GaiaClass:
    def __init__(self, coords, radius=5, epoch="2019-04-01T00:00:00") -> None:
        """Gaia class to handle all of GAIA queries

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): cross-match radius. Defaults
                to 5 arcsec.
            epoch (str, optional): observation epoch. Defaults
                to "2019-04-01T00:00:00".
        """
        self.coords = coords

        # Make a list of coordinates
        try:
            _ = len(coords)
        except TypeError:
            self.coords = coords.flatten()

        self.epoch = Time(epoch, scale="utc")
        self.rad = radius

    def _do_pm_significance(self, gaia):
        """GAIA query does not come with PM significance, so calculate by hand

        Args:
            gaia (astropy.table.table.Table): GAIA crossmatch result

        Returns:
            astropy.table.table.Table: input table with PM error added
        """
        pm = gaia["PM"].data.data * gaia["PM"].unit
        pmra = gaia["pmRA"].data.data * gaia["pmRA"].unit
        pmdec = gaia["pmDE"].data.data * gaia["pmDE"].unit
        pmra_err = gaia["e_pmRA"].data.data * gaia["e_pmRA"].unit
        pmdec_err = gaia["e_pmDE"].data.data * gaia["e_pmDE"].unit
        pmradec_corr = gaia["pmRApmDEcor"]
        pm_err = np.sqrt(
            ((pmra / pm) ** 2) * (pmra_err**2)
            + ((pmdec / pm) ** 2) * (pmdec_err**2)
            + (pmra * pmdec / pm**2) * (pmra_err * pmdec_err) * pmradec_corr
        )
        gaia.add_column(pm_err, name="e_PM", index=gaia.colnames.index("PM") + 1)
        return gaia

    def do_gaia_query(self):
        """Funtion to cross match with GAIA catalog"""
        coords = self.coords

        gaia = Vizier(columns=["all"]).query_region(
            coordinates=coords, radius=self.rad * u.arcsec, catalog=catalogs.GAIA
        )

        if len(gaia) > 0:
            # Then cross matches are found
            gaia = gaia[0]

            gaia_time = Time(["2016-01-01T00:00:00"] * len(gaia), scale="utc")
            gaia_coords = SkyCoord(ra=gaia["RA_ICRS"], dec=gaia["DE_ICRS"])

            # Correct the coordinates now
            pmra = gaia["pmRA"].data.data * gaia["pmRA"].unit
            pmdec = gaia["pmDE"].data.data * gaia["pmDE"].unit
            newcoords = correct_for_proper_motion(
                ra=gaia_coords.ra,
                dec=gaia_coords.dec,
                pmra=pmra,
                pmdec=pmdec,
                ref_time=gaia_time,
                current_time=self.epoch,
            )

            # Add proper motion error to the table as well
            gaia = self._do_pm_significance(gaia=gaia)

            # Add new coordinates
            gaia.add_column(newcoords.ra, name="RA_RACS_epoch")
            gaia.add_column(newcoords.dec, name="DEC_RACS_epoch")
            gaia.add_column(
                newcoords.separation(coords[gaia["_q"] - 1]).arcsec, name="distance"
            )

            gaia.sort(["_q", "distance"])
            self.gaia_res = gaia
            self.gaia_uniq = gaia[np.unique(gaia["_q"].data.data, return_index=True)[1]]

            # Now do basic classification
            self._filter_sources()
        else:
            self.gaia = None

    def _filter_sources(self):
        """Basic filter to classify GAIA sources"""
        coords = self.coords
        gaia_uniq = self.gaia_uniq
        # Store a mask which gives says which sources in our catalog has GAIA matches
        gaia_mask = np.zeros(len(coords)).astype(bool)
        gaia_mask[gaia_uniq["_q"] - 1] = True
        self.has_gaia_match = gaia_mask

        # Look for proper motion and parallax
        pm_mask = gaia_uniq["PM"].data.data / gaia_uniq["e_PM"] > 5
        plx_mask = gaia_uniq["RPlx"].data.data > 5

        star_mask = pm_mask | plx_mask
        is_star = np.zeros(len(coords)).astype(bool)
        is_star[gaia_uniq[star_mask]["_q"] - 1] = True
        self.is_gaia_star = is_star
        return None
