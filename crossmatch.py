import copy
import numpy as np
from psrqpy import QueryATNF
from astropy.time import Time
from astropy.table import Table
from astroquery.vizier import Vizier
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c


from make_wise_cc_plot import wise_color_color_plot

Vizier.TIMEOUT = 300
Vizier.ROW_LIMIT = -1

# Put in all the catalogs


AGN_catalogs = {
    "AllWISE_AGN": "J/ApJS/221/12/table1",
    "LQAC5": "J/A+A/624/A145/lqac5",
    "Milliquas": "VII/290/catalog",
    "WISE_R90": "J/ApJS/234/23/r90cat",
    "SDSS_DR16Q": "VII/289/dr16q",
}


class Catalogtype:
    GAIA = "I/355/gaiadr3"
    WISE = "II/328/allwise"
    AGN = AGN_catalogs


catalogs = Catalogtype


def correct_for_proper_motion(ra, dec, pmra, pmdec, ref_time, current_time):
    """Helper function to correct the catalogs for proper
        motion

    Args:
        ra (astropy.coordinates.sky_coordinate.SkyCoord):
            RA of the sources
        dec (astropy.coordinates.sky_coordinate.SkyCoord):
            Declination of the sources
        pmra (astropy.coordinates.sky_coordinate.SkyCoord):
            Proper motion in RA, includes the cos(dec) terms
        pmdec (astropy.coordinates.sky_coordinate.SkyCoord):
            Proper motion in declination.
        ref_time (astropy.time.core.Time): Time object for reference time
        current_time (astropy.time.core.Time): Time object for current epoch
    """

    dt = current_time - ref_time
    nan_mask = (~np.isnan(pmra)) & (~np.isnan(pmdec)) & (~np.isnan(dt.value))

    oldcoords = SkyCoord(
        ra=ra, dec=dec, pm_ra_cosdec=pmra, pm_dec=pmdec, obstime=ref_time
    )

    # Analytical expressions to correct by hand
    # newra = ra + (pmra / np.cos(dec)) * (dt.value * u.day)
    # newdec = dec + (pmdec) * (dt.value * u.day)

    newcoords = oldcoords.apply_space_motion(dt=dt)
    newcoords.ra[~nan_mask] = oldcoords.ra[~nan_mask]
    newcoords.dec[~nan_mask] = oldcoords.dec[~nan_mask]

    return newcoords


class Gaia:
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

            gaia_time = Time("2016-01-01T00:00:00", scale="utc")
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


class pulsars:
    def __init__(self, coords, radius=5, epoch="2019-04-01T00:00:00") -> None:
        """Pulsar class to handle all of pulsar queries

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): cross-match radius. Defaults
                to 5 arcsec.
            epoch (str, optional): observation epoch. Defaults
                to "2019-04-01T00:00:00".
        """
        self.coords = coords
        self.epoch = Time(epoch, scale="utc")
        self.rad = radius

    def do_pulsar_query(self):
        """Funtion to cross match with known pulsars"""
        coords = self.coords

        # Load pulsar catalog
        psrs = QueryATNF(params=["JNAME", "RAJD", "DECJD", "PMRA", "PMDEC", "POSEPOCH"])
        psrs = psrs.table
        psr_coords = SkyCoord(ra=psrs["RAJD"], dec=psrs["DECJD"])

        posepoch = Time(psrs["POSEPOCH"], format="mjd")
        pmra = psrs["PMRA"].data.data * psrs["PMRA"].unit
        pmdec = psrs["PMDEC"].data.data * psrs["PMDEC"].unit
        print(psr_coords, pmra, pmdec, posepoch)
        newpos = correct_for_proper_motion(
            ra=psr_coords.ra,
            dec=psr_coords.dec,
            pmra=pmra,
            pmdec=pmdec,
            ref_time=posepoch,
            current_time=self.epoch,
        )

        # Cross match it with sources
        ind, d2d, _ = coords.match_to_catalog_sky(newpos)
        match_mask = d2d.arcsec <= 5
        pulsar_match = np.zeros(len(coords)).astype(bool)
        if np.any(match_mask):
            pulsar_match[match_mask] = True
            self.has_pulsar_match = pulsar_match

            matched_pulsars = np.array([None] * len(coords))
            matched_pulsars[match_mask] = psrs["JNAME"][ind[match_mask]]
            self.match_pulsars = matched_pulsars

            match_pulsars_info = psrs[ind[match_mask]]
            match_pulsars_info.add_column(
                np.arange(len(coords))[match_mask], name="_q", index=0
            )
            self.match_pulsars_info = match_pulsars_info
        else:
            self.has_pulsar_match = pulsar_match


class WISE:
    """Wise class to cross match with WISE catalog"""

    def __init__(self, coords, radius=5, plot=False) -> None:
        """WISE class to handle all of WISE queries

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): cross-match radius. Defaults
                to 5 arcsec.
            plot (bool, Optional): Plot WISE color plor? Defaults to False.
        """
        self.coords = coords
        self.rad = radius
        self.plot = plot
        return None

    def do_query(self):
        """Main function that does WISE query"""
        wise = Vizier(columns=["all"]).query_region(
            coordinates=self.coords, radius=self.rad * u.arcsec, catalog=catalogs.WISE
        )
        if len(wise) > 0:
            wise = wise[0]
            # Try to get a classifer based on colors
            w1, w2, w3 = wise["W1mag"], wise["W2mag"], wise["W3mag"]
            e_w1, e_w2, e_w3 = wise["e_W1mag"], wise["e_W2mag"], wise["e_W3mag"]

            # Check the quality flags
            w1_mask = np.array(
                [True if i[0] in ["A", "B"] else False for i in wise["qph"]]
            )
            w2_mask = np.array(
                [True if i[1] in ["A", "B"] else False for i in wise["qph"]]
            )
            w3_mask = np.array(
                [True if i[2] in ["A", "B"] else False for i in wise["qph"]]
            )

            w12 = w1 - w2
            w23 = w2 - w3

            self.colors = np.array([w1, w2, w3])
            self.errors = np.array([e_w1, e_w2, e_w3])
            self.flags = np.array([w1_mask, w2_mask, w3_mask])

            # Calculate errors
            w12_err = np.ones(len(w1)) * np.nan
            w23_err = np.ones(len(w1)) * np.nan
            w12_err = np.sqrt(e_w1**2 + e_w2**2)
            w23_err = np.sqrt(e_w2**2 + e_w3**2)
            w12_upp_lim = np.zeros(len(w1))
            w12_low_lim = np.zeros(len(w1))
            w12_upp_lim[~w1_mask] = 1
            w12_low_lim[~w2_mask] = 1

            w23_upp_lim = np.zeros(len(w2))
            w23_low_lim = np.zeros(len(w2))
            w23_upp_lim[~w2_mask] = 1
            w23_low_lim[~w3_mask] = 1

            self.diff_colors = np.array([w12, w23])
            self.diff_errors = np.array([w12_err, w23_err])
            self.wise = wise

            if self.plot:
                _, ax = wise_color_color_plot()
                ax.errorbar(
                    w23,
                    w12,
                    yerr=w12_err,
                    xerr=w23_err,
                    uplims=w12_upp_lim,
                    lolims=w12_low_lim,
                    xuplims=w23_upp_lim,
                    xlolims=w23_low_lim,
                    capsize=2,
                    barsabove=True,
                    fmt=".",
                    color="k",
                )
                plt.show()
        else:
            self.wise = None


class catalog:
    def __init__(self, coords, radius=5) -> None:
        """Funtion to cross match given coordinates with variuos
            O/IR surveys

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): The cross match radius. Defaults
                to 5 arcsec.
        """
        self.coords = coords
        self.rad = radius

        # We need to correct the positions. But this varies fpr different RACS
        # sources, so define a common time
        self.epoch = Time("2019-05-01T00:00:00", format="isot", scale="utc")

    def crossmatch(self):
        """Main function to crosmacth with various catalogs"""
        coords = self.coords
        rad = self.rad

        # Do GAIA
        gaia = Gaia(coords=coords, radius=rad)
        gaia.do_gaia_query()
        self.gaia = gaia

        # Do ATNF pulsars
        pul = pulsars(coords=coords, radius=rad)
        pul.do_pulsar_query()
        self.pulsars = pul

    def agn_filter(self):
        """Funtion to cross match with AGN catalogs"""
        coords = self.coords
        agn_cat_names = list(catalogs.AGN.keys())
        agn_cats = [catalogs.AGN[i] for i in agn_cat_names]

        has_agn_match = np.zeros(len(coords)).astype(bool)
        agn_match_info = {}
        for i in range(len(agn_cats)):
            agns = Vizier(columns=["all"]).query_region(
                coordinates=coords,
                radius=5 * u.arcsec,
                catalog=agn_cats[i],
            )
            if len(agns) > 0:
                has_agn_match[agns[0]["_q"] - 1] = True
                agn_match_info[agn_cat_names[i]] = agns[0]
            else:
                # which means none of the catalogs have any matches
                agn_match_info[agn_cat_names[i]] = None
        self.has_agn_match = has_agn_match
        self.agn_match_info = agn_match_info


if __name__ == "__main__":
    low = Table.read("../racs_low_racs_mid_initial_cut.parquet")

    coords = SkyCoord(ra=low["ra_deg_cont"], dec=low["dec_deg_cont"])
    # coords = SkyCoord(ra=["17h25m23.45s"], dec=["-30d37m20.61s"])
    cat = catalog(coords=coords[:200])
    # cat = catalog(coords=coords)
    # cat.gaia_filter()
    # cat.pulsar_filter()
    # cat.agn_filter()
    # cat.crossmatch()
    w = WISE(coords=coords[:100], plot=True)
    w.do_query()
