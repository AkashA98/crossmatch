import numpy as np
from astropy.time import Time
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c

Vizier.TIMEOUT = 1000
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
    NVSS = "VIII/65/nvss"
    FIRST = "VIII/92/first14"
    SUMSS = "VIII/81B/sumss212"


catalogs = Catalogtype


def dict_value_to_key(d, value):
    """Helper function to return the key, given a value"""
    keys = np.array(list(d.keys()))
    mask = [True if d[i] == value else False for i in keys]
    if np.any(mask):
        return keys[mask]
    else:
        raise KeyError


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
    nan_mask = (~np.isnan(pmra)) & (~np.isnan(pmdec)) & (~(dt.mask))

    # Astropy doesn't like nan's in time array when applying space motion
    # This is as of v 6.0, so break into two different coordinates

    # Analytical expressions to correct by hand
    # newra = ra + (pmra / np.cos(dec)) * (dt.value * u.day)
    # newdec = dec + (pmdec) * (dt.value * u.day)

    # Apply correction on when pm information is available
    if not np.any(nan_mask):
        return SkyCoord(ra=ra, dec=dec)
    else:
        new_ra, new_dec = ra, dec

        time = Time(np.array(ref_time.mjd[nan_mask]), format="mjd")
        oldcoords = SkyCoord(
            ra=ra[nan_mask],
            dec=dec[nan_mask],
            pm_ra_cosdec=pmra[nan_mask],
            pm_dec=pmdec[nan_mask],
            obstime=time,
        )

        corrected_coords = oldcoords.apply_space_motion(dt=np.array(dt.jd[nan_mask]))
        new_ra[nan_mask] = corrected_coords.ra
        new_dec[nan_mask] = corrected_coords.dec

    return SkyCoord(ra=new_ra, dec=new_dec)


def search_black_cat(tar_coords, text_file: str = "blackcat.txt"):
    """Helper function to crossmatch sources with XRBs using BLACKCAT

    Args:
        tar_coords (astropy.coordinates.sky_coordinate.SkyCoord): Target coordinates
        text_file (str, optional): path to the catalog. Defaults to "blackcat.txt".
    """
    blackcat = Table.read(text_file, format="ascii")
    coords = SkyCoord(
        ra=blackcat["RA"], dec=blackcat["Dec"], unit=(u.hourangle, u.degree)
    )
    ind, d2d, _ = tar_coords.match_to_catalog_sky(coords)
    mask = d2d.arcsec <= 5
    if not np.any(mask):
        return None
    else:
        return mask, blackcat[ind[mask]]
