import copy
import numpy as np
from astropy.coordinates import SkyCoord


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
    nan_mask = (~np.isnan(pmra)) & (~np.isnan(pmdec)) & (~np.isnan(dt.value))

    oldcoords = SkyCoord(
        ra=ra, dec=dec, pm_ra_cosdec=pmra, pm_dec=pmdec, obstime=ref_time
    )

    # Analytical expressions to correct by hand
    # newra = ra + (pmra / np.cos(dec)) * (dt.value * u.day)
    # newdec = dec + (pmdec) * (dt.value * u.day)

    # Apply correction on when pm information is available

    newcoords = copy.deepcopy(oldcoords)

    corrected_coords = oldcoords[nan_mask].apply_space_motion(dt=dt[nan_mask])
    newcoords.ra[nan_mask] = corrected_coords.ra
    newcoords.dec[nan_mask] = corrected_coords.dec

    return newcoords
