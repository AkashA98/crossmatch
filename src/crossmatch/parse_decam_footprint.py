from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Ellipse
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.table import Table
from matplotlib import pyplot as plt


def get_matched_files(path, s, dr2=True):
    """Read the DeCaps image data base file and parse
        the data.

    Args:
        path (str): Path to the imdb.fits file
        s (astropy.coordinates.skycoord.SkyCoord): Coordinates
    """

    hdu = fits.open(path)

    ext1 = hdu[1].data

    c1 = np.hstack((ext1["CORN1RA"], ext1["CORN1DEC"]))
    c2 = np.hstack((ext1["CORN2RA"], ext1["CORN2DEC"]))
    c3 = np.hstack((ext1["CORN3RA"], ext1["CORN3DEC"]))
    c4 = np.hstack((ext1["CORN4RA"], ext1["CORN4DEC"]))

    corners = np.array([c1, c2, c3, c4])

    mask = get_files(corners=corners, s=s)
    match_files = ext1["DTNSANAM"][mask]

    if dr2:
        dr2_mask = [True if "decaps2" in i else False for i in match_files]
        match_files = match_files[dr2_mask]

    cat_files = [i.replace(".fits.fz", ".cat.fits") for i in match_files]
    return match_files, np.array(cat_files)


def get_files(corners, s, size=0.1):
    """Get the corresponding files for which the gives corrdinates
        are in the image

    Args:
        corners (np.array.ndarray): Corners of the image
        s (astropy.coordinates.skycoord.SkyCoord): Coordinates
        size (float, optional): Image size (in degree)
    """

    req = np.array([s.ra.degree, s.dec.degree])[np.newaxis, :]
    diff = corners - req
    isdec_ok = (diff[0][:, 1] + size) * (diff[2][:, 1] - size) < 0
    isra_ok = (diff[0][:, 0] + size) * (diff[1][:, 1] + size) < 0
    isimage_ok = isra_ok & isdec_ok

    return isimage_ok


def get_sources(path, s, rad=15, nneighbor=7.5):
    """Read the DeCaps image data base file and select
        sources within a given cross match radius

    Args:
        path (str): Path to the band merged catalog
        s (astropy.coordinates.skycoord.SkyCoord): Coordinates
        rad (float, optional): cross match radius
    """
    hdu = fits.open(path)

    ext1 = hdu[1].data
    coords = SkyCoord(ra=ext1["ra"] * u.degree, dec=ext1["dec"] * u.degree)
    sep = coords.separation(s)
    sep_mask = sep.arcmin <= rad

    # Select isolated sources (no co crowding, separation alteast 7.5'')
    _, d2d, _ = coords.match_to_catalog_sky(coords, nthneighbor=2)
    iso_mask = d2d.arcsec >= nneighbor

    # Demand atleast 5 band detection
    ndet_mask = ((ext1["nmag"] >= 5).sum(1)) >= 5

    # Demand a minimum SNR of 20
    snr_mask = (ext1["snr_avg"] >= 20).sum(1) >= 5

    # Select against heavily moving
    pos_mask = ext1["posstdev"] <= 0.25

    mask = iso_mask & sep_mask & ndet_mask & snr_mask & pos_mask

    return ext1[mask]


def plot_decaps(img):
    corr = Table.read("decam_vast_match.txt", format="ascii")
    ra_err = corr["ra"].data - corr["ra_deg_cont"].data
    dec_err = corr["dec"].data - corr["dec_deg_cont"].data

    s = SkyCoord("11h42m53.59s -64d57m08.80s")
    s_corr = SkyCoord(
        ra=(s.ra + np.median(ra_err) * u.degree),
        dec=(s.dec + np.median(dec_err) * u.degree),
    )
    hdu = fits.open(img)
    w = WCS(hdu[0].header)
    d = hdu[0].data
    scl = 0.2637
    pos = np.array(w.world_to_pixel(s))

    cut = Cutout2D(d, pos, size=15 * u.arcsec, wcs=w)

    e = Ellipse(
        xy=np.array(cut.wcs.world_to_pixel(s)),
        width=2 * 0.2120 / scl,
        height=2 * 0.2424 / scl,
        angle=13.06,
        fill=False,
        color="r",
    )
    e2 = Ellipse(
        xy=np.array(cut.wcs.world_to_pixel(s)),
        width=2 * 2 * 0.2120 / scl,
        height=2 * 2 * 0.2424 / scl,
        angle=13.06,
        fill=False,
        color="r",
    )
    e3 = Ellipse(
        xy=np.array(cut.wcs.world_to_pixel(s)),
        width=3 * 0.2120 / scl,
        height=3 * 0.2424 / scl,
        angle=13.06,
        fill=False,
        color="r",
    )
    e4 = Ellipse(
        xy=np.array(cut.wcs.world_to_pixel(s)),
        width=2 * 2.5 / scl,
        height=2 * 2.5 / scl,
        angle=13.06,
        fill=False,
        color="r",
    )

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection=cut.wcs)
    ax.imshow(cut.data, cmap="gray_r")

    cut_pos = np.array(cut.wcs.world_to_pixel(s))
    ax.scatter(cut_pos[0], cut_pos[1], marker="x", s=30, c="r", lw=2)

    cut_pos_corr = np.array(cut.wcs.world_to_pixel(s_corr))
    ax.scatter(cut_pos_corr[0], cut_pos_corr[1], marker="x", s=30, c="b", lw=2)
    ax.add_patch(e)
    ax.add_patch(e2)
    ax.add_patch(e3)
    ax.add_patch(e4)
    plt.tight_layout(rect=(0.15, 0.1, 1, 1))
    plt.show()
