import numpy as np
from astropy import units as u
from astropy.table import Table
from astroquery.vizier import Vizier
from matplotlib import pyplot as plt
from matplotlib import rc

rc("font", **{"family": "serif", "serif": ["Computer Modern Roman"]})
rc("text", usetex=True)

from utils import catalogs
from make_wise_cc_plot import wise_color_color_plot


class WISEClass:
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

    def _do_query(self):
        """Main function that does WISE query"""
        wise = Vizier(columns=["all"]).query_region(
            coordinates=self.coords, radius=self.rad * u.arcsec, catalog=catalogs.WISE
        )
        if len(wise) > 0:
            wise = wise[0]

            # Select unique ones
            wise.sort(["_q", "_r"])
            _, ind = np.unique(wise["_q"].data, return_index=True)
            wise = wise[ind]
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
            # w12_err = np.ones(len(w1)) * 0.5
            # w23_err = np.ones(len(w1)) * 0.5
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

            w12_err[(w12_low_lim + w12_upp_lim) > 0] = 0.25
            w23_err[(w23_low_lim + w23_upp_lim) > 0] = 0.25

            self.diff_colors = np.array([w12, w23])
            self.diff_errors = np.array([w12_err, w23_err])
            self.wise = wise

            self.lim_mask = Table(
                [w12_low_lim, w12_upp_lim, w23_low_lim, w23_upp_lim],
                names=[
                    "W12_low_limit",
                    "W12_upp_limit",
                    "W23_low_limit",
                    "W23_upp_limit",
                ],
                dtype=["bool", "bool", "bool", "bool"],
            )

            if self.plot:
                _, ax = wise_color_color_plot()
                for i in range(len(w12)):
                    ax.errorbar(
                        w23[i],
                        w12[i],
                        yerr=w12_err[i],
                        xerr=w23_err[i],
                        uplims=w12_upp_lim[i],
                        lolims=w12_low_lim[i],
                        xuplims=w23_upp_lim[i],
                        xlolims=w23_low_lim[i],
                        capsize=3,
                        barsabove=True,
                        fmt=".",
                        markersize=10,
                        color="k",
                    )
                ax.legend(loc="upper left", fontsize=15)
                plt.show()
        else:
            self.wise = None
