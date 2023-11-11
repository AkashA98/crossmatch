from astropy.time import Time
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord

from gaia import GaiaClass
from pulsars import PulsarClass
from wise import WISEClass
from agns import AGNClass
from tns import TNSClass

Vizier.TIMEOUT = 300
Vizier.ROW_LIMIT = -1


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
        gaia = GaiaClass(coords=coords, radius=rad)
        gaia.do_gaia_query()
        self.gaia = gaia

        # Do ATNF pulsars
        pul = PulsarClass(coords=coords, radius=rad)
        pul.do_pulsar_query()
        self.pulsars = pul

        # Do AGNs
        agn = AGNClass(coords=coords, radius=rad)
        agn._do_crossmatch()
        self.agns = agn

        # Do WISE
        wise = WISEClass(coords=coords, radius=rad, plot=False)
        wise._do_query()
        self.wise = wise


if __name__ == "__main__":
    low = Table.read(
        "../racs_low_mid_catalogs/racs_low_racs_mid_not_in_tgss_initial_cut.parquet"
    )
    # low = Table.read("../tdes/tdes_in_racs.txt", format="ascii")
    # names = low["col0"]
    # names = [i[-3:] if "AT" in i else i[-4:] for i in names]
    coords = SkyCoord(ra=low["ra_deg_cont"], dec=low["dec_deg_cont"])

    cat = catalog(coords=coords)
    cat.crossmatch()
