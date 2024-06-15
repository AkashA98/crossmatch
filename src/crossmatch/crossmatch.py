import os, argparse
from astropy.time import Time
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord

Vizier.TIMEOUT = 1000
Vizier.ROW_LIMIT = -1

from .utils import *
from .gaia import GaiaClass
from .pulsars import PulsarClass
from .wise import WISEClass
from .agns import AGNClass
from .tns import TNSClass
from .simbad import SimbadClass
from .parse_blackcat import get_latest_vserion


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

        self.rad = radius

        # Corrdinates needs to be in the form of an array, so change that
        if len(coords.shape) == 0:
            coords = SkyCoord([coords.to_string("hmsdms")])
        self.coords = coords
        # We need to correct the positions. But this varies for different RACS
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

        # Do Simbad
        sim = SimbadClass(coords=coords, radius=rad)
        sim._do_corssmatch()
        self.simbad = sim

        # Do WISE
        wise = WISEClass(coords=coords, radius=rad, plot=False)
        wise._do_query()
        self.wise = wise

        # Cross match BLACKCAT
        # Make sure to parse blackcat
        get_latest_vserion()
        xrb_res = search_black_cat(coords)
        if xrb_res is not None:
            self.xrb_mask, self.xrbs = xrb_res
        else:
            self.xrbs = xrb_res


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""Script to to targetted searches for VAST data.
        Does force fitting for a given position and if needed searches around the 
        given coordinates to look for sources"""
    )
    parser.add_argument(
        "coords",
        help="""Source coordinates to search for, can be a csv file having the columns
         RA and DEC or a pair of strings with RA and DEC of the source (ascii format).
          For eg. (both of the following work, but beware of the '-' in declination)
         python crossmatch.py 12:54:24.3 +34:12:05.3 (can be hmsdms or degrees) or 
         python crossmatch.py my_list_of_coords.csv""",
        nargs="+",
        type=str,
    )

    args = parser.parse_args()

    coord = args.coords
    try:
        mask = os.path.isfile(args.coords[0])
        # This means a file with RA and DEC are given
        tab = Table.read(args.coords[0], format="ascii")
        coord = parse_coordinates(tab["RA"].data, tab["DEC"].data)
        try:
            names = tab["Name"]
            names = [i.replace(" ", "") for i in names]
        except NameError:
            names = coord.to_string("hmsdms")
    except (FileNotFoundError, TypeError):
        if len(coord) == 2:
            # This means a pair of coords is given
            coord = parse_coordinates(coord[0], coord[1])
        else:
            ra, dec = coord[0].split(" ")
            coord = parse_coordinates(ra, dec)
        names = coord.to_string("hmsdms")
    coords = SkyCoord("11h15m04.42s -61d09m39.80s")

    cat = catalog(coords=coords)
    cat.crossmatch()
