import numpy as np
from psrqpy import QueryATNF
from astropy.time import Time
from astropy.coordinates import SkyCoord
from utils import correct_for_proper_motion


class PulsarClass:
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
