import os, json, requests
from astropy import units as u
from astropy.table import Table


class TNSClass:
    """TNS class object to handle TNS queries"""

    tns_api_url = "https://www.wis-tns.org/api/get"

    def __init__(self, coords, radius=5) -> None:
        """Funtion to cross match given coordinates with TNS server.

        Args:
            coords (astropy.coordinates.sky_coordinate.SkyCoord):
                The coordinates that needs to be cross matched
            radius (float, optional): The cross match radius. Defaults
                to 5 arcsec.
        """
        self.coords = coords
        self.rad = radius

        # Get the bot id
        bot_id = os.environ["TNS_BOT_ID"]
        bot_name = os.environ["TNS_BOT_NAME"]
        bot_api_key = os.environ["TNS_API_KEY"]

        tns_bot_sign = (
            'tns_marker{"tns_id":'
            + bot_id
            + ', "type": "bot", "name": '
            + bot_name
            + "}"
        )
        self.tns_bot_sign = tns_bot_sign
        self.tns_api_key = bot_api_key

    def _make_search_objs(self):
        """Helper function to make search objects"""
        ra = [i.to_string(u.hour, sep=":") for i in self.coords.ra]
        dec = [i.to_string(u.degree, sep=":") for i in self.coords.dec]
        search_obj = [
            [
                ("ra", f"{ra[i]}"),
                ("dec", f"{dec[i]}"),
                ("radius", "5"),
                ("units", "arcsec"),
                ("objname", ""),
                ("objname_exact_match", 0),
                ("internal_name", ""),
                ("internal_name_exact_match", 0),
                ("objid", ""),
                ("public_timestamp", ""),
            ]
            for i in range(len(ra))
        ]
        self.search_obj = search_obj

    def _search(self):
        search_url = self.tns_api_url + "/search"
        tns_marker = self.tns_bot_sign
        headers = {"User-Agent": tns_marker}
        search_objs = self.search_obj

        # Store cross matches
        obj_name = []
        obj_id = []
        obj_type = []

        for obj in search_objs:
            search_data = {"api_key": self.tns_api_key, "data": json.dumps(obj)}
            response = requests.post(search_url, headers=headers, data=search_data)
            try:
                res = self._validate_search(response)
                res = res["data"]["reply"]
                if res != []:
                    obj_id.append(res["objid"])
                    obj_name.append(res["objname"])
                    obj_type.append(res["prefix"])
                else:
                    obj_id.append([""])
                    obj_name.append([""])
                    obj_type.append([""])
            except ConnectionAbortedError:
                obj_id.append([""])
                obj_name.append([""])
                obj_type.append([""])
        tab = Table([obj_id, obj_name, obj_type], names=["ID", "NAme", "Category"])
        self.res = tab

    def _validate_search(self, req):
        """Validate a given TNS request

        Args:
            req (requests.models.Response): TNS response
        """
        if req.status_code == 200:
            result = json.loads(req.text)
            return result
        else:
            raise ConnectionAbortedError
