"""Build the nodes.geojson file."""

from importlib import resources

import geopandas as gpd
import requests


def build_geojson_base_layer(
    resolution: int = 10, year: int = 2021, crs: int = 4326
) -> None:
    """Download and reshape the NUTS regions to meet ESM chart needs.

    Parameters
    ----------
    resolution
        The resolution of the Polygons geometries.
    year
        The year of geometry publication.
    crs
        The coordinate reference system.

    Returns
    -------
    :
        Writes the reshaped geojson to the data directory.

    Raises
    ------
    HttpError
        For failed GET requests to the server.
    """
    scheme = "https"
    server = "gisco-services.ec.europa.eu"
    path = "distribution/v2/nuts/geojson"
    file = f"NUTS_RG_{resolution}M_{year}_{crs}.geojson"
    url = f"{scheme}://{server}/{path}/{file}"

    response = requests.get(url)
    response.raise_for_status()
    gdf = gpd.GeoDataFrame.from_features(response.json(), crs=crs)

    def nuts_id_to_esm_node(nuts_id: str) -> str:
        """Transform NUTS3 region codes.

        Parameters
        ----------
        nuts_id
            The NUTS3 code.

        Returns
        -------
        :
            The NUTS2 regions for Austria, the NUTS1 regions for Germany
            and the country code for all other countries.
        """
        if nuts_id.startswith("AT"):
            return nuts_id[:4]
        elif nuts_id.startswith("DE"):
            return nuts_id[:3]
        return nuts_id[:2]

    gdf["ESM_ID"] = gdf["NUTS_ID"].apply(nuts_id_to_esm_node)

    nuts_to_esm_mapping = {
        "DE1": "DE0 1",  # Baden-Würthemberg, DE0 0
        "DE2": "DE0 0",  # Bavaria, DE0 1
        "DE3": "DE0 3",  # Berlin -> DE - Mideast
        "DE4": "DE0 3",  # Brandenburg -> DE - Mideast
        "DE5": "DE0 4",  # Bremen -> DE - North
        "DE6": "DE0 4",  # Hamburg -> DE - North
        "DE7": "DE0 2",  # Hessen -> DE - Midwest
        "DE8": "DE0 3",  # Mecklenburg-Vorpommern -> DE - Mideast
        "DE9": "DE0 4",  # Niedersachsen -> DE - North
        "DEA": "DE0 2",  # Nordrhein-Westfalen -> DE - Midwest
        "DEB": "DE0 2",  # Rheinland-Pfalz -> DE - Midwest
        "DEC": "DE0 2",  # Saarland -> DE - Midwest
        "DED": "DE0 3",  # Sachsen -> DE - Mideast
        "DEE": "DE0 3",  # Sachsen-Anhalt -> DE - Mideast
        "DEF": "DE0 4",  # Schleswig-Holstein -> DE - North
        "DEG": "DE0 3",  # Thüringen -> DE - Mideast
        "AT11": "AT0 0",  # Burgenland
        "AT12": "AT0 1",  # Niederösterreich
        "AT13": "AT0 2",  # Wien
        "AT21": "AT0 3",  # Kärnten
        "AT22": "AT0 4",  # Steiermark
        "AT31": "AT0 5",  # Oberösterreich
        "AT32": "AT0 6",  # Salzburg
        "AT33": "AT0 8",  # Tirol
        "AT34": "AT0 7",  # Vorarlberg
    }
    # gdf = gdf.replace({"ESM_ID": nuts_to_esm_mapping})

    nodes = gdf.dissolve(by="ESM_ID")
    duplicates = [s for s in nodes.index if s.startswith(("AT", "DE")) and len(s) != 5]
    nodes = nodes.drop(duplicates)

    # some countries are not modeled
    nodes = nodes.drop(["TR", "IS", "CY"])

    file_path = resources.files("evals") / "data" / "nodes.geojson"
    nodes["geometry"].to_file(file_path, driver="GeoJSON")


if __name__ == "__main__":
    build_geojson_base_layer()
