"""""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """''
--------------------------------
Author: Jacob Romeo             |
Updated: 09/05/2025             |
Querying using Vizier Catalog   |
--------------------------------
       __.-._  ____ " Fear is the path to the dark side.
      '-._"7' /             Fear leads to anger.
       /'.-c                      Anger leads to hate.
       |  /T                             Hate leads to suffering." - Yoda, the wise coder
      _)_/LI

-> Outputs the queried catalog as a csv and parquet file of member stars in given clusters.
""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" """""" ""
import time
from astroquery.vizier import Vizier


# Vizier catalog query
def vizierQuery(
    cat="J/A+A/686/A42",
    input=-1,  # -1 is default for all rows
    col=[
        "Name",
        "ID",
        "GaiaDR3",
        "Prob",
        "RA_ICRS",
        "DE_ICRS",
        "pmRA",
        "pmDE",
        "Gmag",
        "BP-RP",
        "dist50",
        "logAge50",
        "Mass50",
    ],
):
    """
    Quering the Vizier catalog for clusters and members data
    specified catalog (J/A+A/686/A42), specifying the columns with 'cols'.
    """
    start_time = time.time()
    try:
        V = Vizier(columns=col)
        catalog = V.find_catalogs(cat)
        V.ROW_LIMIT = input  # or a large number
        print("Getting catalogs...")
        catalogs = V.get_catalogs(catalog.keys())
        print("Finished, converting and saving as csv and parquet files...")
        members = catalogs[1]
        members = members.to_pandas()
        members.to_csv("Data/members.csv", index=False)
        members.to_parquet("Data/members.parquet", index=False)
        print(members.columns)
    except Exception as e:
        print(f"Query failed: {e}")

    end = time.time()
    elapsed = (end - start_time) / 60
    print(f"Query completed in {elapsed:.2f} minutes, found {len(members)} stars...")


if "__main__" == __name__:
    vizierQuery()
