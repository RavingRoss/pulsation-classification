from astroquery.mast import Catalogs
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from scipy.stats import skew
import astropy.units as u
import lightkurve as lk
import pandas as pd
import numpy as np
import time
import os


class GetData:
    """
    Class to get the TESS lightcurve data for the candidate roAp stars
    """

    def __init__(self, stars=pd.read_parquet("Data/members.parquet")):
        if not os.path.exists("Data/Candidates"):
            os.makedirs("Data/Candidates")
        self.stars = stars

    def getShortestCadence(self, ra, dec):
        """
        Using the RA and DEC, found from the Gaia query in 'GetData' class,
        of the star to get its lightcurve data from the shortest
        cadence available for that star, using the lightkurve library.
        """
        coords = coord.SkyCoord(ra, dec, unit="deg")
        ticResult = Catalogs.query_region(coords, radius=5 * u.arcsec, catalog="TIC")  # type: ignore
        tic_id = ticResult["ID"][0]

        if len(ticResult) == 0:
            print(f"No TIC match found near RA={ra}, DEC={dec}")
            return None, None, None, None, None, None

        search = lk.search_lightcurve(
            f"TIC {tic_id}", author="TESS-SPOC", exptime=200
        )  # SPOC

        if len(search) == 0:
            print(
                f"No light curves found for TIC {tic_id}, (RA={ra}, DEC={dec})"
            )  # Debugging
            return None, None, None, None, None, None

        print("Printing original:\n", search)

        # Step 1: Find the shortest cadence
        shortest_cadence = min(
            [entry.exptime.value for entry in search if entry.exptime is not None]
        )
        print(f"Shortest cadence found: {shortest_cadence} s")

        # Step 2: Filter entries with that exact cadence and valid target name
        valid_entries = []
        for entry in search:
            try:
                if (
                    entry.exptime is not None
                    and entry.exptime.value == shortest_cadence
                ):
                    lc = entry.download()
                    if not np.all(np.isnan(lc.flux.value)):  # type: ignore
                        valid_entries.append(entry)
                        print(f"Valid entry added: {entry.exptime.value}")
                    else:
                        print(f"Skipped NaN flux entry: {entry.exptime.value}")
            except Exception as e:
                print(f"Error checking entry: {e}")
                continue

        if not valid_entries:
            print(f"No valid light curves found with cadence {shortest_cadence} s")
            return None, None, None, None, None, None

        # Step 3: Download and stitch all valid entries
        try:
            downloaded = lk.LightCurveCollection(
                [entry.download() for entry in valid_entries]
            )
            print(f"Number of light curves downloaded: {len(downloaded)}")
            if len(downloaded) == 1:
                stitched = downloaded[0]
                norm_lc = stitched.normalize()  # type: ignore
                print("No stitching needed...")
            else:
                stitched = downloaded.stitch()
                norm_lc = stitched.normalize()
                print("Stitched light curve...")

            pg = stitched.flatten().to_periodogram(normalization="amplitude")  # type: ignore
            # For Finding the roataion period
            pg_low = stitched.to_periodogram(normalization="amplitude", maximum_frequency=5.0)  # type: ignore
            print(f"Converted to periodogram...")
            pg.show_properties()

            peaks, _ = find_peaks(pg.power, height=1e-4)

            if peaks.shape == (0,):  # Filtering out the stars with no features
                print("No features to analyze...")
                return None, None, None, None, None, None

            pgskewed = skew(pg.power, nan_policy="omit")
            print("PG Skew value:", pgskewed)

            flux = np.array(norm_lc.flux.to_value())
            lcskewed = skew(flux, nan_policy="omit")
            print("LC Skew value:", lcskewed)

            # To find the frequency of max peak, want to do for lower frequency spectrum to get rotation period
            fmax = pg_low.frequency_at_max_power.to_value()
            rot_period = 2 / fmax  # Might need to change this depending on the star

            return norm_lc, pg, tic_id, pgskewed, lcskewed, rot_period

        except Exception as e:
            print(f"Download failed for TIC {tic_id}: {e}")
            return None, None, None, None, None, None

    def getTESSData(self):
        """
        Using the RA and DEC, found from the Gaia query in 'GetData' class,
        of the star to get its lightcurve data using the lightkurve library.
        """
        cands = self.stars
        roap_cands = []

        print("Starting TESS data download...")
        i = 0
        for row in cands.itertuples(index=False):
            ra, dec = row.RA_ICRS, row.DE_ICRS  # type: ignore
            i += 1
            print("-" * 40, f"\nIteration #{i}")

            lc, pg, tic_id, pgskewed, lcskewed, rot_period = self.getShortestCadence(
                ra, dec
            )

            if lc is None:
                continue

            print(f"Downloaded data for TIC {tic_id}")

            id = pg.targetid
            nyquist = pg.nyquist.value
            maxPwFreq = pg.frequency_at_max_power.value
            maxPw = pg.max_power.value

            # Adding data to the stars list
            if (ra in cands["RA_ICRS"].values) and (dec in cands["DE_ICRS"].values):
                roap_cands.append(
                    {
                        "ID": row.ID,  # type: ignore
                        "TargetID": f"TIC {id}",
                        "GaiaDR3": row.GaiaDR3,  # type: ignore
                        "Solar Mass": row.Mass50,  # type: ignore
                        "RA (deg)": ra,
                        "DEC (deg)": dec,
                        "pmRA": row.pmRA,  # type: ignore
                        "pmDEC": row.pmDE,  # type: ignore
                        "Nyquist Frequency (1/d)": round(nyquist, 6),
                        "Max Power (ppm)": round(maxPw, 6),
                        "Frequency at max Power (1/d)": round(maxPwFreq, 4),
                        "Rotation Period (d)": round(rot_period, 4),
                        "GMAG": row.GMAG0,  # type: ignore
                        "BP_RP": row.BP_RP0,  # type: ignore
                        "PG Skew": pgskewed,
                        "LC Skew": lcskewed,
                    }
                )
                print(f"Added TIC {id} to candidates list.")

    def plotAll(self, lc, pg, tic_id):
        """
        Plotting the lightcurve and periodogram of the star
        """
        if lc is None or pg is None:
            print(f"No data to plot for TIC {tic_id}")
            return

        plt.figure(figsize=(14, 6))

        # Plotting the light curve
        plt.subplot(1, 2, 1)
        lc.plot()
        plt.title(f"TIC {tic_id} Light Curve")

        # Plotting the periodogram
        plt.subplot(1, 2, 2)
        pg.plot()
        plt.title(f"TIC {tic_id} Periodogram")

        plt.tight_layout()
        plt.show()
