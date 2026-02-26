from astropy.timeseries import LombScargle
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

    # def __init__(self, stars=pd.read_parquet("Data/dan_dataset")):
    #     if not os.path.exists("Data/Candidates"):                           # <----- NEED TO CHANGE TO INPUT DAN'S STAR DATA
    #         os.makedirs("Data/Candidates")
    #     self.stars = stars

    def BTJDtoFrequency(self, time):
        """
        Converting the BTJD time input from the large A/F type
        star catalog to frequency in cycles/day which is used to
        plot the periodograms.
        Input:
            - Time (BTJD)
        Returns:
            - Frequencies (cycles/day)
        """

        # Calculate time span and sampling
        btjdTime = np.array(time)
        timeSpan = np.max(btjdTime) - np.min(btjdTime)
        N = len(btjdTime)

        # Calculating nyquist frequency from the average cadence in days
        avgCadence = timeSpan / (N - 1)
        nyquist = 1.0 / (2.0 * avgCadence)

        # Creating frequency array from 0
        resolution = 1.0 / timeSpan
        frequencies = np.arange(10, 200, resolution)

        print(f"Time span: {timeSpan:.2f} days")
        print(f"Average cadence: {avgCadence*24:.2f} hours")
        print(f"Nyquist frequency: ~{nyquist:.6f} cycles/day")
        print(f"Frequency resolution: {resolution:.6f} cycles/day")

        return frequencies

    def processStars(
        self, frequency, flux
    ):  # <----- Change to input data from Dan and Dan's giant catalog of A/F type stars
        """
        Using the time and flux values from large A/F type
        star catalog with short cadence data (200s exposures).
        This function goes through the stars and uses conditions
        to choose high, mid, and low-tier roAp candidates.
        Inputs:
            - Frequency (cycles/day)
            - Flux (electrons/s)
        Conditions:
            - High pulsation frequency (>=60 cycles/day).
            - Large, positive skew (High: >=5, Mid: <5 and >0, Low: ~0 or <0).
        Returns:
            - Normalized LC
            - Phase-Folded Diagram
            - High and Low Periodogram
            - Rotation Period
        """
        try:

            # To find the frequency of max peak, want to do for lower frequency spectrum to get rotation period
            fmax = np.max(frequency)
            rot_period = 2 / fmax

            pg = stitched.flatten().to_periodogram(normalization="amplitude")
            pg_low = stitched.to_periodogram(
                normalization="amplitude", maximum_frequency=5.0
            )  # For Finding the roataion period
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
            print(f"Analysis Failed: {e}")
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
            ra, dec = row.RA_ICRS, row.DE_ICRS
            i += 1
            print("-" * 40, f"\nIteration #{i}")

            lc, pg, tic_id, pgskewed, lcskewed, rot_period = self.processStars(ra, dec)

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
                        "ID": row.ID,
                        "TargetID": f"TIC {id}",
                        "GaiaDR3": row.GaiaDR3,
                        "Solar Mass": row.Mass50,
                        "RA (deg)": ra,
                        "DEC (deg)": dec,
                        "pmRA": row.pmRA,
                        "pmDEC": row.pmDE,
                        "Nyquist Frequency (1/d)": round(nyquist, 6),
                        "Max Power (ppm)": round(maxPw, 6),
                        "Frequency at max Power (1/d)": round(maxPwFreq, 4),
                        "Rotation Period (d)": round(rot_period, 4),
                        "GMAG": row.GMAG0,
                        "BP_RP": row.BP_RP0,
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


if __name__ == "__main__":
    GD = GetData()

    tic_id = "TIC 310817678"
    author = "TESS-SPOC"

    lc = lk.search_lightcurve(tic_id, exptime=200, author=author)
    stitch = lc.download_all().stitch()
    time = stitch.time.value
    flux = stitch.flux.value

    ls = LombScargle(time, flux)

    pg = stitch.flatten().to_periodogram(normalization="amplitude")

    f = GD.BTJDtoFrequency(time)
    p = ls.power(f)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 8))
    ax1, ax2 = axs.flatten()

    ax1.plot(f, p)
    ax2.plot(pg.frequency.value, pg.power.value)
    plt.xlabel("Frequency (cycles/day)")
    plt.ylabel("Power")
    ax1.set_title("Lomb-Scargle Periodogram")
    ax2.set_title("Lightkurve Periodogram")
    plt.show()
