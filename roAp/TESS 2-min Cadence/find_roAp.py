"""
------------------------------------
Author: Jacob Romeo                 |
Updated: 02/16/2026                 |
TESS Cycle 2 Obs with 2-min Cadence |
------------------------------------
       __.-._  ____ " Fear is the path to the dark side.
      '-._"7' /             Fear leads to anger.
       /'.-c                      Anger leads to hate.
       |  /T                             Hate leads to suffering." - Yoda, the wise coder
      _)_/LI

-> UPDATES <-
✔️ Copied most of Test.py since the plotting tools are very good for that. But need
edit so that input data is 'time' and 'flux' (Data from sarek1).
✔️ Using the automation tool from Test.py to find peaks in each periodogram.
    + NOT stitching the LCs as they look weird, so only using one LC download instead of all (if more than one were found).
    + Need to come up with a way to find the best possible LC-set if there are > 1 LCs when searching.
    + Or maybe since we are given flux and time, when using "LightCurve()" it only returns a single LC-set.
- Everything is working, now need to add classification system for low, mid, and high end candidates.
"""

import os

import lightkurve as lk
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import skew

# --------------------------------------------------------------------------------------------------------------
# Where the fun begins:


class find_roAp:
    """
    Class to find roAp pulsations in TESS data using given data from sarek1 server which
    outputs the Flux and Time of stars. With this we search through the database to find
    roAp candidates.
    - Functions:
        - __init__(self, flux, time): Initializes the class with the flux and time data to
                                      get the LightCurve.
        - cpd_to_mhz(self): Converts cycles/day to mHz.
        - mhz_to_cpd(self): Converts mHz to cycles/day.
        - find_peaks(self): Finds peaks in the power spectrum.
        - calculate_skewness(self): Calculates the skewness of the power spectrum.
        - plot_lightcurve(self): Plots the high and low power spectra with the Normalized
                                 and phase-folded light curves.
    """

    def __init__(self, flux, time):
        """
        Initializes the class with the flux and time data to get the LightCurve.

        Also gives the rotation period, conversions to mmags, and the conversion
        factor for cycles/day to mHz.
        """
        self.flux = flux
        self.time = time
        self.light_curve = lk.LightCurve(time=time, flux=flux)  # type: ignore

        try:
            if len(self.light_curve) == 0:
                raise ValueError("Light curve is empty")

            self.pg = self.light_curve.flatten().to_periodogram(
                normalization="amplitude"
            )
            self.pg.show_properties()

            self.pg_low = self.light_curve.to_periodogram(
                normalization="amplitude", maximum_frequency=5.0
            )

            # To find the frequency of max peak, want to do for lower frequency spectrum to get rotation period
            fmax = self.pg_low.frequency_at_max_power.to_value()
            self.rot_period = 2 / fmax

            # Threshold to find peaks
            self.mmag_lps = (
                0.0010857362047581294 * 1e6 * self.pg_low.power.value
            )  # Convert to mmag for lower frequency Periodogram
            self.mmag_ps = (
                0.0010857362047581294 * 1e6 * self.pg.power.value
            )  # Convert to mmag for main Periodogram

            # Define conversion functions
            self.factor_mHz = 1 / 86400 * 1000  # cycles/day → mHz
        except ValueError as e:
            print(e)

    def cpd_to_mhz(self, x):
        """
        Converts cycles/day to mHz.
        - Returns:
            - mHz (float): Frequency in mHz.
        """
        return x * self.factor_mHz

    def mhz_to_cpd(self, x):
        """
        Converts mHz to cycles/day.
        - Returns:
            - cpd (float): Frequency in cycles/day.
        """
        return x / self.factor_mHz

    def find_peaks(self):
        """
        Finds peaks in the power spectrum.
        - Returns:
            - num (int): Arbituary number to find threshold.
            - threshold (int): Threshold value for peak detection.
            - freq (list): List of frequencies of the peaks.
            - peaks (list): List of peak heights.
            - spacing (list): List of spacings between consecutive peaks.
        """
        freq = []
        mean_power = np.nanmean(self.mmag_ps)
        max_power = np.nanmax(self.mmag_ps)
        num = 5
        threshold = num * mean_power
        print(
            f"Max Power: {max_power} mmag | Mean Power: {mean_power} mmag | Threshold set to: {threshold} mmag"
        )

        indices, properties = find_peaks(self.mmag_ps, height=float(threshold))
        peaks = properties["peak_heights"]

        # Get the frequencies of the peaks
        peak_freqs = self.pg.frequency[indices].to_value()

        # Calculate the spacing between consecutive peaks
        spacing = np.diff(peak_freqs)

        for p in indices:
            f = self.pg.frequency[p].to_value()
            freq.append(f)

        return num, threshold, freq, peaks, spacing

    def calculate_skewness(self):
        """
        Calculates the skewness of the power spectrum.
        - Returns:
            - skew (float): Skewness of the power spectrum.
        """
        return skew(self.mmag_ps, nan_policy="omit")

    def classify_candidates(self):
        """
        Classifies the candidates based on their skewness.
        - Returns:
            - classification (str): Classification of the candidates.
        """

        if not os.path.exists("Data/Candidates"):
            os.makedirs("Data/Candidates")

        if (
            self.calculate_skewness() >= 5
        ):  # Flag for human inspection (found from https://arxiv.org/pdf/2312.04199)
            if not os.path.exists("Data/Candidates/Tier_1"):
                os.makedirs("Data/Candidates/Tier_1")
            return "Tier I"
        elif (self.calculate_skewness() > 0) & (self.calculate_skewness() < 5):
            if not os.path.exists("Data/Candidates/Tier_2"):
                os.makedirs("Data/Candidates/Tier_2")
            return "Tier II"
        else:
            if not os.path.exists("Data/Candidates/Tier_3"):
                os.makedirs("Data/Candidates/Tier_3")
            return "Tier III"

    def plot_lightcurve(self, tic_id):
        """
        Plots the high and low power spectra with the Normalized
        and phase-folded light curves.
        - Returns:
            - 4 sub-plots
        """

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))
        ax1, ax2, ax3, ax4 = axs.flatten()

        # Lower frequency periodogram
        ax1.plot(self.pg_low.frequency.value, self.mmag_lps, color="k", linewidth=1.0)
        ax1.set_xlim(0, max(self.pg_low.frequency.value) / 2)
        ax1.set_ylim(0, max(self.mmag_lps) * 1.05)

        # Add secondary x-axis for mHz
        secax = ax1.secondary_xaxis("top", functions=(self.cpd_to_mhz, self.mhz_to_cpd))
        secax.set_xlabel("Frequency (mHz)")
        secax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

        ax1.set_title("Low-Frequency Periodogram")
        ax1.set_xlabel(r"Frequency ($d^{-1}$)")
        ax1.set_ylabel("Amplitude (mmag)")
        ax1.set_title("Low-Frequency Periodogram", fontweight="bold")

        # Plotting the main periodogram with peaks and threshold
        ax3.plot(self.pg.frequency.value, self.mmag_ps, color="k", linewidth=1.0)

        ymax = ax3.get_ylim()[1]  # Get current y-axis maximum

        num, threshold, freq, peaks, spacing = self.find_peaks()
        skew = self.calculate_skewness()

        # Find the index of the central peak (highest peak)
        for i, (f, h) in enumerate(zip(freq, peaks)):
            if h == max(self.mmag_ps):  # Central peak or very far from peak
                ax3.vlines(
                    x=f,
                    ymin=h,
                    ymax=ymax * 1.1,
                    linestyle=":",
                    color="red",
                    label="Central Peak",
                    zorder=0,
                )  # Central
            else:
                ax3.vlines(
                    x=f,
                    ymin=h,
                    ymax=ymax * 1.1,
                    linestyle="--",
                    color="lightgrey",
                    linewidth=0.8,
                    zorder=0,
                )  # Side lobes
        ax3.axhline(
            y=threshold,
            color="blue",
            linestyle="--",
            label=f"Threshold ({num}" + r"$\sigma$)",
        )  # Threshold line

        # Add more ticks
        ax3.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))  # ~10 ticks on x-axis
        ax3.yaxis.set_major_locator(ticker.MaxNLocator(nbins=10))  # ~10 ticks on y-axis

        xmin, xmax = min(freq) - 2, max(freq) + 2
        ax3.set_xlim(xmin, xmax)
        ax3.set_ylim(0, ymax * 1.1)
        ax3.set_xlabel(r"Frequency ($d^{-1}$)")
        ax3.set_ylabel("Amplitude (mmag)")

        # Add secondary x-axis for mHz
        secax = ax3.secondary_xaxis("top", functions=(self.cpd_to_mhz, self.mhz_to_cpd))
        secax.set_xlabel("Frequency (mHz)")
        secax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

        ax3.plot([], [], color="none", label=f"Skew = {skew:.3f} mmag")

        ax3.legend(loc="best")
        ax3.set_title("High Frequency Periodogram", fontweight="bold")

        lc = self.light_curve.normalize()

        # Light Curve plot
        mmag_lc = -2.5 * np.log10(lc.flux / np.mean(lc.flux)) * 1000
        ax2.scatter(lc.time.value, mmag_lc, s=1, color="k")
        ax2.set_xlabel("Time (BTJD)")
        ax2.set_ylabel("Δ Magnitude (mmag)")
        ax2.invert_yaxis()  # Invert y-axis for magnitude
        ax2.set_title("Normalized Light Curve", fontweight="bold")

        # Bin the folded light curve to reduce scatter
        # Bin: specify bin size in *days*, so convert desired phase bin width to time
        folded_lc = lc.fold(period=self.rot_period)
        desired_phase_bin = 0.0005  # 1% of the phase
        bin_size_days = desired_phase_bin * self.rot_period

        binned = folded_lc.bin(time_bin_size=bin_size_days)

        # Convert relative flux to mmag
        mmag_fpd = -2.5 * np.log10(binned.flux / np.mean(binned.flux)) * 1000
        phase = (
            binned.time.value % self.rot_period
        ) / self.rot_period  # Phase from 0 to 1
        # mmag_fpd = -2.5 * np.log10(folded_lc.flux / np.mean(folded_lc.flux)) * 1000
        # phase = (folded_lc.time.value%rot_period) / rot_period  # Phase from 0 to 1
        ax4.scatter(phase, mmag_fpd, s=1, color="k")
        ax4.plot([], [], color="none", label=f"Rot Period: {self.rot_period:.3f} d")
        ax4.set_xlabel("Phase")
        ax4.set_ylabel("Δ Magnitude (mmag)")
        ax4.set_title(
            "Phase-folded Light Curve",
            fontweight="bold",
        )
        ax4.invert_yaxis()  # Invert y-axis for magnitude
        ax4.legend(loc="best")  # loc="center right", bbox_to_anchor=(1, 0.93)

        plt.suptitle(f"Analysis of {tic_id}", fontweight="bold")
        plt.tight_layout()

        if self.classify_candidates() == "Tier I":
            plt.savefig(f"Data/Candidates/Tier_1/analysis_{tic_id}.png", dpi=300)
            print(f"Saved {tic_id} as Tier I")
        elif self.classify_candidates() == "Tier II":
            plt.savefig(f"Data/Candidates/Tier_2/analysis_{tic_id}.png", dpi=300)
            print(f"Saved {tic_id} as Tier II")
        else:
            plt.savefig(f"Data/Candidates/Tier_3/analysis_{tic_id}.png", dpi=300)
            print(f"Saved {tic_id} as Tier III")

        # print(f"Peak Spacing: {spacing}")


if __name__ == "__main__":
    tic_list = [
        "TIC 387115314",
        "TIC 310817678",
        "TIC 467074220",
        "TIC 134631231",
    ]

    # tic_id = "TIC 387115314"
    # author = "TESS-SPOC"
    for id in tic_list:
        lc = lk.search_lightcurve(
            id, cadence="short"
        )  # , author="TESS-SPOC", sector=sectors)
        print(lc)

        try:
            if lc is None or len(lc) == 0:
                print(f"No light curves found for {id}.")
                exit()
            else:
                stitch = lc.download()
                # stitch = lc.download_all().stitch()
                clean_stitch = stitch.remove_nans()  # type: ignore
                print(
                    f"Stitched light curve for {id} with {len(clean_stitch)} data points."
                )
        except Exception as e:
            print(f"Error downloading or stitching light curves for {id}: {e}")
            exit()

        # When I get the dataset from Dan this will be replaced with the file input.
        time, flux = clean_stitch.time.value, clean_stitch.flux.value
        fr = find_roAp(flux, time)

        fr.plot_lightcurve(id)

    # print(fr.classify_candidates())
