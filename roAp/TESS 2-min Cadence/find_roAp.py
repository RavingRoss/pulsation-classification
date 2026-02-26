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
- Copied most of Test.py since the plotting tools are very good for that. But need
edit so that input data is 'time' and 'flux' (Data from sarek1).
- Using the automation tool from Test.py to find peaks in each periodogram.
    ~ I think it works, still testing.
"""

import lightkurve as lk

# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# import matplotlib.ticker as ticker
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import skew

# --------------------------------------------------------------------------------------------------------------


tic_id = "TIC 387115314"  # GOOD, fig 11
author = "TESS-SPOC"

lc = lk.search_lightcurve(
    tic_id, cadence="short"
)  # , author="TESS-SPOC", sector=sectors)
print(lc)

try:
    if lc is None or len(lc) == 0:
        print(f"No light curves found for {tic_id}.")
        exit()
    else:
        # stitch = lc.download()
        stitch = lc.download_all().stitch()
        clean_stitch = stitch.remove_nans()
        print(f"Stitched light curve for {tic_id} with {len(stitch)} data points.")
except Exception as e:
    print(f"Error downloading or stitching light curves for {tic_id}: {e}")
    exit()


time, flux = clean_stitch.time.value, clean_stitch.flux.value

pg = clean_stitch.flatten().to_periodogram(normalization="amplitude")
pg.show_properties()

pg_low = clean_stitch.to_periodogram(normalization="amplitude", maximum_frequency=5.0)

# To find the frequency of max peak, want to do for lower frequency spectrum to get rotation period
fmax = pg_low.frequency_at_max_power.to_value()
rot_period = 2 / fmax

# Threshold to find peaks
mmag_lps = (
    0.0010857362047581294 * 1e6 * pg_low.power.value
)  # Convert to mmag for lower frequency Periodogram
mmag_ps = (
    0.0010857362047581294 * 1e6 * pg.power.value
)  # Convert to mmag for main Periodogram

# Define conversion functions
factor_mHz = 1 / 86400 * 1000  # cycles/day → mHz


# --------------------------------------------------------------------------------------------------------------
# Where the fun begins:


class find_roAp:
    """
    Class to find roAp pulsations in TESS data using given data from sarek1 server which
    outputs the Flux and Time of stars. With this we search through the database to find
    roAp candidates.
    - Functions:
        - cpd_to_mhz(self): Converts cycles/day to mHz.
        - mhz_to_cpd(self): Converts mHz to cycles/day.
        - find_peaks(self): Finds peaks in the power spectrum.
        - calculate_skewness(self): Calculates the skewness of the power spectrum.
        - plot_lightcurve(self): Plots the high and low power spectra with the Normalized
                                 and phase-folded light curves.
    """

    def __init__(self, flux, time):
        """
        Initializes the class with the flux and time data.
        """
        self.flux = flux
        self.time = time

    def cpd_to_mhz(self, x):
        """
        Converts cycles/day to mHz.
        - Returns:
            - mHz (float): Frequency in mHz.
        """
        return x * factor_mHz

    def mhz_to_cpd(self, x):
        """
        Converts mHz to cycles/day.
        - Returns:
            - cpd (float): Frequency in cycles/day.
        """
        return x / factor_mHz

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
        mean_power = np.nanmean(mmag_ps)
        max_power = np.nanmax(mmag_ps)
        num = 5
        threshold = num * mean_power
        print(
            f"Max Power: {max_power} mmag | Mean Power: {mean_power} mmag | Threshold set to: {threshold} mmag"
        )

        indices, properties = find_peaks(mmag_ps, height=float(threshold))
        peaks = properties["peak_heights"]

        # Get the frequencies of the peaks
        peak_freqs = pg.frequency[indices].to_value()

        # Calculate the spacing between consecutive peaks
        spacing = np.diff(peak_freqs)

        for p in indices:
            f = pg.frequency[p].to_value()
            freq.append(f)

        return num, threshold, freq, peaks, spacing

    def calculate_skewness(self):
        """
        Calculates the skewness of the power spectrum.
        - Returns:
            - skew (float): Skewness of the power spectrum.
        """
        return skew(mmag_ps, nan_policy="omit")

    def plot_lightcurve(self):
        """
        Plots the high and low power spectra with the Normalized
        and phase-folded light curves.
        - Returns:
            - 4 plots
        """
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))
        ax1, ax2, ax3, ax4 = axs.flatten()

        # Lower frequency periodogram
        ax1.plot(pg_low.frequency.value, mmag_lps, color="k", linewidth=1.0)
        ax1.set_xlim(0, max(pg_low.frequency.value) / 2)
        ax1.set_ylim(0, max(mmag_lps) * 1.05)

        # Add secondary x-axis for mHz
        secax = ax1.secondary_xaxis("top", functions=(self.cpd_to_mhz, self.mhz_to_cpd))
        secax.set_xlabel("Frequency (mHz)")
        secax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

        ax1.set_title("Low-Frequency Periodogram")
        ax1.set_xlabel(r"Frequency ($d^{-1}$)")
        ax1.set_ylabel("Amplitude (mmag)")
        ax1.set_title("Low-Frequency Periodogram", fontweight="bold")

        # Plotting the main periodogram with peaks and threshold
        ax3.plot(pg.frequency.value, mmag_ps, color="k", linewidth=1.0)

        ymax = ax3.get_ylim()[1]  # Get current y-axis maximum

        num, threshold, freq, peaks, spacing = self.find_peaks()
        skew = self.calculate_skewness()

        # Find the index of the central peak (highest peak)
        for i, (f, h) in enumerate(zip(freq, peaks)):
            if h == max(mmag_ps):  # Central peak or very far from peak
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

        ax3.legend()
        ax3.set_title(f"Periodogram (Skew = {skew:.3f})", fontweight="bold")

        lc = clean_stitch.normalize()

        # Light Curve plot
        mmag_lc = -2.5 * np.log10(lc.flux / np.mean(lc.flux)) * 1000
        ax2.scatter(lc.time.value, mmag_lc, s=1, color="k")
        ax2.set_xlabel("Time (BTJD)")
        ax2.set_ylabel("Δ Magnitude (mmag)")
        ax2.invert_yaxis()  # Invert y-axis for magnitude
        ax2.set_title("Normalized Light Curve", fontweight="bold")

        # Bin the folded light curve to reduce scatter
        # Bin: specify bin size in *days*, so convert desired phase bin width to time
        folded_lc = lc.fold(period=rot_period)
        desired_phase_bin = 0.0005  # 1% of the phase
        bin_size_days = desired_phase_bin * rot_period

        binned = folded_lc.bin(time_bin_size=bin_size_days)

        # Convert relative flux to mmag
        mmag_fpd = -2.5 * np.log10(binned.flux / np.mean(binned.flux)) * 1000
        phase = (binned.time.value % rot_period) / rot_period  # Phase from 0 to 1
        # mmag_fpd = -2.5 * np.log10(folded_lc.flux / np.mean(folded_lc.flux)) * 1000
        # phase = (folded_lc.time.value%rot_period) / rot_period  # Phase from 0 to 1
        ax4.scatter(phase, mmag_fpd, s=1, color="k")
        ax4.set_xlabel("Phase")
        ax4.set_ylabel("Δ Magnitude (mmag)")
        ax4.set_title(
            f"Phase-folded Light Curve (Period = {(rot_period):.3f} d)",
            fontweight="bold",
        )
        ax4.invert_yaxis()  # Invert y-axis for magnitude

        plt.suptitle(f"Analysis of {tic_id}", fontweight="bold")
        plt.tight_layout()
        # plt.show()
        plt.savefig(f"analysis_{tic_id}.png", dpi=300)

        print(f"Peak Spacing: {spacing}")


if __name__ == "__main__":
    fr = find_roAp(flux, time)

    fr.plot_lightcurve()
    print("Finished compiling")
