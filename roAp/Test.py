from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import skew
import lightkurve as lk
import numpy as np

"""tic_id = 'TIC 233200244'
author = "TESS-SPOC"
sectors = 57,56,60,58,59"""

# tic_id = 'TIC 165052884' # WEIRD

# tic_id = 'TIC 101624823' # NONE at 200s

'''tic_id = 'TIC 259017938'
author = "TESS-SPOC"'''

"""tic_id = 'TIC 335457083' # NOISE
author = "TESS-SPOC"
sectors = 73#60"""

'''tic_id = 'TIC 445493624' # Lots of noise
author="TESS-SPOC"'''

"""tic_id = 'TIC 467074220' # Lots of noise
author="TESS-SPOC"
sectors = 73#, 59"""

tic_id = "TIC 310817678"  # GOOD, fig 8
author = "TESS-SPOC"

'''tic_id = 'TIC 387115314' # GOOD, fig 11
author = "TESS-SPOC"'''

lc = lk.search_lightcurve(tic_id, exptime=200, author=author)  # , sector=sectors)
print(lc)

try:
    if lc is None or len(lc) == 0:
        print(f"No light curves found for {tic_id}.")
        exit()
    else:
        # stitch = lc.download()
        stitch = lc.download_all().stitch()
        print(f"Stitched light curve for {tic_id} with {len(stitch)} data points.")
except Exception as e:
    print(f"Error downloading or stitching light curves for {tic_id}: {e}")
    exit()

time = stitch.time.value

pg = stitch.flatten().to_periodogram(normalization="amplitude")
pg.show_properties()

pg_low = stitch.to_periodogram(normalization="amplitude", maximum_frequency=10.0)

# To find the frequency of max peak, want to do for lower frequency spectrum to get rotation period
pmax = pg_low.frequency.value[np.argmax(pg_low.power.value)]
rot_period = 2 / pmax
print(f"Rotation Period: {1/pmax} days")

# Threshold to find peaks
mmag_lps = (
    0.0010857362047581294 * 1e6 * pg_low.power.value
)  # Convert to mmag for lower frequency Periodogram
mmag_ps = (
    0.0010857362047581294 * 1e6 * pg.power.value
)  # Convert to mmag for main Periodogram

# Define conversion functions
factor_mHz = 1 / 86400 * 1000  # cycles/day → mHz


def cpd_to_mhz(x):
    return x * factor_mHz


def mhz_to_cpd(x):
    return x / factor_mHz


mean_power = np.nanmean(mmag_ps)
num = 5
threshold = num * mean_power
print(f"Mean Power: {mean_power}, Threshold set to: {threshold}")

indices, properties = find_peaks(mmag_ps, height=float(threshold))
peaks = properties["peak_heights"]

# Get the frequencies of the peaks
peak_freqs = pg.frequency[indices].to_value()

# Calculate the spacing between consecutive peaks
spacing = np.diff(peak_freqs)

skew_val = skew(mmag_ps, nan_policy="omit")

freq = []
for p in indices:
    f = pg.frequency[p].to_value()
    freq.append(f)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 8))
ax1, ax2, ax3, ax4 = axs.flatten()

# Lower frequency periodogram
ax1.plot(pg_low.frequency.value, mmag_lps, color="k", linewidth=1.0)
ax1.set_xlim(0, max(pg_low.frequency.value) / 2)
ax1.set_ylim(0, max(mmag_lps) * 1.05)

# Add secondary x-axis for mHz
secax = ax1.secondary_xaxis("top", functions=(cpd_to_mhz, mhz_to_cpd))
secax.set_xlabel("Frequency (mHz)")
secax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

ax1.set_title("Low-Frequency Periodogram")
ax1.set_xlabel(r"Frequency ($d^{-1}$)")
ax1.set_ylabel("Amplitude (mmag)")
ax1.set_title("Low-Frequency Periodogram", fontweight="bold")

# Plotting the main periodogram with peaks and threshold
ax3.plot(pg.frequency.value, mmag_ps, color="k", linewidth=1.0)

# pg.plot(ax=ax3, color='black', zorder=2)
ymax = ax3.get_ylim()[1]  # Get current y-axis maximum
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
    """elif i < len(spacing) and spacing[i] > 0.5:
        ax1.vlines(x=f, ymin=h, ymax=ymax, linestyle=':', color='lightgrey', zorder=0)"""  # Central (by spacing)

ax3.axhline(
    y=threshold, color="blue", linestyle="--", label=f"Threshold ({num}" + r"$\sigma$)"
)  # Threshold line

# Plot spacing arrows and labels
"""for i in range(len(freq) - 1):
    x0, x1 = freq[i], freq[i+1]
    y = ymax * 0.95  # Place arrow near the top
    ax1.annotate(
        '', xy=(x1, y), xytext=(x0, y),
        arrowprops=dict(arrowstyle='<->', color='red', lw=1)
    )
    offset = 0.1 * (x1 - x0)  # 10% of the spacing to the right of x1
    ax1.text(x1 + offset, y * 1.06, f"{spacing[i]:.3f}", color='red', ha='center', va='bottom', fontsize=8)"""

# Add more ticks
ax3.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))  # ~10 ticks on x-axis
ax3.yaxis.set_major_locator(ticker.MaxNLocator(nbins=10))  # ~10 ticks on y-axis

xmin, xmax = min(freq) - 2, max(freq) + 2
ax3.set_xlim(xmin, xmax)
ax3.set_ylim(0, ymax * 1.1)
ax3.set_xlabel(r"Frequency ($d^{-1}$)")
ax3.set_ylabel("Amplitude (mmag)")

# Add secondary x-axis for mHz
secax = ax3.secondary_xaxis("top", functions=(cpd_to_mhz, mhz_to_cpd))
secax.set_xlabel("Frequency (mHz)")
secax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

ax3.legend()
ax3.set_title(f"Periodogram (Skew = {skew_val:.3f})", fontweight="bold")

stitch = stitch.normalize()

# Light Curve plot
mmag_lc = -2.5 * np.log10(stitch.flux / np.mean(stitch.flux)) * 1000
ax2.scatter(stitch.time.value, mmag_lc, s=1, color="k")
ax2.set_xlabel("Time (BTJD)")
ax2.set_ylabel("Δ Magnitude (mmag)")
ax2.invert_yaxis()  # Invert y-axis for magnitude
# lc_plot = stitch.plot(ax=ax2, color='black', zorder=2)
ax2.set_title("Normalized Light Curve", fontweight="bold")

# Bin the folded light curve to reduce scatter
# Bin: specify bin size in *days*, so convert desired phase bin width to time
folded_lc = stitch.fold(period=rot_period)
desired_phase_bin = 0.0005  # 1% of the phase
bin_size_days = desired_phase_bin * rot_period

binned = folded_lc.bin(time_bin_size=bin_size_days)

# Convert relative flux to mmag
mmag_fpd = -2.5 * np.log10(folded_lc.flux / np.mean(folded_lc.flux)) * 1000
ax4.scatter(folded_lc.phase.value, mmag_fpd, s=1, color="k")
ax4.set_xlabel("Phase")
ax4.set_ylabel("Δ Magnitude (mmag)")
# folded_lc.plot(ax=ax4, color='black', zorder=2)
ax4.set_title(
    f"Phase-folded Light Curve (Period = {(1/pmax):.3f} d)", fontweight="bold"
)
ax4.invert_yaxis()  # Invert y-axis for magnitude

plt.suptitle(f"Analysis of {tic_id}", fontweight="bold")
plt.tight_layout()
plt.show()

print(f"Peak Spacing: {spacing}")
