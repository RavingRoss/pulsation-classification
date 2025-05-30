from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.stats import skew
import lightkurve as lk
import numpy as np
# TIC 335457083, weird periodogram at 200 s exptime
tic_id = 'TIC 445493624'
lc = lk.search_lightcurve(tic_id, exptime=120, author='SPOC')#, author='SPOC')
print(lc)
lc = lc.download_all()
stitch = lc.stitch().normalize(unit='ppm')
#stitch = lc.normalize(unit='ppm')
time = stitch.time.value

pg = stitch.flatten().to_periodogram(normalization='amplitude')
pg.show_properties()

# Automated threshold to find peaks
threshold = float(0.45 * np.max(pg.power))
print(f'Threshold: {threshold}')
indices, properties = find_peaks(pg.power, prominence=3e-5)
peaks = properties['prominences']

# Get the frequencies of the peaks
peak_freqs = pg.frequency[indices].to_value()

# Calculate the spacing between consecutive peaks
spacing = np.diff(peak_freqs)

skew = skew(pg.power, nan_policy='omit')
print('Skew value:', skew)

freq = []
for p in indices:
    f = pg.frequency[p].to_value()
    freq.append(f)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 5))
ax1, ax2 = axs.flatten()

pg.plot(ax=ax1, color='black', zorder=2)
ymax = ax1.get_ylim()[1]  # Get current y-axis maximum
# Find the index of the central peak (highest peak)
for i, (f, h) in enumerate(zip(freq, peaks)):
    if i == len(freq) // 2:
        ax1.vlines(x=f, ymin=h, ymax=ymax, linestyle='--', color='lightgrey', zorder=0)  # Central
    else:
        ax1.vlines(x=f, ymin=h, ymax=ymax, linestyle=':', color='lightgrey', zorder=0)   # Side lobes
        
# Plot spacing arrows and labels
'''for i in range(len(freq) - 1):
    x0, x1 = freq[i], freq[i+1]
    y = ymax * 0.95  # Place arrow near the top
    ax1.annotate(
        '', xy=(x1, y), xytext=(x0, y),
        arrowprops=dict(arrowstyle='<->', color='red', lw=1)
    )
    offset = 0.1 * (x1 - x0)  # 10% of the spacing to the right of x1
    ax1.text(x1 + offset, y * 1.06, f"{spacing[i]:.3f}", color='red', ha='center', va='bottom', fontsize=8)'''

ax1.set_xlim(min(freq)-2, max(freq)+2)
ax1.set_title(f'Periodogram (Skew: {skew:.3f})')

stitch.plot(ax=ax2, color='black', zorder=2)
ax2.set_xlim(min(time), min(time) + 10)  # Adjust x-axis limits as needed
ax2.set_title(f'Light Curve')

plt.suptitle(f'Analysis of {tic_id}')
plt.show()