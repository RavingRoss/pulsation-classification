import matplotlib.patches as patches
import matplotlib.pyplot as plt
import astropy.units as u 
import lightkurve as lk
import pandas as pd
import numpy as np
import shutil
import os

cands = pd.read_csv('Data/TESS Data/All_Candidates.csv')

power = cands['Max Power (ppm)']
skewed = cands['PG Skew']

'''
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8,6))
ax1.hist(np.log10(power), color='red')
ax2.hist(np.log10(skewed), color='blue')
plt.show()
'''
# Make sure folders exist
folders = {
    'low_end': 'Data/TESS Data/Graphs/Low_End',
    'mid_end': 'Data/TESS Data/Graphs/Mid_End',
    'high_end': 'Data/TESS Data/Graphs/High_End',
}

for folder in folders.values():
    os.makedirs(folder, exist_ok=True)

low_end = []
mid_end = []
high_end = []

for _, row in cands.iterrows():
    target_id = row['TargetID']
    pgskew = row['PG Skew']
    lcskew = row['LC Skew']
    
    '''filename = f"Final Plot of {target_id}.png"
    source_path = os.path.join("Data/TESS Data/Graphs", filename)  # your original plot location
    
    if not os.path.exists(source_path):
        print(f"File not found: {source_path}")
        continue'''

    if np.log10(pgskew) < 1:
        low_end.append(row)
        #shutil.move(source_path, os.path.join(folders['low_end'], filename))
    if (np.log10(pgskew) > 1) & (np.log10(pgskew) < 1.5):
        mid_end.append(row)
        #shutil.move(source_path, os.path.join(folders['mid_end'], filename))
    if (np.log10(pgskew) > 1.5) & (np.log10(pgskew) < 2):
        high_end.append(row)
        #shutil.move(source_path, os.path.join(folders['high_end'], filename))
        
print(f'Low end: {len(low_end)/len(cands)*100:.3f}%')
print(f'Mid end: {len(mid_end)/len(cands)*100:.3f}%')
print(f'High end: {len(high_end)/len(cands)*100:.3f}%')
print(f'Total Number: {len(cands)}')

middle = pd.DataFrame(mid_end)
best = pd.DataFrame(high_end)
last = pd.DataFrame(low_end)

all = pd.concat((best, middle), ignore_index=True)

best = best.sort_values(by='Frequency at max Power (1/d)', ascending=False)
middle = middle.sort_values(by='Frequency at max Power (1/d)', ascending=False)
all = all.sort_values(by='PG Skew', ascending=False)

fig, ax = plt.subplots(figsize=(8,6))

# Filter rows where mass > 3
filtered = all[all['Solar Mass'] > 2.5].reset_index(drop=True)
best = best[best['Solar Mass'] > 2.5].reset_index(drop=True)
middle = middle[middle['Solar Mass'] > 2.5].reset_index(drop=True)

print(f'High/Mid-Tier Skew: {len(all)/len(cands)*100:.3f}%')
print(f'Candidates > 3 Solar Masses: {len(filtered)/len(all)*100:.3f}%')
filtered.to_csv('Data/TESS Data/Massive_Candidates.csv', index=False)

best_freq = best['Frequency at max Power (1/d)']
mid_freq = middle['Frequency at max Power (1/d)']
best_mass = best['Solar Mass']
mid_mass = middle['Solar Mass']

ax.scatter(best_freq, best_mass, c='m', label='Best Tier', zorder=2)
ax.scatter(mid_freq, mid_mass, c='c', label='Mid Tier', zorder=2)

ymin, ymax = ax.get_ylim()
rect1 = patches.Rectangle((0.5, ymin), 4 - 0.5, ymax - ymin,
linewidth=2, edgecolor='black', facecolor='black', alpha=0.1, zorder=3, label=r'$\gamma$ Dor')
rect2 = patches.Rectangle((5, ymin), 40 - 5, ymax - ymin,
linewidth=2, edgecolor='green', facecolor='green', alpha=0.1, zorder=3, label=r'$\delta$ Scuti')
ax.add_patch(rect1)
ax.add_patch(rect2)

ax.set_ylabel('Solar Mass', fontsize=14)
ax.set_xlabel('Max Frequency [1/d]', fontsize=14)
ax.set_title(r'BSS Candidates Above $2.5M_{\odot}$', fontsize=18)
ax.grid(True, zorder=0) 
ax.legend(loc='best')
plt.savefig('Data/TESS Data/Mass vs Freq')
#plt.show()

best_cands = pd.read_csv('Data/TESS Data/Graphs/Best/Best_Candidates.csv') 

print(f"Best out of 180: {len(best_cands)/len(cands)*100:.3f}%")