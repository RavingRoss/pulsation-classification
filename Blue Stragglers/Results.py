import matplotlib.pyplot as plt
import astropy.units as u 
import lightkurve as lk
import pandas as pd
import numpy as np
import shutil
import os

cands = pd.read_csv('Data/TESS Data/All_Candidates.csv')

power = cands['Max Power (ppm)']
skew = cands['PG Skew']

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(8,6))

ax1.hist(np.log10(power), color='red')
ax2.hist(np.log10(skew), color='blue')
#plt.show()

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
    skew = row['PG Skew']
    
    '''filename = f"Final Plot of {target_id}.png"
    source_path = os.path.join("Data/TESS Data/Graphs", filename)  # your original plot location
    
    if not os.path.exists(source_path):
        print(f"File not found: {source_path}")
        continue'''

    if np.log10(skew) < 1:
        print('log10(skew) < 1:')
        print(target_id, skew)
        low_end.append(row)
        #shutil.move(source_path, os.path.join(folders['low_end'], filename))
    if (np.log10(skew) > 1) & (np.log10(skew) < 1.5):
        print('log10(skew) > 1 and < 1.5:')
        print(target_id, skew)
        mid_end.append(row)
        #shutil.move(source_path, os.path.join(folders['mid_end'], filename))
    if (np.log10(skew) > 1.5) & (np.log10(skew) < 2):
        print('log10(skew) > 1.5 and < 2:')
        print(target_id, skew)
        high_end.append(row)
        #shutil.move(source_path, os.path.join(folders['high_end'], filename))
        
print(f'Low end: {len(low_end)/len(cands)*100:.3f}%')
print(f'Mid end: {len(mid_end)/len(cands)*100:.3f}%')
print(f'High end: {len(high_end)/len(cands)*100:.3f}%')
        