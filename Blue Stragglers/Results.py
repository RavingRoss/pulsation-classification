import matplotlib.pyplot as plt
import astropy.units as u 
import lightkurve as lk
import pandas as pd

'''

Looking at the Candidates in more detail, doing another lightcurve search and looking at more parameters
to conclude if they are indeed a Blue Straggler, another type of pulsator, or neither.

'''
# The 'best-candidates.csv' file was made manually by inspecting the results from BlueStragglers.py
cands = pd.read_csv('Data/TESS Data/Graphs/best-candidates.csv')

for i in range(0, len(cands)):
    
    tic_id = cands['TargetID'][i]
    ra, dec = cands['RA'][i], cands['DEC'][i]
    
    print("Starting Light Curve search....")
    
    search = lk.search_lightcurve(f"{tic_id}", author='SPOC', cadence='short')
        
    print('Printing results:\n', search)
    
    lc = search.download().normalize(unit='ppm')
    lc.plot()
    plt.title(f'{tic_id} Light Curve')
    plt.show()
    
    lc_all = search.download_all()
    lc_stitched = lc_all.stitch().normalize(unit='ppm')
    
    print(f"Downloaded data for TIC {tic_id}")
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    
    pg = lc.flatten().to_periodogram(normalization='amplitude')
    pg.plot(ax=ax1)
    max_lim = pg.frequency_at_max_power+(50*(1/u.day))
    ax1.set_xlim(0, max_lim.to_value())
    
    pg_all = lc_stitched.flatten().to_periodogram(normalization='amplitude')
    pg_all.show_properties()
    pg_all.plot(ax=ax2)
    ax2.set_xlim(0, max_lim.to_value())
    plt.suptitle(f'{tic_id} Power Spectrum')
    plt.show()