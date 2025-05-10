import read_mist_models as rmm

from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from scipy.stats import skew
import astropy.units as u 
import lightkurve as lk
import pandas as pd
import numpy as np
import os

'''

Looking at the Candidates in more detail, doing another lightcurve search and looking at more parameters
to conclude if they are indeed a Blue Straggler, another type of pulsator, or neither.

'''
# The 'best-candidates.csv' file was made manually by inspecting the results from BlueStragglers.py

def peaks_test():
    
    path = 'Data/TESS Data/Graphs/CM-best-candidates.csv'
    cands = pd.read_csv(path)
    
    for i in range(9, 10):
        
        tic_id = cands['TargetID'][i]
        ra, dec = cands['RA'][i], cands['DEC'][i]
        
        print("Starting Light Curve search....")
        
        search = lk.search_lightcurve(f"{tic_id}", author='SPOC', cadence='short')
            
        print('Printing results:\n', search)
        
        lc = search.download().normalize(unit='ppm')
        time = lc.time.value
        xmin, xmax = min(time), (min(time)+10)
        lc.plot()
        plt.xlim(xmin, xmax)
        plt.title(f'{tic_id} Light Curve')
        plt.show()
        
        lc_all = search.download_all()
        lc_stitched = lc_all.stitch().normalize(unit='ppm')
        
        print(f"Downloaded data for TIC {tic_id}")
        
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
        
        pg = lc.flatten().to_periodogram(normalization='amplitude')
        
        x = pg.power
        peaks, _ = find_peaks(x, height=1e-4)
        print(peaks.shape)
        if peaks.shape == (0,): # Figure this out
            print("No features to analyze...")
        else:
            print("Plotting...")
            pg.plot(ax=ax1)
            max_lim = pg.frequency_at_max_power+(40*(1/u.day))
            ax1.set_xlim(0, max_lim.to_value())
            
            pg_all = lc_stitched.flatten().to_periodogram(normalization='amplitude')
            pg_all.show_properties()
            pg_all.plot(ax=ax2)
            ax2.set_xlim(0, max_lim.to_value())
            plt.suptitle(f'{tic_id} Power Spectrum')
            plt.show()

def skew_test():
    path = 'Data/TESS Data/Graphs/CM-best-candidates.csv'
    cands = pd.read_csv(path)
    
    if 'Skew LC' in cands.columns:
        print(f'Skewed values exist, loading results...')
        skewed = cands['Skew PG'].values
        skewlc = cands['Skew LC'].values
        maxpw = cands['Max Power'].values
        cand = cands.drop('Skew', axis=1)
        cand.to_csv(path, index=False)
        
        # Calculate the 2D density
        xy = np.vstack([skewlc, maxpw])
        z = gaussian_kde(xy)(xy)
        
        fig, ax = plt.subplots(figsize=(10,6))
        sc = ax.scatter(skewlc, maxpw, c=z, s=20, cmap='viridis', zorder=2)
       
        plt.colorbar(sc, ax=ax, label='Density')
        plt.grid(zorder=0)
        plt.show()
    
    else:
        xsigmas = []
        ysigmas = []
        print(f'Skewed values do not exist, creating new...')
        for i in range(0, len(cands)):
            
            tic_id = cands['TargetID'][i]
            ra, dec = cands['RA'][i], cands['DEC'][i]
            
            print("Starting Light Curve search....")
            
            search = lk.search_lightcurve(f"{tic_id}", author='SPOC', cadence='short')
                
            print('Printing results:\n', search)
            
            lc_all = search.download_all()
            lc_stitched = lc_all.stitch().normalize(unit='ppm')
            
            print(f"Downloaded data for TIC {tic_id}")
                    
            pg = lc_stitched.flatten().to_periodogram(normalization='amplitude')
            pg.show_properties()
            
            x = pg.power
            y = lc_stitched.flux.to_value()
            y = np.array(y)
            print(x)
            xsigma = skew(x, nan_policy='omit', keepdims=False)
            ysigma = skew(y, nan_policy='omit', keepdims=False)
            print('PG Skew value:', xsigma)
            xsigmas.append(xsigma)
            print('LC Skew value:', ysigma)
            ysigmas.append(ysigma)
        
        cands['PG Skew'] = xsigmas
        cands['LC Skew'] = ysigmas
        cands.to_csv(path, index=False)
        
        # Calculate the 2D density
        xy = np.vstack([ysigmas, cands['Max Power']])  # x and y are the same
        z = gaussian_kde(xy)(xy)
        
        fig, ax = plt.subplots()
        sc = ax.scatter(ysigmas, cands['Max Power'], c=z, s=20, cmap='viridis', zorder=2)
       
        plt.colorbar(sc, ax=ax, label='Density')
        plt.grid(zorder=0)
        plt.show()
    
def isochrone():
    
    isocmd = rmm.ISOCMD('Data/MIST_iso_reddened.iso.cmd')
    
    print ('version: ', isocmd.version)
    print ('photometric system: ', isocmd.photo_sys)
    print ('abundances: ', isocmd.abun)
    print ('rotation: ', isocmd.rot)
    print ('ages: ', [round(x,2) for x in isocmd.ages])
    print ('number of ages: ', isocmd.num_ages)
    print ('available columns: ', isocmd.hdr_list)
    
    age_ind = isocmd.age_index(8.5) #returns the index for the desired age
    G = isocmd.isocmds[age_ind]['Gaia_G_DR2Rev']
    BP = isocmd.isocmds[age_ind]['Gaia_BP_DR2Rev']
    RP = isocmd.isocmds[age_ind]['Gaia_RP_DR2Rev']
    plt.scatter(cands['BP_RP'], cands['GMAG'], label='Cands', c='r')
    plt.plot(BP-RP, G, label='Isochrone') 
    plt.xlabel('Bessell B - Bessell R')
    plt.ylabel('Bessell V')
    plt.gca().invert_yaxis()
    plt.legend(loc='best')
    plt.show()
    
if __name__ == '__main__':
    #peaks_test()
    skew_test()
   #isochrone()