'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
--------------------------------
Author: Jacob Romeo             |         
Updated: 05/02/2025             |   
Blue Stargler Classification Alg|    
-> Cont. of .ipynb file         |
--------------------------------
       __.-._  ____ " Fear is the path to the dark side. 
      '-._"7' /             Fear leads to anger. 
       /'.-c                      Anger leads to hate.
       |  /T                             Hate leads to suffering." - Yoda, the wise coder
      _)_/LI

-> UPDATES <-
- Making the original .ipynb file into a .py file for easier access and readability.
    - This file will be used to classify roAp stars based on their time-dependent and magnitude properties.
- Still need to add Time-Dependent properties to the classification algorithm.
    - Possibly using TESS data to determine the time-dependent properties of the stars.
- Added TESS using lightkurve package, still testing with different clusters.
- Need to find specific clusters to look for Blue Stragglers in.
    - Used the J/A+A/686/A42 catalog from Vizier to find clusters and members instead of doing gaia query.
    - Then cross matched the results with the starhorse catalog (I/354) to get the de-reddened data.
- Added a function to plot the CMD, LC, and PWS of the TESS stars.
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# This is the Vizier_Catalog.py file that contains the functions to get the data from the clusters and members files
import VizierCatalog as VC

# Importing packages as needed
from astroquery.mast import Catalogs
import matplotlib.patches as patches
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
import lightkurve as lk
import pandas as pd
import numpy as np
import glob
import os

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')
        
class GetData:
    """
    Getting the data from the clusters and members files from Vizier 
    cross matching with starhorse, and returning the data as a 
    pandas dataframe to plot the CMD.
    """
    def plotData(self, roApFile, clustersFile, target):
        
        try:
            clusters = pd.read_parquet(clustersFile)
            #clusters = pd.read_csv(clustersFile)
            roAp = pd.read_csv(roApFile)
            clusters = clusters[clusters['ID'] == target]
            c = pd.DataFrame(clusters)
            #members = pd.read_csv(membersFile)
            # Getting the data from the clusters file
            gmag = clusters['GMAG0'].values
            bp_rp = clusters['BP-RP0'].values
            
            # Getting the data from the roAp file
            gmag_r = roAp['GMAG0'].values
            bp_rp_r = roAp['BP-RP0'].values
            
            '''gmag = clusters['Gmag'].values
            bp_rp = clusters['BPmag'].values - clusters['RPmag'].values'''
            
            fig, ax = plt.subplots()
            ax.scatter(bp_rp, gmag, s=5, c='m', label='Cluster', zorder=2)
            ax.scatter(bp_rp_r, gmag_r, s=5, c='r', label='roAp', zorder=2)
            ax.invert_yaxis()
            plt.title(f"Color Magnitude Diagram of ID {target}", fontweight='bold')
            ax.set_xlabel("BP-RP (Temperature)")
            ax.set_ylabel("GMAG (Absolute Magnitude)")
            ax.grid(True, zorder=0)
            
            # Define the region of interest (ROI)
            xmin, xmax = min(bp_rp), 0.2 # x-axis range
            ymin, ymax = min(gmag), 2   # y-axis range
            
            # Add a rectangle to highlight the region
            rect = patches.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                                      linewidth=2, edgecolor='blue', facecolor='blue', alpha=0.3, zorder=3)
            ax.add_patch(rect)
            
            plt.legend(loc='best')
            # Uncomment 'show' if you want to see the plot before selecting the region
            plt.savefig(f"Data/Candidates/Graphs/CMD of ID {target}", dpi=300, bbox_inches='tight')
            #plt.show()
            
            # Apply the condition to filter rows after the plot is closed
            if xmax is not None and ymax is not None:
                condition = ((gmag < ymax)) & ((bp_rp < xmax))
                print(f'\nxmin:{xmin},xmax:{xmax},\nymin:{ymin},ymax:{ymax}')
                df = c.loc[condition]

                if df.empty:
                    print(f"No data found in the selected region for target {target}.")
                else:
                    # Save the DataFrame back to the CSV file / Change based on what the input file is
                    newpath = f'Data/Candidates/stragler-candidates-ID{target}.csv'
                    # newpath = f'Data/1000 at 1 deg Candidates/stragler-candidates-ID{target}.csv'
                    # newpath = f'Data/300 at 0_5 deg Candidates/stragler-candidates-ID{target}.csv'
                    df.to_csv(newpath, index=False)
                    print(f'Selected data saved to {newpath} of length {len(df)} vs original length {len(c)}')
            else:
                print('No region selected.')
        except Exception as e:
            print(f"An error occurred for target ID {target}: {e}")

class TimeDataTESS :
    """
    Class to get the TESS data from the lightkurve library 
    and plot the CMD, LC, and PWS of the TESS stars
    to analyze their variablity.
    """
    def getShortestCadence(self, ra, dec, plot=False, local_file='Data/targets_120s.parquet'):
        """
        Using the RA and DEC, found from the Gaia query in 'GetData' class,
        of the star to get its lightcurve data from the shortest
        cadence available for that star, using the lightkurve library.
        """
        coords = coord.SkyCoord(ra, dec, unit="deg")
        ticResult = Catalogs.query_region(coords, radius=5 * u.arcsec, catalog="TIC")
        tic_id = ticResult['ID'][0]
        
        # Leaving this commented incase want to try again in the future, doing this returned 0 results
        '''# Reading the small cadence raw data from Dan Hey's GIT
        url = 'https://raw.githubusercontent.com/danhey/tess-atl/refs/heads/main/catalog/targets_120s.csv'

        # Check if the file exists locally
        if os.path.exists(local_file):
            print(f"Loading known targets from local file: {local_file}")
            known = pd.read_parquet(local_file)
        else:
            print(f"Downloading known targets from URL: {url}")
            known = pd.read_csv(url)
            # Save the file locally for future use
            os.makedirs(os.path.dirname(local_file), exist_ok=True)  # Create directory if it doesn't exist
            known.to_parquet(local_file, index=False)
            print(f"Known targets saved locally to: {local_file}")
        
        print(f'Cross referencing TIC {tic_id} with known data...')
        if tic_id in known['TICID']:'''

        if len(ticResult) == 0:
            print(f"No TIC match found near RA={ra}, DEC={dec}")
            return None, None

        search = lk.search_lightcurve(f"TIC {tic_id}", cadence='short')#author='SPOC', mission='TESS', exptime=120)
        
        if len(search) == 0:
            print(f"No light curves found for TIC {tic_id}, (RA={ra}, DEC={dec})") # Debugging
            return None, None

        print('Printing original:\n', search)
        
        # Step 1: Find the shortest cadence
        shortest_cadence = min([entry.exptime.value for entry in search if entry.exptime is not None])
        print(f"Shortest cadence found: {shortest_cadence} s")

        # Step 2: Filter entries with that exact cadence and valid target name
        valid_entries = []
        for entry in search:
            try:
                if entry.exptime is not None and entry.exptime.value == shortest_cadence:
                    lc = entry.download()
                    if not np.all(np.isnan(lc.flux.value)):
                        valid_entries.append(entry)
                        print(f"Valid entry added: {entry.exptime.value}")
                    else:
                        print(f"Skipped NaN flux entry: {entry.exptime.value}")
            except Exception as e:
                print(f"Error checking entry: {e}")
                continue

        if not valid_entries:
            print(f"No valid light curves found with cadence {shortest_cadence} s")
            return None, None

        # Step 3: Download and stitch all valid entries
        try:
            downloaded = lk.LightCurveCollection([entry.download() for entry in valid_entries])
            print(f"Number of light curves downloaded: {len(downloaded)}")
            stitched = downloaded.stitch()
            stitched = stitched.normalize(unit='ppm')
            print('Stitched and Normalized light curve...')

            if plot:
                stitched.plot()
            
            return stitched, tic_id

        except Exception as e:
            print(f"Download failed for TIC {tic_id}: {e}")
            return None, None
        '''else:
            print(f"TIC {tic_id} not found in known data.")
            return None, None'''

    def getTESSData(self, targetID):
        """
        Using the RA and DEC, found from the Gaia query in 'GetData' class,
        of the star to get its lightcurve data using the lightkurve library.
        """        
        # Reading the candidates from the csv file
        cands = pd.read_csv(f"Data/Candidates/stragler-candidates-ID{targetID}.csv")
        # Rename the column 'OldName' to 'NewName'
        cands.rename(columns={'BP-RP0': 'BP_RP0'}, inplace=True)
        stars = []
        
        print('Starting TESS data download...')
        i = 0
        for row in cands.itertuples(index=False):
            ra, dec = row.RA_ICRS, row.DE_ICRS
            prob = row.Prob
            i += 1
            print( '-'*40, f'\nIteration #{i}')
            
            # Skipping if the candidate has a low probability
            # This is a placeholder condition, adjust as needed
            if prob < 0.5:
                print(f"Skipping candidate with low probability: {prob}")
                continue
            
            lc, tic_id = self.getShortestCadence(ra, dec)
            
            if lc is None:
                continue

            # I am using Gaia's 'Table' because using extend or append does not work with how the data is formatted
            tbl = Table([lc.flux, lc.time], names=('Flux', 'Time (BKJD)'))
            tbl.write(f'Data/TESS Data/Light Curves/Light Curve TIC {tic_id}.csv', 
                      format='ascii.csv', overwrite=True)

            pg = lc.flatten().to_periodogram(normalization='amplitude')
            print(f"Downloaded data for TIC {tic_id}")
            pg.show_properties()

            id = pg.targetid
            nyquist = pg.nyquist.value
            maxPwFreq = pg.frequency_at_max_power.value
            maxPw = pg.max_power.value

            # I am using Gaia's 'Table' because using extend or append does not work with how the data is formatted
            tbl = Table([pg.frequency, pg.power], names=('Frequency', 'Power'))
            tbl.write(f'Data/TESS Data/Power Spectrums/Power Spectrum TIC {id}.csv', 
                      format='ascii.csv', overwrite=True)
            
            # Adding data to the stars list
            if (ra in cands['RA_ICRS'].values) and (dec in cands['DE_ICRS'].values):
                stars.append({
                    'Cluster Name': row.Name,
                    'ID': row.ID,
                    'TargetID': f"TIC {id}",
                    'GaiaDR3': row.GaiaDR3,
                    'Probability': prob,
                    'RA': ra,
                    'DEC': dec,
                    'Nyquist Frequency': round(nyquist, 6),
                    'Max Power': round(maxPw, 6),
                    'Frequency at max Power': round(maxPwFreq, 4),
                    'GMAG': row.GMAG0,
                    'BP_RP': row.BP_RP0
                })
                    
            # Returning the stars data as a pandas dataframe
            stars = pd.DataFrame(stars)
            print(stars)
            stars.to_csv(f"Data/TESS Data/Clusters/TESS Cluster {stars['Cluster Name'][0]}.csv")
            return stars
        
    def plotTESSData(self, targetID):
        """
        Plotting, literally everything, from the original CMD, now including
        the TESS object explicitely, to the power spectrum of the the TESS object.
        """
        tess_results = self.getTESSData(targetID)
        
        # Data from the roAp dataset
        file1 = 'Data/field_roAp_gaiaNstarhorse.csv'
        roAp = pd.read_csv(file1)
        gmag0 = roAp['GMAG0']
        bp_rp0 = roAp['BP-RP0']
        # Data from the cluster candidates
        file2 = f'Data/Candidates/stragler-candidates-ID{targetID}.csv'
        cands = pd.read_csv(file2)
        gmag1 = cands['GMAG0']
        bp_rp1 = cands['BP-RP0']
        # Data from the TESS stars
        stars = tess_results
            
        for i in range(len(stars)):
            # Getting the data from the TESS stars file
            maxPwFreq = stars['Frequency at max Power'][i]
            nyquist = stars['Nyquist Frequency'][i]
            id = stars['TargetID'][i]
            gmag = stars['GMAG'][i]
            bp_rp = stars['BP_RP'][i]
            name = stars['Cluster Name'][i]
            
            print(f'Plotting {id}...')
            
            # Making the subplots and setting the figsize
            fig, axs = plt.subplots(2, 2, figsize=(10, 12))
            
            # Plotting the CMD
            axs[0,0].scatter(bp_rp, gmag, c='black', label=f'{id}', s=30, zorder=3, marker='*')
            axs[0,0].scatter(bp_rp0, gmag0, c='m', label='roAps', s=10, zorder=1)
            axs[0,0].scatter(bp_rp1, gmag1, c='c', label=f'{name}', s=10, zorder=1)
            #axs[0].scatter(bp_rp1, gmag1, c='c', label='roAp', s=2, zorder=2)
            '''axs[0,0].set_xlim(bp_rp-0.2, bp_rp+0.2)
            axs[0,0].set_ylim(gmag-0.5, gmag+0.5)'''
            axs[0,0].set_title(f'CMD of TESS star with Cluster and roAp Stars', fontweight='bold')
            axs[0,0].set_ylabel('Absolute Mag [GMAG]', fontweight='bold')
            axs[0,0].set_xlabel('Color/Temp [Bp-Rp]', fontweight='bold')
            axs[0,0].grid(zorder=0)
            axs[0,0].legend(loc='best')
            axs[0,0].invert_yaxis()
            
            # Read the data from the light curve graphing profile
            file4 = f'Data/TESS Data/Light Curves/Light Curve {id}.csv'
            plot0 = Table.read(file4, format='ascii.csv')
            # Getting flux and time from the graphing profile
            flux = plot0['Flux']
            time = plot0['Time (BKJD)']
            
            # Plotting the light curve
            axs[1,0].plot(time, flux, c='k', zorder=3)
            axs[1,0].set_title(f'Light Curve of {id}', fontweight='bold')
            axs[1,0].set_ylabel('Flux [e/s]', fontweight='bold')
            axs[1,0].set_xlabel('Time [BKJD]', fontweight='bold')
            axs[1,0].grid(zorder=0)
            
            # Read the data from the power spectrum graphing profile
            file5 = f'Data/TESS Data/Power Spectrums/Power Spectrum {id}.csv'
            plot1 = Table.read(file5, format='ascii.csv')
            # Getting frequency and power from the graphing profile
            freq = plot1['Frequency']
            power = plot1['Power']
            
            # Plotting the power spectrum
            axs[1,1].plot(freq, power, c='k', zorder=3)
            axs[1,1].set_title(f'Power Spectrum of {id}', fontweight='bold')
            axs[1,1].set_xlabel('Frequency [1/d]', fontweight='bold')
            axs[1,1].grid(zorder=0)
            
            # Plotting the power spectrum with limited x-axis to zoom in on features
            axs[0,1].plot(freq, power, c='k', zorder=3, label=f'nyquist={nyquist} 1/d')
            axs[0,1].set_title(f'Limited Power Spectrum of {id}', fontweight='bold')
            axs[0,1].set_ylabel('Power [ppm]', fontweight='bold')
            axs[0,1].grid(zorder=0)
            axs[0,1].set_xlim(0, maxPwFreq + 5) # For limited graph
            #axs[0,1].set_xlim(0, 5) # For limited graph
            axs[0,1].legend(loc='best')
            
            plt.tight_layout()
            # Save the figure with a specific filename
            plt.savefig(f"Data/TESS Data/Graphs/CMD, LC, and PWS of {id}", dpi=300, bbox_inches='tight')
            #plt.show()    
        
class Main:
    """
    Main class to run the program and get the data from the clusters and members files
    from Vizier cross matching with starhorse, and returning the data as a pandas dataframe
    to plot the CMD, PWS, and LC.
    """
    def __init__(self):
        self.roAp = 'Data/field_roAp_gaiaNstarhorse.csv'
        self.pfile = 'Data/cross-referenced-results.parquet'
        self.cfile = 'Data/cross-referenced-results.csv'
        self.file_pattern = "Data/TESS Data/Clusters/TESS Cluster *.csv"
    
    def concatenate_results(self, output_file='Data/TESS Data/All_Candidates.csv'):
        """
        Concatenate all TESS Cluster result files into a single file.
        """
        # Get all CSV files in the directory
        file_pattern = self.file_pattern
        all_files = glob.glob(file_pattern)

        if not all_files:
            print("No result files found to concatenate.")
            return

        print(f"Found {len(all_files)} clusters to concatenate.")

        # Read and concatenate all files
        dataframes = []
        for file in all_files:
            try:
                df = pd.read_csv(file)
                # Drop the Unnamed column if it exists
                if 'Unnamed: 0' in df.columns:
                    df = df.drop(columns=['Unnamed: 0'])
                dataframes.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")
                continue

        # Combine all DataFrames into one
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df = combined_df.sort_values(by='ID')

        # Save the combined DataFrame to a single CSV file
        combined_df.to_csv(output_file, index=False)
        print(f"All candidates saved to {output_file}.")
        return combined_df
    
    # 'start' is where you left off and 'stop' variable is how many clusters there are in the file if queried all rows
    def run(self, start=1, stop=7167):
        for l in range(start-1, stop):
            start += 1
            targetID = start - 1
            # Running the GetData class to get the data from the cross-referenced results
            try:
                print(f"\nProcessing target ID {targetID}...")
                GetData().plotData(self.roAp, self.pfile, targetID)
                TimeDataTESS().plotTESSData(targetID)
            except Exception as e:
                print(f"An error occurred for target ID {targetID}: {e}")
                continue  # Skip to the next iteration
        print(f"Finished processing all targets up to ID {targetID}.")
        
        # Concatenate the data from all the clusters' TESS results
        results = self.concatenate_results()
        print(results.head(len(results)))

if __name__ == "__main__":
    # Only run if you want to cross-reference your cluster with starhorse, if not, import your own cluster file
    # from self.pfile or  self.cfile, parquet or csv resp, in the Main class __init__ 
    # If doing so, may need to change RA and DEC, GMAG and BP-RP, and other parameter column names in
    # TimeDataTESS class functions.
    print('Running the Vizier_Catalog.py file...')
    VC.GetData.main(input=-1) # -1 for all rows during the Vizier query, if testing set to lower number
    print('Finished cross-referencing, finding candidates and getting TESS data to analyze them...')
    
    # Either set stop=some number, or stop=stop_clust to run through all of the clusters
    clusters = pd.read_csv('Data/clusters.csv') 
    stop_clust = len(clusters) # (will take a while, like 2 days, unless you change input value for Vizier query)     
    
    Main().run(start=1, stop=stop_clust) # Change the stop variable to the target ID you want to start from, default is 1 and 7167, resp.
    print('Done! Hope you got data lol!')