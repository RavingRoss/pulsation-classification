'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
--------------------------------
Author: Jacob Romeo             |         
Updated: 05/05/2025             |   
Blue Straggler Classification Alg|    
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

# This is the VizierCatalog.py file that contains the functions to get the data from the clusters and members files
import VizierCatalog as VC
# This is the read_mist_models.py file from MIST which gives theoretical isochrones
import read_mist_models as rmm

# Importing packages as needed
from astroquery.mast import Catalogs
import matplotlib.patches as patches
from scipy.stats import gaussian_kde
from astroquery.vizier import Vizier 
from scipy.signal import find_peaks
import astropy.coordinates as coord 
import matplotlib.pyplot as plt
from astropy.table import Table 
from scipy.stats import skew
import astropy.units as u 
import lightkurve as lk
import pandas as pd
import numpy as np
import traceback
import time
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
    def plotData(self, roApFile, membersFile, clustersFile, target):
        
        try:
            cluster_mem = pd.read_parquet(membersFile)
            cluster = pd.read_parquet(clustersFile)
            #clusters = pd.read_csv(clustersFile)
            #roAp = pd.read_csv(roApFile)
            cluster_mem = cluster_mem[cluster_mem['ID'] == target]
            cluster = cluster[cluster['ID'] == target]
            c = pd.DataFrame(cluster_mem)
            #members = pd.read_csv(membersFile)
            # Getting the data from the clusters file
            gmag = cluster_mem['GMAG0'].values
            bp_rp = cluster_mem['BP-RP0'].values
            
            name = cluster['Name'].iloc[0]
            age = float(cluster['logAge50'].values)
            print('Cluster age:', age)
            
            isocmd = rmm.ISOCMD('Data/MIST_iso_reddened.iso.cmd') # Importing ISOCMD image for isochrones
            age_ind = isocmd.age_index(age) #returns the index for the desired age
            G = isocmd.isocmds[age_ind]['Gaia_G_DR2Rev']
            BP = isocmd.isocmds[age_ind]['Gaia_BP_DR2Rev']
            RP = isocmd.isocmds[age_ind]['Gaia_RP_DR2Rev']
            
            # Getting the data from the roAp file
            '''gmag_r = roAp['GMAG0'].values
            bp_rp_r = roAp['BP-RP0'].values'''
            
            fig, ax = plt.subplots()
            ax.scatter(bp_rp, gmag, s=5, c='m', label='Cluster', zorder=2)
            ax.plot(BP-RP, G+10, label=f'Isochrone Age: {age:.3f}', c='k', zorder=1) 
            #ax.scatter(bp_rp_r, gmag_r, s=5, c='r', label='roAp', zorder=2)
            ax.invert_yaxis()
            plt.title(f"Color Magnitude Diagram of {name}", fontweight='bold')
            ax.set_xlabel("BP-RP (Temperature)")
            ax.set_ylabel("GMAG (Absolute Magnitude)")
            ax.grid(True, zorder=0)
            
            # Define the region of interest (other values commented, uncomment if using de-reddened data from starhorse)
            xmin, xmax = min(bp_rp), 0.5 # 0.2 # x-axis range
            ymin, ymax = min(gmag), 13 # 2  # y-axis range
            
            # Add a rectangle to highlight the region
            rect = patches.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                                      linewidth=2, edgecolor='blue', facecolor='blue', alpha=0.3, zorder=3)
            ax.add_patch(rect)
            
            '''ax.set_xlim(xmin-1, max(bp_rp)+1)
            ax.set_ylim(max(gmag)+1, ymin+1)'''
            
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
                    df.to_csv(newpath, index=False)
                    print(f'Selected data saved to {newpath} of length {len(df)} vs original length {len(c)}')
            else:
                print('No region selected.')
        except Exception as e:
            print(f"\nException for target ID {target}: {type(e).__name__} — {e}")
            traceback.print_exc()
            
    # Vizier catalog query
    def vizierQuery(self, cat='J/A+A/686/A42', input=-1, # -1 is default for all rows
    col=['Name', 'ID', 'GaiaDR3', 'Prob', 'RA_ICRS', 'DE_ICRS','pmRA', 'pmDE', 'Gmag', 'BP-RP', 'dist50', 'logAge50', 'Mass50']): 
        """
        Quering the Vizier catalog for clusters and members data
        specified catalog (J/A+A/686/A42), specifying the columns with 'cols'.
        """
        start_time = time.time()
        try:
            V = Vizier(columns=col)
            catalog = V.find_catalogs(cat)
            V.ROW_LIMIT = input # or a large number
            print("Getting catalogs...")
            catalogs = V.get_catalogs(catalog.keys())
            print("Finished, converting and saving as csv files...")
            clusters = catalogs[0]
            clusters = clusters.to_pandas()
            clusters.to_csv('Data/clusters.csv', index=False)
            members = catalogs[1]
            members = members.to_pandas()
            members.to_csv('Data/members.csv', index=False)
            print(f'{len(clusters)} Clusters')
            print(clusters.columns)
            print(f'{len(members)} Members of clusters')
            print(members.columns)
        except Exception as e:
            print(f"Query failed: {e}")
            
        end = time.time()
        elapsed = (end - start_time) / 60
        print(f"Query completed in {elapsed:.2f} minutes, moving on...")
        #return clusters, members

class TimeDataTESS :
    """
    Class to get the TESS data from the lightkurve library 
    and plot the CMD, LC, and PWS of the TESS stars
    to analyze their variablity.
    """
    def getShortestCadence(self, ra, dec, local_file='Data/targets_120s.parquet'):
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
            return None, None, None, None

        search = lk.search_lightcurve(f"TIC {tic_id}", author='SPOC', cadence='short')#author='SPOC', mission='TESS', exptime=120)
        
        if len(search) == 0:
            print(f"No light curves found for TIC {tic_id}, (RA={ra}, DEC={dec})") # Debugging
            return None, None, None, None

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
            return None, None, None, None

        # Step 3: Download and stitch all valid entries
        try:
            downloaded = lk.LightCurveCollection([entry.download() for entry in valid_entries])
            print(f"Number of light curves downloaded: {len(downloaded)}")
            stitched = downloaded.stitch()
            stitched = stitched.normalize(unit='ppm')
            print('Stitched and Normalized light curve...')
            
            pg = lc.flatten().to_periodogram(normalization='amplitude')
            print(f"Converted to periodogram...")
            pg.show_properties()
            
            peaks, _ = find_peaks(pg.power, height=1e-4)
            
            if peaks.shape == (0,): # Filtering out the stars with no features
                print("No features to analyze...")
                return None, None, None, None
            
            skewed = skew(pg.power, nan_policy='omit')
            print('Skew value:', skewed)
            
            return stitched, pg, tic_id, skewed

        except Exception as e:
            print(f"Download failed for TIC {tic_id}: {e}")
            return None, None, None, None
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
            
            lc, pg, tic_id, skewed = self.getShortestCadence(ra, dec)
            
            if lc is None:
                continue

            # I am using Gaia's 'Table' because using extend or append does not work with how the data is formatted
            tbl = Table([lc.flux, lc.time], names=('Flux', 'Time (BKJD)'))
            tbl.write(f'Data/TESS Data/Light Curves/Light Curve TIC {tic_id}.csv', 
                      format='ascii.csv', overwrite=True)

            print(f"Downloaded data for TIC {tic_id}")

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
                    'pmRA' : row.pmRA,
                    'pmDEC' : row.pmDE,
                    'Nyquist Frequency': round(nyquist, 6),
                    'Max Power': round(maxPw, 6),
                    'Frequency at max Power': round(maxPwFreq, 4),
                    'GMAG': row.GMAG0,
                    'BP_RP': row.BP_RP0,
                    'Skew' : skewed
                })
                    
            # Returning the stars data as a pandas dataframe
            stars = pd.DataFrame(stars)
            print(stars)
            stars.to_csv(f"Data/TESS Data/Clusters/TESS Cluster {stars['Cluster Name'][0]}.csv")
            return stars, lc, pg
        
    def plotTESSData(self, targetID):
        """
        Plotting, literally everything, from the original CMD, now including
        the TESS object explicitely, to the power spectrum of the the TESS object.
        """
        tess_results, lc, pg = self.getTESSData(targetID)
        
        # Data from the cluster members dataset
        file1 = 'Data/members.parquet'
        clust_mem = pd.read_parquet(file1)
        file10 = 'Data/clusters.parquet'
        clust = pd.read_parquet(file10)
        clust_mem = clust_mem[clust_mem['ID'] == targetID]
        clust = clust[clust['ID'] == targetID]
        gmag0 = clust_mem['GMAG0']
        bp_rp0 = clust_mem['BP-RP0']
        prob0 = clust_mem['Prob']
        age = float(clust['logAge50'].values)
        # Data from the TESS stars
        stars = tess_results
            
        for i in range(len(stars)):
            # Getting the data from the TESS stars file
            nyquist = stars['Nyquist Frequency'][i]
            id = stars['TargetID'][i]
            gmag = stars['GMAG'][i]
            bp_rp = stars['BP_RP'][i]
            name = stars['Cluster Name'][i]
            prob = stars['Probability'][i]
            skewed = stars['Skew'][i]
            
            print(f'Plotting {id}...')
            
            # Making the subplots and setting the figsize
            fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 12))
            ax1, ax2, ax3, ax4 = axs.flatten()
            
            isocmd = rmm.ISOCMD('Data/MIST_iso_reddened.iso.cmd') # Importing ISOCMD image for isochrones
            age_ind = isocmd.age_index(age) #returns the index for the desired age
            G = isocmd.isocmds[age_ind]['Gaia_G_DR2Rev']
            BP = isocmd.isocmds[age_ind]['Gaia_BP_DR2Rev']
            RP = isocmd.isocmds[age_ind]['Gaia_RP_DR2Rev']
            
            # Plotting the CMD
            sc = ax1.scatter(bp_rp0, gmag0, c=prob0, cmap='viridis', label=f'{name}', s=10, zorder=1)
            ax1.scatter(bp_rp, gmag, c=prob, cmap='viridis', label=f'{id}', s=80, zorder=3, marker='*')
            ax1.plot(BP-RP, G+10, label=f'Isochrone Age: {age:.3f}', c='m', zorder=1) 
            ax1.plot([], [], ':', label=f'Skew: {skewed:.3f}')
            #ax1.plot([], [], '.', label=f'log10(age): {age:.2f}')
            fig.colorbar(sc, ax=ax1, label='Probability')
            ax1.set_title(f'CMD of TESS star with Cluster', fontweight='bold')
            ax1.set_ylabel('Absolute Mag [GMAG]', fontweight='bold')
            ax1.set_xlabel('Color/Temp [Bp-Rp]', fontweight='bold')
            ax1.set_xlim(min(bp_rp0)-0.5, max(bp_rp0)+0.5)
            ax1.set_ylim(max(gmag0)+1, min(gmag0)-1)
            ax1.grid(zorder=0)
            ax1.legend(loc='best')
            #ax1.invert_yaxis()
            
            # Plotting the light curve
            time = lc.time.value
            xmin, xmax = min(time), (min(time)+10)
            lc.plot(ax=ax3, label='')
            ax3.set_xlim(xmin, xmax)
            ax3.set_title(f'Limited Light Curve of {id}', fontweight='bold')
            ax3.set_ylabel('Flux [e/s]', fontweight='bold')
            ax3.set_xlabel('Time [BKJD]', fontweight='bold')
            ax3.grid(zorder=0)
            
            # Plotting the power spectrum
            pg.plot(ax=ax4, label='')
            ax4.set_title(f'Power Spectrum of {id}', fontweight='bold')
            ax4.set_xlabel('Frequency [1/d]', fontweight='bold')
            ax4.set_ylabel('')
            ax4.grid(zorder=0)
            
            # Plotting the power spectrum with limited x-axis to zoom in on features
            max_lim = pg.frequency_at_max_power+(40*(1/u.day))
            pg.plot(ax=ax2, label=f'Nyquist: {nyquist}')
            ax2.set_title(f'Limited Power Spectrum of {id}', fontweight='bold')
            ax2.set_ylabel('Power [ppm]', fontweight='bold')
            ax2.set_xlabel('')
            ax2.grid(zorder=0)
            ax2.set_xlim(0, max_lim.to_value())
            ax2.legend(loc='best')
            
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
        '''self.pfile = 'Data/cross-referenced-results.parquet'
        self.cfile = 'Data/cross-referenced-results.csv' '''
        self.cmfile = 'Data/members.csv'
        self.pmfile = 'Data/members.parquet'
        self.cfile = 'Data/clusters.csv'
        self.pfile = 'Data/clusters.parquet'
        self.file_pattern = "Data/TESS Data/Clusters/TESS Cluster *.csv"
        
        if os.path.exists(self.pfile) & os.path.exists(self.pmfile):
            print("Getting parquet from dir")
            self.pfile = self.pfile
            self.pmfile = self.pmfile
        else:
            print("Creating parquet files for faster iterating")
            par_m = pd.read_csv(self.cmfile)
            par_m.rename(columns={'Gmag': 'GMAG0'}, inplace=True)
            par_m.rename(columns={'BP-RP': 'BP-RP0'}, inplace=True)
            par_m.to_parquet(self.pmfile, index=False)
            self.pmfile = self.pmfile
            
            par_c = pd.read_csv(self.cfile)
            par_c.to_parquet(self.pfile, index=False)
            self.pfile = self.pfile
    
    # 'start' is where you left off and 'stop' variable is how many clusters there are in the file if queried all rows
    def run(self, start=1, stop=7167):
        start_time = time.time()
        try:
            for l in range(start-1, stop):
                start += 1
                targetID = start - 1
                # Running the GetData class to get the data from the cross-referenced results
                try:
                    print(f"\nProcessing target ID {targetID}...")
                    GetData().plotData(self.roAp, self.pmfile, self.pfile, targetID)
                    TimeDataTESS().plotTESSData(targetID)
                except Exception as e:
                    print(f"\nException for target ID {targetID}: {e}")
                    continue  # Skip to the next iteration
        except KeyboardInterrupt:
            print("\nCaught KeyboardInterrupt!")
        finally:
            print(f"Finished processing up to ID {targetID} out of {stop}.")
            # Concatenate the data from all the clusters' TESS results and graph the skew density plot
            self.plot_skew()
            end = time.time()
            elapsed = (end - start_time) / 3600
            print(f'Done in {elapsed:.2f} hours! Hope you got data lol!')
            # Reraise the KeyboardInterrupt to stop the program
            raise KeyboardInterrupt
    
    def plot_skew(self, path='Data/TESS Data/Skewed Density Plot'):
        """
        Plotting the skewed values for each candidate to generate
        a desnity plot in determining the best candidates.
        """
        results = self.concatenate_results()
        skewed = results['Skew']
        power = results['Max Power']
        # Calculate the 2D density
        xy = np.vstack([skewed, power])  # x and y are the same
        z = gaussian_kde(xy)(xy)
        
        fig, ax = plt.subplots(figsize=(10,6))
        sc = ax.scatter(skewed, power, c=z, s=20, cmap='viridis', label='Skew', zorder=2)
    
        cbar = plt.colorbar(sc, ax=ax, label='Density')
        cbar.set_ticks([])  # Remove the tick marks and numbers
        plt.xlabel('Skew value')
        plt.xlabel('Max Power [ppm]')
        plt.grid(zorder=0)
        plt.title("Density Plot to Determine 'Best' Candidates")
        plt.legend(loc='best')
        plt.savefig(path, dpi=300, bbox_inches='tight')
        #plt.show()
        
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
                
        
if __name__ == "__main__":
    # Only run if you want to cross-reference your cluster with starhorse, if not, import your own cluster file
    # from self.pfile or  self.cfile, parquet or csv resp, in the Main class __init__ function
    # If doing so, may need to change RA and DEC, GMAG and BP-RP, and other parameter column names in
    # TimeDataTESS class functions.
    start_time = time.time()
    '''print('Running the Vizier_Catalog.py file...')
    VC.GetData.main(input=-1) # -1 for all rows during the Vizier query, if testing set to lower number
    print('Finished cross-referencing, finding candidates and getting TESS data to analyze them...')
    '''
    # Change the rows var to the number of rows you want your input files to be, if testing set to lower number.
    #GetData().vizierQuery(input=-1)
    
    # Either set stop=some number, or stop=stop_clust to run through all of the clusters (found from vizierQuery function)
    clusters = pd.read_csv('Data/clusters.csv') 
    stop_clust = len(clusters) # will take a while, like 2 days, unless you change input value for Vizier query.

    #-------> CURRENTLY USING MEMBERS.CSV <-------
    # Stop variable to the target ID you want to start from, default is 1 and 7167, resp.
    Main().run(start=1593, stop=stop_clust) # Stopped at ID 1593