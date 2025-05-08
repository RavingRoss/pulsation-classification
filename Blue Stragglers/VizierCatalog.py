from requests.exceptions import ChunkedEncodingError, ReadTimeout 
from astropy.coordinates import SkyCoord 
from astroquery.vizier import Vizier 
from astropy import units as u 
import pandas as pd 
import time
import os

class GetData:
    """
    Class to get data from catalog J/A+A/686/A42 using Vizier
    to cross-reference stars from the StarHorse catalog.
    """
    # Vizier catalog query
    def vizierQuery(self, cat='J/A+A/686/A42', input=-1): # -1 is default for all rows
        """
        Quering the Vizier catalog for clusters and members data
        specified catalog (J/A+A/686/A42).
        """
        V = Vizier()
        
        catalog = V.find_catalogs(cat)
        V.ROW_LIMIT = input # or a large number
        catalogs = V.get_catalogs(catalog.keys())
        clusters = catalogs[0]
        clusters = clusters.to_pandas()
        clusters.to_csv('Data/clusters.csv', index=False)
        members = catalogs[1]
        members = members.to_pandas()
        members.to_csv('Data/members.csv', index=False)
        print(f'{len(members)} Clusters')
        print(f'{len(members)} Members of clusters')
        return clusters, members

    # Query StarHorse catalog
    def queryStarhorse(self, ra, dec, radius=1):
        """
        Query StarHorse catalog using RA, Dec, and radius (in degrees).
        """
        V = Vizier()
        V.ROW_LIMIT = -1  # No row limit
        # Create a SkyCoord object for the coordinates
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        
        print(f'Querying StarHorse catalog for RA: {ra}, Dec: {dec}, Radius: {radius} degrees...')
        
        maxRetries = 3
        retryDelay = 10  # seconds
        for attempt in range(1, maxRetries+1):
            try: 
                result = V.query_region(coord, radius=f"{radius} deg", catalog="I/354")
                if len(result) > 0:
                    starhorse_data = result[0].to_pandas()
                    
                    # Save the result incrementally to a CSV file
                    output_file = 'Data/starhorse_incremental.csv'
                    starhorse_data.to_csv(output_file, mode='a', header=not os.path.exists(output_file), index=False)
                    
                    return starhorse_data
                else:
                    print(f"No data found for RA: {ra}, Dec: {dec}. Moving on...")
                    return pd.DataFrame()
                break
            except (TimeoutError, ReadTimeout, ChunkedEncodingError, ConnectionError) as e:
                print(f"TimeoutError: {e}. Retrying after {retryDelay} seconds...")
                time.sleep(retryDelay)
        # If all retries fail, log the failure and move on
        print(f"Failed after {maxRetries} attempts for RA: {ra}, Dec: {dec}. Moving on...")
        return pd.DataFrame()

    # Cross-reference clusters and members with Gaia and StarHorse
    def crossReference(self, clusters, members):
        starhorse_results_file = 'Data/starhorse_incremental.csv'
        cross_referenced_file = 'Data/cross-referenced-results.csv'
        final_starhorse_file = 'Data/starhorse.csv'

        # Remove the cross-referenced results file if it exists (optional, for fresh runs)
        '''if os.path.exists(cross_referenced_file):
            os.remove(cross_referenced_file)'''

        # Load members and clean RA/DEC columns
        members = pd.DataFrame(members)
        members['RA_ICRS'] = pd.to_numeric(members['RA_ICRS'], errors='coerce')
        members['DE_ICRS'] = pd.to_numeric(members['DE_ICRS'], errors='coerce')

        # Create SkyCoord object for members
        members_coords = SkyCoord(ra=members['RA_ICRS'], 
                                dec=members['DE_ICRS'],
                                unit=(u.deg, u.deg), 
                                frame='icrs')

        # Check if either file exists
        if not os.path.exists(starhorse_results_file) and not os.path.exists(final_starhorse_file):
            print("No StarHorse files found. Querying StarHorse catalog...")
            # Query StarHorse for each cluster
            for _, row in clusters.iterrows():
                ra, dec = row['RA_ICRS'], row['DE_ICRS']
                starhorse_data = self.queryStarhorse(ra, dec)

                if not starhorse_data.empty:
                    # Save the result incrementally to the final file
                    starhorse_data.to_csv(final_starhorse_file, mode='a', 
                                        header=not os.path.exists(final_starhorse_file), index=False)

            # Combine results into a single DataFrame (optional, only if needed later)
            '''if os.path.exists(final_starhorse_file):
                starhorse_results = pd.read_csv(final_starhorse_file)
                print(f"StarHorse data loaded from {final_starhorse_file} of length {len(starhorse_results)}.")
            else:
                starhorse_results = pd.DataFrame()
                print("No StarHorse data retrieved.")'''

        # Process StarHorse data in chunks for cross-referencing
        print('Processing StarHorse data in chunks for cross-referencing...')
        file_to_process = starhorse_results_file if os.path.exists(starhorse_results_file) else final_starhorse_file
        chunk_size = 0
        try:
            for chunk in pd.read_csv(file_to_process, chunksize=1000000):
                # Clean RA/DEC columns in the chunk
                chunk['RA_ICRS'] = pd.to_numeric(chunk['RA_ICRS'], errors='coerce')
                chunk['DE_ICRS'] = pd.to_numeric(chunk['DE_ICRS'], errors='coerce')

                # Create SkyCoord object for the chunk
                starhorse_coords = SkyCoord(ra=chunk['RA_ICRS'], 
                                            dec=chunk['DE_ICRS'],
                                            unit=(u.deg, u.deg), 
                                            frame='icrs')

                # Perform cross-matching
                chunk_size += chunk
                print(f'Cross-Referencing chunk {len(chunk_size)}...')
                idx, d2d, _ = members_coords.match_to_catalog_sky(starhorse_coords)

                # Add the angular separation to the members DataFrame
                members['match_idx'] = idx
                members['separation'] = d2d.arcsec

                # Filter matches within a threshold (e.g., 1 arcsecond)
                match_threshold = 1.0  # in arcseconds
                matched_members = members[members['separation'] <= match_threshold].copy()

                # Add matched StarHorse data to the members DataFrame
                matched_starhorse = chunk.iloc[matched_members['match_idx']].reset_index(drop=True)
                merged_results = pd.concat([matched_members.reset_index(drop=True), matched_starhorse], axis=1)
                
                # Save the merged results incrementally
                merged_results.to_csv(cross_referenced_file, mode='a', 
                                    header=not os.path.exists(cross_referenced_file), index=False)
        except Exception as e:
            print(f"Error processing StarHorse data in chunks: {e}")

    # Main function to run the class
    def main(self, input=-1): # -1 for all rows
        # Using the Vizier catalog to get clusters and members
        clustersdf = pd.read_csv('Data/clusters.csv') if os.path.exists('Data/clusters.csv') else pd.DataFrame()
        membersdf = pd.read_csv('Data/members.csv') if os.path.exists('Data/members.csv') else pd.DataFrame()

        if clustersdf.empty and membersdf.empty:
            print('Querying Vizier members and cluster catalogs...')
            clusters, members = self.vizierQuery(input=input)
        else:
            clusters = clustersdf
            members = membersdf
            print('Clusters and members data loaded from CSV files...')

        # Cross-reference clusters and members with StarHorse
        print('Getting StarHorse data for cross-referencing...')
        self.crossReference(clusters, members)

        cross_referenced_file = 'Data/cross-referenced-results.csv'

        if os.path.exists(cross_referenced_file):
            print("Cross-referencing completed successfully, converting to Parquet...")
            df = pd.read_csv(cross_referenced_file)
            # Sorting by ID
            df.sort_values(by='ID')
            # Dropping duplicates based on 'GaiaDR3' column
            df.drop_duplicates(subset='GaiaDR3', keep='first', inplace=True)
            # Converting to parquet file for efficient storage
            df.to_parquet('Data/cross-referenced-results.parquet', index=False)
        else:
            print("No matches found or failed...")