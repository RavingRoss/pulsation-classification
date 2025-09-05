# Pulsation Classification Repository 
Includes my findings and past/current research on pulsators which include Blue Stragglers and roAps and will eventually include delta scuti, gamma dors, and more in the future!

---

# Acknowledgments 

I would like to acknowledge the team who developed the lightkurve package, as this project would not be where it is without them! I am excited to see what the future holds for that team of devs and for the future of Asteroseismology!

---

# Dependencies

Easiest way to install the packages is in a virtual environment, which can be made the following way, ```python -m venv .venv```. Then activate the environment, ```source .venv/bin/activate```; it can be deactivated by typing ```deactivate```. The necessary packages can be found in the ```requirements.txt``` file; to install these, run ```pip install -r requirements.txt``` inside the environment. If you would like to generate this file as well, you can do this by running ```pip freeze > requirements.txt``` in the environment. 

The repository uses 'pre-commit' and 'black' for organization and double checking the commit. After editing the code and wanting to commit your changes, 'pre-commit' will run the checks then fail, as it organized the code accordingly. To commit, simply re-add the changes, ``` git add .```, then commit again, ``` git commit -m "message here"```, then push the chagnes, ```git push```. If there is a commit error, the message will show and need to be fixed before committing and pushing the changes. 

# 1. Blue Straggler Candidates

This project identifies and filters potential BSS candidates from stellar cluster data using color-magnitude diagrams (CMDs) and their Power Spectrums (PS). It uses membership probabilities and photometric criteria to select candidates, visualize them, using the lightkurve package, and save the results. 

---

## File Structure

- Blue Stragglers/
  - .py files to run
  - Notes
  - Constraints for pip
  - Data/
    - Candidates/
      - CMD of cluster w/ selected region
      - .csv's of candidates
    - TESS Data/
      - Clusters/
        - .csv's of candidates in their respective clusters
      - Graphs/
        - Final plots of candidates 
      - Light Curves/
        - .csv's of flux and BKJD (time)
      - Power Spectrums/
        - .csv's of amplitude and frequencies
    - Input catalogs (.csv and .parquet)
  - Presentation/
    - Slides
    - Paper
    - Some figures
    
---

## Features

- Filters cluster members based on `Gmag` and `BP-RP` color-magnitude selection and their Power Spectrums.
- Uses membership probability to colormap cluster members.
- Highlights individual candidate stars on top of the full cluster along with a theoretical isochrones of the cluster's ~age.
- Automatically skips or removes results if no candidates match.
- Supports matplotlib subplots with selective colorbars using lightkurve package.
- End result showcases the CMD, lightcurve, and PS of the selected candidate.

---

## Example Plot

> CMD with color-mapped membership probability and a highlighted candidate, including its lightcurve and PS:

![Example of Final Plot](https://github.com/RavingRoss/pulsation-classification/blob/main/Blue%20Stragglers/Data/TESS%20Data/Graphs/High_End/Final%20Plot%20of%20TIC%20157567602.png)

---

# 2. roAp Candidates

This project is focused on building up on the current catalogs of roAp stars. Starting off with running tests on known roAps, and eventually searching for new candidates using TESS. 

---

## File Structure

- roAp/
  - main roAp.py file (not up-to-date)
  - test.py file for creating test plots
  - TIC test plots

--- 

## Features

- Selects stars only with high frequency's.
- Runs selected stars through extra checks:
  - Skew check will classify high/low priority.
  - Manually check test plots.
- Will organize priority once script is stopped.
- Exports data using .parquet files to minimize storage.

---

## Example Plot 

> Includes: lower-end PSD, normalized LC, high-end PSD, and normalized phase-folded LC.

![Example of Test Plot](https://github.com/RavingRoss/pulsation-classification/blob/46f3d8bbd85ed0f864dffd0dd3100ded25569515/roAp/TIC%20310817678.png)
