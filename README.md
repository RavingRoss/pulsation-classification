# Pulsation Classification
Includes my findings and past/current research on pulsators which include roAps, Blue Stragglers, delta scuti, gamma dors, and maybe more in the future!

---

# Acknowledgments 

I would like to acknowledge the team who developed the lightkurve package, as this project would not be where it is without them! I am excited to see what the future holds for that team of devs and for the future of Asteroseismology!

---

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

![Example of Final Plot](https://github.com/RavingRoss/pulsation-classification/blob/main/Blue%20Stragglers/Data/TESS%20Data/Graphs/Final%20Plot%20of%20TIC%20157567602.png)
