# Radio Mosaic Analysis Codes

![Python](https://img.shields.io/badge/python-3.8%2B-blue) ![License](https://img.shields.io/badge/license-MIT-green)

This repository contains a set of codes for creating mosaics from radio images. The current implementation uses **uGMRT Band 2 images** centered at a single frequency of **147.4 MHz**. Detailed information about the data and methodology can be found in **Elahi et al. (2025), MNRAS (submitted)**.

---

## Table of Contents
- [Dependencies](#dependencies)
- [Capabilities](#capabilities)
- [Data Availability](#data-availability)
- [Usage Instructions](#usage-instructions)
- [Contact](#contact)

---

## Dependencies

The following packages are required:

- **Astropy** – used extensively for FITS handling, WCS, and other astronomical calculations.  
- **PyBDSF** – for background RMS estimation, source finding, and catalogue generation.  
- **CASA** – the `imsmooth` task is used for PSF matching before mosaicing.  
- **Aegean** – specifically the `AeRes` tool, used for adding sources to the residual image for completeness corrections in source counts.  

---

## Capabilities

These codes allow you to:

- Create mosaics from PSF-matched images.  
- Generate source catalogues.  
- Classify sources as point-like or extended.  
- Compute source counts, including necessary corrections for **false detection rate, completeness, and visibility area**.  

---

## Data Availability

- Required images can be requested from the developer or collaborators.  
- **TGSS, GLEAM, and GLEAM-X** data are publicly available via the **[Vizier Astronomical Database](https://vizier.u-strasbg.fr/)**.  
- The **uGMRT catalogue** is included in this repository and will soon be available on Vizier.  

---

## Usage Instructions

1. **Start with the FITS images**: `*.SP2B.PBCOR.FITS`.  

2. **PSF Matching**  
   Run CASA with the provided script:  
   ```bash
   casa -c imsmooth.py
   ```

3. **Mosaicing**  
   Open `mosaicing.ipynb` to:  
   - Read PSF-matched images and generate RMS maps using PyBDSF (used as weights).  
   - Create the mosaic using the PSF-matched images and their weights.  
   - Generate the RMS map of the mosaic.  
   - Compute the median RMS of the mosaic (used in the paper abstract).  

4. **Plotting**  
   `plots.ipynb` generates various plots, including:  
   - Observation strategy  
   - Multiple pointing centres (PCs)  
   - Mosaic overview  
   - Triple panel plots  
   - RMS maps  
   - Completeness curves  

5. **Source Classification**  
   `pointextended.ipynb` classifies sources as point-like or extended.  

6. **Crossmatching and Astrometry**  
   `crossmatch.ipynb` performs astrometry checks, flux comparison, and calculates NVSS–uGMRT spectral indices.  

7. **Catalogue Generation**  
   `catalogue.ipynb` prints a sample catalogue of detected sources.  

8. **Bright Source Analysis**  
   `crossmatch-bright-spind.ipynb` identifies bright sources and extracts their information from multiple catalogues.  

9. **Image Summary**  
   `image_summary.ipynb` provides a summary of the mosaic and associated images.

10. **Soure Counts**
    `sourecounts_corrections.ipynb` makes the corrections for **false detection rate, completeness, and visibility area**

---

## Contact

For questions, requests for data, or collaboration inquiries, please contact: **Asif Elahi**.
