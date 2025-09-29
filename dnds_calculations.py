import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import convolve
from astropy.io import fits
from astropy.wcs import WCS
import bdsf
import joblib
from joblib import Parallel, delayed

####### SECTION 1: GENERATE MOCK SOURCES #######

# Power-law function for fitting
def power_law(S, A, alpha):
    """Power-law function: A * S^alpha"""
    return A * S**alpha

# Generate mock sources following dN/dS ∝ S^-1.6
def generate_mock_sources(num_sources, flux_min, flux_max, extended_fraction, seed=None):
    """Generate mock sources with power-law distributed flux densities."""
    if seed is not None:
        np.random.seed(seed)
    
    alpha = -1.6  # For dN/dS ∝ S^-1.6
    
    # Generate uniform random numbers
    u = np.random.uniform(0, 1, num_sources)

    # Sample from the power-law distribution using inverse transform sampling
    flux_densities = flux_min * (((flux_max / flux_min) ** (alpha + 1) - 1) * u + 1) ** (1 / (alpha + 1))

    # Mark some sources as extended
    is_extended = np.random.rand(num_sources) < extended_fraction

    return flux_densities, is_extended

# Compute differential source count dN/dS
def compute_dn_ds(flux_densities, flux_bins):
    """Compute differential source count dN/dS."""
    hist, bin_edges = np.histogram(flux_densities, bins=flux_bins)
    
    # Compute bin centers and widths
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = np.diff(bin_edges)

    # Compute dN/dS (normalize counts by bin width)
    dn_ds = hist / bin_widths

    return bin_centers, dn_ds

# Fit a power law to dN/dS
def fit_power_law(flux_densities, flux_bins):
    """Fit a power law to dN/dS."""
    flux_vals, dn_ds_vals = compute_dn_ds(flux_densities, flux_bins)

    # Filter out invalid values
    valid = (dn_ds_vals > 0) & (flux_vals > 0)
    log_flux_vals = np.log10(flux_vals[valid])
    log_dn_ds_vals = np.log10(dn_ds_vals[valid])

    # Perform curve fitting on log-log scale
    popt, _ = curve_fit(lambda logS, logA, alpha: logA + alpha * logS, log_flux_vals, log_dn_ds_vals, p0=[3, -1.6])
    
    # Extract parameters
    logA_fit, alpha_fit = popt
    A_fit = 10**logA_fit  # Convert log(A) back to A

    return A_fit, alpha_fit

# Plot dN/dS with power-law fit
def plot_dn_ds(flux_densities, flux_bins):
    """Plot the differential source count dN/dS with power-law fit."""
    flux_vals, dn_ds_vals = compute_dn_ds(flux_densities, flux_bins)
    
    # Fit power law
    A_fit, alpha_fit = fit_power_law(flux_densities, flux_bins)

    # Plot results
    plt.figure(figsize=(8, 5))
    plt.loglog(flux_vals, dn_ds_vals, marker='o', linestyle='none', label='Generated Sources')
    plt.loglog(flux_vals, power_law(flux_vals, A_fit, alpha_fit), linestyle='--', label=f'Fit: $S^{{{alpha_fit:.2f}}}$')

    plt.xlabel("Flux Density (Jy)")
    plt.ylabel(r"$dN/dS$")
    plt.title("Differential Source Counts")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

    print(f"Fitted power-law index: {alpha_fit:.2f} (Expected: -1.6)")



####### SECTION 4: RUN PYBDSF #######

def run_pybdsf(image_filename):
    """Run PyBDSF source finder on the image."""
    img = bdsf.process_image(image_filename, psf_vary_do=True, adaptive_rms_box = True, rms_box=(150, 30), rms_box_bright=(50, 10), clobber=True, quiet=True)
    img.write_catalog(format='fits', catalog_type='srl', clobber=True)

def run_pybdsf_residual(image_filename):
    """Run PyBDSF source finder on the image."""
    img = bdsf.process_image(image_filename, psf_vary_do=True, adaptive_rms_box = True, rms_box=(150, 30), rms_box_bright=(50, 10), clobber=True, quiet=True)
    img.export_image(img_type='gaus_resid', clobber=True)


def run_pybdsf_neg(image_filename):
    """Run PyBDSF source finder on the image."""
    img = bdsf.process_image(image_filename, adaptive_rms_box = True, rms_box=(150, 30), rms_box_bright=(50, 10), clobber=True, quiet=True)
    img.write_catalog(format='fits', catalog_type='srl', clobber=True)

def read_pybdsf_catalog(catalog_fits):
    """Read PyBDSF FITS catalog and return total flux density."""
    with fits.open(catalog_fits) as hdul:
        return np.array(hdul[1].data['Total_flux'])

# def run_pybdsf(image_filename):
#     """Run PyBDSF source finder on the image and return total flux densities."""
#     img = bdsf.process_image(image_filename, psf_vary_do=True, adaptive_rms_box=True, 
#                              rms_box=(150, 30), rms_box_bright=(50, 10), clobber=True, quiet=True)
    
#     # Extract total flux density directly from the catalog
#     total_flux = np.array(img.srl_catalog['Total_flux'])

#     return total_flux

# def run_pybdsf(image_filename):
#     """Run PyBDSF source finder on the image and return total flux densities."""
#     img = bdsf.process_image(image_filename, psf_vary_do=False, adaptive_rms_box=True, 
#                              rms_box=(150, 30), rms_box_bright=(50, 10), clobber=True, quiet=True)
    
#     # Check if sources were detected
#     if hasattr(img, 'srl'):
#         total_flux = np.array(img.srl['Total_flux'])
#     else:
#         total_flux = np.array([])  # Return an empty array if no sources are found

#     return total_flux


####### SECTION 3: INJECT SOURCES INTO IMAGE #######


def create_psf(header, pixel_scale, is_extended):
    """Generate a CLEAN beam PSF kernel using BMAJ, BMIN, and BPA from the FITS header."""
    
    # Check if required beam parameters exist in the header
    if not all(key in header for key in ('BMAJ', 'BMIN', 'BPA')):
        raise ValueError("FITS header missing BMAJ, BMIN, or BPA.")

    # Convert beam major/minor axes from degrees to pixels
    bmaj = header['BMAJ'] * 3600 / pixel_scale  
    bmin = header['BMIN'] * 3600 / pixel_scale  
    bpa = np.radians(header['BPA'])  # Convert position angle to radians

    # Scale up beam size for extended sources
    if is_extended:
        # scale_factor = 2.0  # Modify this factor as needed
        scale_factor = 1.2  # Modify this factor as needed
        bmaj *= scale_factor
        bmin *= scale_factor

    # Define PSF kernel size (ensure it's odd and sufficiently large)
    kernel_size = max(int(6 * max(bmaj, bmin)), 5)
    kernel_size += 1 - kernel_size % 2  # Force odd size

    # Generate coordinate grid centered at zero
    y, x = np.meshgrid(
        np.arange(kernel_size) - kernel_size // 2,
        np.arange(kernel_size) - kernel_size // 2,
        indexing='ij'
    )

    # Apply correct astronomical convention for rotation
    x_rot = x * np.cos(-bpa) + y * np.sin(-bpa)
    y_rot = -x * np.sin(-bpa) + y * np.cos(-bpa)

    # Define elliptical Gaussian PSF
    psf_kernel = np.exp(-0.5 * ((x_rot / bmaj) ** 2 + (y_rot / bmin) ** 2))

    # Normalize the PSF to ensure flux conservation
    psf_kernel /= np.sum(psf_kernel)
    beam_area = (np.pi * bmaj * bmin) / (4 * np.log(2))
    
    return beam_area, psf_kernel


def inject_single_source(x, y, flux, beam_area, psf_kernel, shape):
    """Inject a single source convolved with the PSF."""
    kernel_size = psf_kernel.shape[0]
    half_size = kernel_size // 2

    # Skip if out of bounds
    if not (0 <= x < shape[1] and 0 <= y < shape[0]):
        return np.zeros(shape)

    # Define patch region
    x_min, x_max = max(0, x - half_size), min(shape[1], x + half_size + 1)
    y_min, y_max = max(0, y - half_size), min(shape[0], y + half_size + 1)

    patch_x_min, patch_x_max = half_size - (x - x_min), half_size + (x_max - x)
    patch_y_min, patch_y_max = half_size - (y - y_min), half_size + (y_max - y)

    # Create a point source at the center
    patch = np.zeros((kernel_size, kernel_size))
    patch[half_size, half_size] = flux 

    # Convolve with the PSF
    convolved_patch = convolve(patch, psf_kernel) * beam_area # Convert Jy/pixel to Jy/beam

    # Inject only within valid region
    injected_patch = np.zeros(shape)
    injected_patch[y_min:y_max, x_min:x_max] = convolved_patch[patch_y_min:patch_y_max, patch_x_min:patch_x_max]

    return injected_patch


def inject_sources(image_filename, flux_densities, is_extended_list, output_filename, n_jobs=20, batch_size=20):
    with fits.open(image_filename, mode='readonly') as hdu:
        header = hdu[0].header
        wcs = WCS(header, naxis=2)
        original_data = hdu[0].data.astype(np.float64)

    if original_data.ndim == 4:
        base_data = original_data[0, 0, :, :]
    elif original_data.ndim == 3:
        base_data = original_data[0, :, :]
    elif original_data.ndim == 2:
        base_data = original_data.copy()
    else:
        raise ValueError(f"Unexpected FITS data dimensions: {original_data.shape}")

    pixel_scale_arcsec = np.abs(header['CDELT1']) * 3600

    # Identify valid (non-NaN) pixel locations
    valid_pixels = np.argwhere(~np.isnan(base_data))
    valid_ra, valid_dec = wcs.all_pix2world(valid_pixels[:, 1], valid_pixels[:, 0], 0)

    # Sample RA, Dec positions only from the valid region
    num_sources = len(flux_densities)
    selected_indices = np.random.choice(len(valid_ra), num_sources, replace=False)
    ra_values, dec_values = valid_ra[selected_indices], valid_dec[selected_indices]

    x_pix, y_pix = wcs.all_world2pix(ra_values, dec_values, 0)

    new_data = base_data.copy()
    for i in range(0, num_sources, batch_size):
        batch_x = x_pix[i:i + batch_size]
        batch_y = y_pix[i:i + batch_size]
        batch_flux = flux_densities[i:i + batch_size]
        batch_is_extended = is_extended_list[i:i + batch_size]

        injected_patches = Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(inject_single_source)(
                int(round(batch_x[j])), int(round(batch_y[j])), batch_flux[j],
                *create_psf(header, pixel_scale_arcsec, batch_is_extended[j]),
                base_data.shape
            )
            for j in range(len(batch_x))
        )

        for patch in injected_patches:
            new_data += patch

    hdu_modified = fits.PrimaryHDU(new_data, header=header)
    hdu_modified.writeto(output_filename, overwrite=True)

    return ra_values, dec_values



def invert_fits_image(input_fits, output_fits):
    """Invert only the non-NaN pixels in the FITS image."""
    with fits.open(input_fits) as hdul:
        data = hdul[0].data

        # Invert only non-NaN values
        mask = np.isnan(data)  # Identify NaN values
        data = -data  # Invert all values
        data[mask] = np.nan  # Restore NaN values

        # Save the new FITS file
        hdul[0].data = data
        hdul.writeto(output_fits, overwrite=True)

    print(f"Inverted FITS image saved to {output_fits}")


def compute_fdr(original_fluxes, negative_fluxes, flux_bins):
    """Estimate the false detection rate per flux bin."""
    # flux_bins = np.logspace(np.log10(min(original_fluxes)), np.log10(max(original_fluxes)), bins)
    # flux_bins = np.logspace(np.log10(0.01), np.log10(20), bins)
    # bin_centers = (flux_bins[:-1] + flux_bins[1:]) / 2
    
    pos_counts, _ = np.histogram(original_fluxes, bins=flux_bins)
    neg_counts, _ = np.histogram(np.abs(negative_fluxes), bins=flux_bins)  # Use abs() for negative fluxes
    
    fdr = np.divide(neg_counts, pos_counts + neg_counts, out=np.zeros_like(neg_counts, dtype=float), where=(pos_counts + neg_counts) > 0)
    # fdr = np.divide(pos_counts-neg_counts, pos_counts, out=np.zeros_like(neg_counts, dtype=float), where=(pos_counts + neg_counts) > 0)
    # hale et al. 19, Sagar et al 2025. 
    return pos_counts, neg_counts


#############
def compute_histogram(flux_densities, bins):
    flux_min, flux_max = min(flux_densities), max(flux_densities) # (0.02, 21.58725416794234)
    flux_min = 0.02 if (flux_min<0.02) else flux_min
    
    flux_bins = np.logspace(np.log10(flux_min), np.log10(flux_max), bins)
    bin_centers = (flux_bins[:-1] + flux_bins[1:]) / 2
    counts, bin_edges = np.histogram(flux_densities, bins=flux_bins)
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    return flux_bins, bin_centers, counts, bin_edges, bin_widths



def visibility_function(flux_bins, sigma=5., rms_fits_file="final_mosaic.pybdsf_rms.fits"):
    """
    Computes the visibility function correction V(S) by determining the survey area 
    over which sources at different flux levels can be detected.

    Parameters:
    - flux_bins (numpy array): Flux density bins (in Jy).
    - sigma (float): Detection threshold in units of sigma (default: 5σ).
    - rms_fits_file (str): Path to the RMS map FITS file.

    Returns:
    - total_valid_area (float): The total survey area (in square degrees).
    - effective_image_area (numpy array): Effective area corresponding to each flux bin.
    """

    # Load RMS map
    with fits.open(rms_fits_file) as hdul_rms:
        rms_data = hdul_rms[0].data[0, :, :]
        hdr = hdul_rms[0].header

    # Compute pixel area in square degrees
    CDELT1 = abs(hdr.get("CDELT1", 1))
    CDELT2 = abs(hdr.get("CDELT2", CDELT1))
    pixel_area_deg2 = CDELT1 * CDELT2  

    # Mask NaN values
    valid_mask = np.isfinite(rms_data)
    total_valid_area = np.sum(valid_mask) * pixel_area_deg2  # Total survey area

    # Flatten RMS values in valid regions
    rms_flat = rms_data[valid_mask].flatten()

    # Compute detection thresholds per pixel
    detection_thresholds = rms_flat * sigma  

    # # Effective survey area for each flux bin
    # effective_image_area = np.array([np.sum(flux_bins[i] >= detection_thresholds) * pixel_area_deg2 
    #                                  for i in range(len(flux_bins))])

    effective_image_area = np.array([
        np.sum((flux_bins[i] >= detection_thresholds) & (flux_bins[i] < flux_bins[i+1])) * pixel_area_deg2 
        for i in range(len(flux_bins) - 1)])

    deg2_to_sr = (np.pi / 180) ** 2  # Conversion factor from deg² to sr
    effective_image_area_sr = effective_image_area * deg2_to_sr  # Convert to steradians
    total_valid_area_sr = total_valid_area * deg2_to_sr  # Convert to steradians

    return total_valid_area_sr, effective_image_area_sr


# --- Functions ---

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
import shutil
import os
import warnings

# import dnsn_calculations as dnds


# Suppress warnings
# warnings.filterwarnings("ignore", category=UserWarning)


def generate_synthetic_catalog(input_catalog, output_catalog, flux_densities, is_extended, valid_ra, valid_dec, seed=None):
    if seed is not None:
        np.random.seed(seed)
    with fits.open(input_catalog) as hdul:
        data = hdul[1].data

    # Beam parameters (arcsec)
    beam_maj = data['Maj'] * 3600
    beam_min = data['Min'] * 3600
    mean_maj, mean_min = beam_maj.mean(), beam_min.mean()

    num = len(flux_densities)
    a_new = np.where(is_extended, 2 * mean_maj, mean_maj)
    b_new = np.where(is_extended, 2 * mean_min, mean_min)
    new_peak = flux_densities * (mean_maj*mean_min) / (a_new*b_new)

    idx = np.random.choice(len(valid_ra), num, replace=False)
    new_ra, new_dec = valid_ra[idx], valid_dec[idx]
    pa = np.random.uniform(14, 25, num)

    df = pd.DataFrame({
        'ra': new_ra,
        'dec': new_dec,
        'peak_flux': new_peak,
        'a': a_new,
        'b': b_new,
        'pa': pa
    })
    Table.from_pandas(df).write(output_catalog, format='votable', overwrite=True)
    print(f"Generated {num} synthetic sources")
    return flux_densities, new_ra, new_dec


def match_recovered_sources(inj_ra, inj_dec, rec_ra, rec_dec, max_sep=40.0):
    inj = SkyCoord(inj_ra*u.deg, inj_dec*u.deg)
    rec = SkyCoord(rec_ra*u.deg, rec_dec*u.deg)
    idx, sep2d, _ = match_coordinates_sky(inj, rec)
    mask = sep2d.arcsec < max_sep
    return mask, idx

from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree
import numpy as np

def match_recovered_sources(inj_ra, inj_dec, rec_ra, rec_dec, max_sep=30.0):
    """
    Match injected sources to recovered sources within a given separation threshold.

    Parameters:
    - inj_ra, inj_dec : arrays of RA and Dec for injected sources [degrees]
    - rec_ra, rec_dec : arrays of RA and Dec for recovered sources [degrees]
    - max_sep : matching radius [arcsec]

    Returns:
    - inj_unique : indices of matched injected sources
    - rec_unique : corresponding indices of matched recovered sources (1-to-1)
    """

    inj = SkyCoord(inj_ra, inj_dec, unit='deg')
    rec = SkyCoord(rec_ra, rec_dec, unit='deg')

    # Convert to 3D Cartesian coordinates (unit sphere) for efficient matching
    inj_cart = np.vstack(inj.cartesian.xyz).T
    rec_cart = np.vstack(rec.cartesian.xyz).T

    # KD-tree for fast spatial matching
    tree = cKDTree(rec_cart)
    max_rad = (max_sep * u.arcsec).to(u.rad).value  # convert to radians
    dist, idx = tree.query(inj_cart, distance_upper_bound=max_rad)

    # Identify valid matches (i.e., within max_sep)
    matched = dist != np.inf
    inj_matched = np.where(matched)[0]
    rec_matched = idx[matched]

    # Remove duplicate recovered matches (enforce 1-to-1 mapping)
    # Keep only the first instance of each recovered index
    _, unique_idx = np.unique(rec_matched, return_index=True)
    inj_unique = inj_matched[unique_idx]
    rec_unique = rec_matched[unique_idx]

    return inj_unique, rec_unique



