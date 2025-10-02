#!/usr/bin/env python
# coding: utf-8

# In[6]:


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs


# In[7]:


def mosaic_fits_paper(image_files, weight_files, output_file="final_mosaic.fits"):
    """Create a mosaic of multiple FITS images using inverse-variance weight maps and preserving metadata."""
    
    images = []
    weights = []
    wcs_list = []
    primary_header = None  # Store the full header from the first image

    # Load images and WCS
    for i, file in enumerate(image_files):
        with fits.open(file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header, naxis=2)  # Extract only RA/DEC WCS

            # Copy the header from the first file
            if primary_header is None:
                primary_header = header.copy()  # Preserve all metadata

            # Extract only the 2D intensity map (first frequency & Stokes)
            while data.ndim > 2:
                data = data[0]

            images.append(data)
            wcs_list.append(wcs)

        # Load weight map (assumed to be RMS maps)
        with fits.open(weight_files[i]) as hdul:
            weight_data = hdul[0].data

            # Extract only 2D weight map
            while weight_data.ndim > 2:
                weight_data = weight_data[0]

            # **Fix:** Avoid division by zero or negative values
            weight_data = np.where(weight_data > 0, 1.0 / (weight_data**2), 0)

            weights.append(weight_data)

    # Determine optimal WCS to cover all images
    mosaic_wcs, mosaic_shape = find_optimal_celestial_wcs([(img, wcs) for img, wcs in zip(images, wcs_list)])

    # Initialize arrays for summing images and weights
    mosaic_sum = np.zeros(mosaic_shape)
    weight_sum = np.zeros(mosaic_shape)

    # Reproject and combine all images
    for i, (img, wcs) in enumerate(zip(images, wcs_list)):
        print(f"Reprojecting {image_files[i]}...")
        reprojected, footprint = reproject_interp((img, wcs), mosaic_wcs, shape_out=mosaic_shape)
        weight_reprojected, _ = reproject_interp((weights[i], wcs_list[i]), mosaic_wcs, shape_out=mosaic_shape)

        # Ensure valid_pixels is a 2D mask matching `reprojected`
        valid_pixels = footprint > 0  # Boolean 2D mask

        # Apply the mask properly
        reprojected[~valid_pixels] = np.nan  # Mask invalid areas
        weight_reprojected[~valid_pixels] = 0  # Mask invalid areas in weight

        # Accumulate weighted sum
        mosaic_sum += np.where(valid_pixels, np.nan_to_num(reprojected) * weight_reprojected, 0)
        weight_sum += np.where(valid_pixels, weight_reprojected, 0)

    # Avoid division by zero in regions with no data
    with np.errstate(divide='ignore', invalid='ignore'):
        mosaic_final = np.where(weight_sum > 0, mosaic_sum / weight_sum, np.nan)  # Normalize only valid regions

    # Mask areas with no footprint
    mosaic_final[weight_sum == 0] = np.nan  # Set undefined regions to NaN

    # Update the header with the new WCS but keep all other metadata
    primary_header.update(mosaic_wcs.to_header())  # Only update WCS keys

    # Save the mosaic with the full preserved header
    fits.writeto(output_file, mosaic_final[np.newaxis, :, :], header=primary_header, overwrite=True)
    print(f"Mosaic saved as {output_file} with full metadata.")

# from reproject import reproject_interp, find_optimal_celestial_wcs
# from astropy.io import fits
# from astropy.wcs import WCS
# import numpy as np

def mosaic_fits(image_files, weight_files=None, output_file="final_mosaic.fits"):
    """
    Create a mosaic of multiple FITS images using inverse-variance weight maps (optional).
    If no weights are given, uniform weights (=1) are used.
    Preserves metadata from the first FITS header.
    """
    
    images = []
    weights = []
    wcs_list = []
    primary_header = None  # Store the full header from the first image

    # Load images and WCS
    for i, file in enumerate(image_files):
        with fits.open(file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header, naxis=2)  # Extract only RA/DEC WCS

            # Copy the header from the first file
            if primary_header is None:
                primary_header = header.copy()  # Preserve all metadata

            # Extract only the 2D intensity map (first frequency & Stokes)
            while data.ndim > 2:
                data = data[0]

            images.append(data)
            wcs_list.append(wcs)

        # --- Load weight maps if provided ---
        if weight_files is not None and i < len(weight_files):
            with fits.open(weight_files[i]) as hdul:
                weight_data = hdul[0].data
                while weight_data.ndim > 2:
                    weight_data = weight_data[0]

                # Inverse-variance weights (1 / rms^2), safe against zeros
                weight_data = np.where(weight_data > 0, 1.0 / (weight_data**2), 0)
        else:
            # Use uniform weight = 1 everywhere
            weight_data = np.ones_like(data)

        weights.append(weight_data)

    # Determine optimal WCS to cover all images
    mosaic_wcs, mosaic_shape = find_optimal_celestial_wcs([(img, wcs) for img, wcs in zip(images, wcs_list)])

    # Initialize arrays for summing images and weights
    mosaic_sum = np.zeros(mosaic_shape)
    weight_sum = np.zeros(mosaic_shape)

    # Reproject and combine all images
    for i, (img, wcs) in enumerate(zip(images, wcs_list)):
        print(f"Reprojecting {image_files[i]}...")
        reprojected, footprint = reproject_interp((img, wcs), mosaic_wcs, shape_out=mosaic_shape)
        weight_reprojected, _ = reproject_interp((weights[i], wcs), mosaic_wcs, shape_out=mosaic_shape)

        valid_pixels = footprint > 0

        reprojected[~valid_pixels] = np.nan
        weight_reprojected[~valid_pixels] = 0

        mosaic_sum += np.where(valid_pixels, np.nan_to_num(reprojected) * weight_reprojected, 0)
        weight_sum += np.where(valid_pixels, weight_reprojected, 0)

    # Normalize weighted mosaic
    with np.errstate(divide='ignore', invalid='ignore'):
        mosaic_final = np.where(weight_sum > 0, mosaic_sum / weight_sum, np.nan)

    mosaic_final[weight_sum == 0] = np.nan

    # Update header with new WCS
    primary_header.update(mosaic_wcs.to_header())

    # Save final mosaic
    fits.writeto(output_file, mosaic_final[np.newaxis, :, :], header=primary_header, overwrite=True)
    print(f"Mosaic saved as {output_file} with full metadata (weights={'given' if weight_files else 'uniform=1'}).")

from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from reproject import reproject_interp, reproject_exact, reproject_adaptive
from reproject.mosaicking import find_optimal_celestial_wcs

def reproject_mosaic(image_files, weight_files, output_file="final_mosaic.fits", method="interp"):
    """Create a mosaic of multiple FITS images using inverse-variance weight maps and preserving metadata.
    
    Parameters:
        image_files (list): List of FITS image file paths.
        weight_files (list): List of corresponding weight map file paths.
        output_file (str): Filename for the final mosaic output.
        method (str): Reprojection method ('interp', 'exact', 'adaptive').
    """
    
    images = []
    weights = []
    wcs_list = []
    primary_header = None  # Store metadata from the first image
    
    # Load images and WCS
    for i, file in enumerate(image_files):
        with fits.open(file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header, naxis=2)
            
            if primary_header is None:
                primary_header = header.copy()
            
            while data.ndim > 2:
                data = data[0]
            
            images.append(data)
            wcs_list.append(wcs)
        
        # Load weight map (assumed to be RMS maps)
        with fits.open(weight_files[i]) as hdul:
            weight_data = hdul[0].data
            while weight_data.ndim > 2:
                weight_data = weight_data[0]
            
            weight_data = np.where(weight_data > 0, 1.0 / (weight_data**2), 0)
            weights.append(weight_data)
    
    # Determine optimal WCS
    mosaic_wcs, mosaic_shape = find_optimal_celestial_wcs([(img, wcs) for img, wcs in zip(images, wcs_list)])
    
    mosaic_sum = np.zeros(mosaic_shape)
    weight_sum = np.zeros(mosaic_shape)
    
    # Choose reprojection method
    reprojection_methods = {
        "interp": reproject_interp,
        "exact": reproject_exact,
        "adaptive": reproject_adaptive
    }
    
    if method not in reprojection_methods:
        raise ValueError("Invalid reprojection method. Choose from 'interp', 'exact', or 'adaptive'.")
    
    reproject_fn = reprojection_methods[method]
    
    # Reproject and combine images
    for i, (img, wcs) in enumerate(zip(images, wcs_list)):
        print(f"Reprojecting {image_files[i]} using {method}...")
        reprojected, footprint = reproject_fn((img, wcs), mosaic_wcs, shape_out=mosaic_shape)
        weight_reprojected, _ = reproject_fn((weights[i], wcs_list[i]), mosaic_wcs, shape_out=mosaic_shape)
        
        valid_pixels = footprint > 0
        reprojected[~valid_pixels] = np.nan
        weight_reprojected[~valid_pixels] = 0
        
        mosaic_sum += np.where(valid_pixels, np.nan_to_num(reprojected) * weight_reprojected, 0)
        weight_sum += np.where(valid_pixels, weight_reprojected, 0)
    
    # Normalize and mask invalid regions
    with np.errstate(divide='ignore', invalid='ignore'):
        mosaic_final = np.where(weight_sum > 0, mosaic_sum / weight_sum, np.nan)
    
    mosaic_final[weight_sum == 0] = np.nan
    
    # Update header and save output
    primary_header.update(mosaic_wcs.to_header())
    fits.writeto(output_file, mosaic_final[np.newaxis, :, :], header=primary_header, overwrite=True)
    print(f"Mosaic saved as {output_file} with full metadata.")


from scipy.ndimage import gaussian_filter

def cosine_taper(shape, taper_fraction=0.1):
    """Create a 2D cosine taper mask for image edges."""
    ny, nx = shape
    y = np.hanning(int(ny * taper_fraction))
    x = np.hanning(int(nx * taper_fraction))
    taper = np.outer(y, x)
    
    # Expand to full size by padding with ones
    full_taper = np.ones((ny, nx))
    full_taper[: y.size, : x.size] *= taper
    full_taper[: y.size, -x.size :] *= taper[::-1, :]
    full_taper[-y.size :, : x.size] *= taper[:, ::-1]
    full_taper[-y.size :, -x.size :] *= taper[::-1, ::-1]
    
    return full_taper

def mosaic_fits_smooth(image_files, weight_files, output_file="final_mosaic.fits", smooth_sigma=3, taper_fraction=0.1):
    """Create a mosaic of multiple FITS images using inverse-variance weight maps, smoothing overlaps."""
    images = []
    weights = []
    wcs_list = []
    primary_header = None  # Store the full header from the first image
    
    # Load images and WCS
    for i, file in enumerate(image_files):
        with fits.open(file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header, naxis=2)  # Extract only RA/DEC WCS

            if primary_header is None:
                primary_header = header.copy()  # Preserve metadata
            
            while data.ndim > 2:
                data = data[0]  # Extract 2D intensity map

            # Apply cosine taper
            data *= cosine_taper(data.shape, taper_fraction)
            
            images.append(data)
            wcs_list.append(wcs)
        
        # Load weight map (assumed to be RMS maps)
        with fits.open(weight_files[i]) as hdul:
            weight_data = hdul[0].data
            while weight_data.ndim > 2:
                weight_data = weight_data[0]  # Extract 2D weight map

            # Convert RMS map to inverse-variance weights, avoiding division by zero
            weight_data = np.where(weight_data > 0, 1.0 / (weight_data**2), 0)
            
            # Apply Gaussian smoothing to weight map
            weight_data = gaussian_filter(weight_data, sigma=smooth_sigma)
            
            weights.append(weight_data)
    
    # Determine optimal WCS for the mosaic
    mosaic_wcs, mosaic_shape = find_optimal_celestial_wcs([(img, wcs) for img, wcs in zip(images, wcs_list)])
    
    mosaic_sum = np.zeros(mosaic_shape)
    weight_sum = np.zeros(mosaic_shape)
    
    # Reproject and combine all images
    for i, (img, wcs) in enumerate(zip(images, wcs_list)):
        print(f"Reprojecting {image_files[i]}...")
        reprojected, footprint = reproject_interp((img, wcs), mosaic_wcs, shape_out=mosaic_shape)
        weight_reprojected, _ = reproject_interp((weights[i], wcs_list[i]), mosaic_wcs, shape_out=mosaic_shape)

        valid_pixels = footprint > 0  # Mask invalid areas
        reprojected[~valid_pixels] = np.nan
        weight_reprojected[~valid_pixels] = 0

        mosaic_sum += np.where(valid_pixels, np.nan_to_num(reprojected) * weight_reprojected, 0)
        weight_sum += np.where(valid_pixels, weight_reprojected, 0)
    
    # Normalize by weights
    with np.errstate(divide='ignore', invalid='ignore'):
        mosaic_final = np.where(weight_sum > 0, mosaic_sum / weight_sum, np.nan)
    
    # Mask undefined regions
    mosaic_final[weight_sum == 0] = np.nan
    
    # Update header with new WCS
    primary_header.update(mosaic_wcs.to_header())
    
    # Save the mosaic with full metadata
    fits.writeto(output_file, mosaic_final[np.newaxis, :, :], header=primary_header, overwrite=True)
    print(f"Mosaic saved as {output_file} with full metadata.")



"""
# Example usage
input_images = ["AAAA.SP2B.PBCOR.FITS", "BBBB.SP2B.PBCOR.FITS", "CCCC.SP2B.PBCOR.FITS", "DDDD.SP2B.PBCOR.FITS", 
                "EEEE.SP2B.PBCOR.FITS", "FFFF.SP2B.PBCOR.FITS", "GGGG.SP2B.PBCOR.FITS"]
weights = ["AAAA.SP2B.PBCOR.pybdsf_rms.fits", "BBBB.SP2B.PBCOR.pybdsf_rms.fits", "CCCC.SP2B.PBCOR.pybdsf_rms.fits", "DDDD.SP2B.PBCOR.pybdsf_rms.fits", 
                "EEEE.SP2B.PBCOR.pybdsf_rms.fits", "FFFF.SP2B.PBCOR.pybdsf_rms.fits", "GGGG.SP2B.PBCOR.pybdsf_rms.fits"]
# input_images = ["AAAA.SP2B.PBCOR.FITS", "BBBB.SP2B.PBCOR.FITS", "CCCC.SP2B.PBCOR.FITS", "DDDD.SP2B.PBCOR.FITS", 
#                 "EEEE.SP2B.PBCOR.FITS", "FFFF.SP2B.PBCOR.FITS", "GGGG.SP2B.PBCOR.FITS"]
mosaic_fits(input_images, weights)
"""

# In[ ]:




