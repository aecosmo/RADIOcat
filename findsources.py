#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

# Load GLEAM catalog
gleam_catalog = "GLEAM_EGC_v2.fits.gz"
tgss_catalog = "TGSSADR1_7sigma_catalog.fits"
pybdsf_catalog= "final_mosaic.pybdsf.srl.fits"
nvss_catalog= "nvss.dat.gz.fits"
nvss_catalog= "VIII_65_nvss.dat.gz.fits.gz"


def find_sources_old(catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='gleam'):
    if name == 'gleam':
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAJ2000"]   # RA in degrees
            gleam_dec = hdul[1].data["DEJ2000"]  # DEC in degrees
            gleam_flux = hdul[1].data["int_flux_wide"]  # Wideband peak flux density
            gleam_flux = hdul[1].data["int_flux_151"]  # Wideband peak flux density
            gleam_flux_err = hdul[1].data["err_int_flux_wide"]  # Flux error (check column name)
            gleam_flux_err = hdul[1].data["err_int_flux_151"]  # Flux error (check column name)
            # print(gleam_ra, gleam_dec)
    elif name == 'gleamx1':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAdeg"]   # RA in degrees
            gleam_dec = hdul[1].data["DEdeg"]  # DEC in degrees
            gleam_flux = hdul[1].data["Fintwide"]  # Total flux density
            gleam_flux = hdul[1].data["Fint151"]  # Total flux density
            gleam_flux_err = hdul[1].data["e_Fintwide"]  # Flux error (check column name)
            gleam_flux_err = hdul[1].data["e_Fint151"]  # Flux error (check column name)
    elif name == 'tgss':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RA"]   # RA in degrees
            gleam_dec = hdul[1].data["DEC"]  # DEC in degrees
            gleam_flux = hdul[1].data["Total_flux"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Total_flux"]  # Flux error (check column name)


            

    elif name == 'nvss':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra1 = hdul[1].data["RAJ2000"]   # RA in HH MM SS
            gleam_dec1 = hdul[1].data["DEJ2000"]  # DEC in DD MM SS
            gleam_flux = hdul[1].data["S1.4"]  # Total flux density
            gleam_flux_err = hdul[1].data["e_S1.4"]  # Flux error (check column name)
            
            # Convert to SkyCoord (automatically converts to degrees)
            nvss_coords = SkyCoord(gleam_ra1, gleam_dec1, unit=(u.hourangle, u.deg), frame='icrs')
            
            # Convert to degrees
            gleam_ra = nvss_coords.ra.deg  # RA in degrees
            gleam_dec = nvss_coords.dec.deg  # Dec in degrees


            
            # print(gleam_ra, gleam_dec)

    else:  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RA"]   # RA in degrees
            gleam_dec = hdul[1].data["DEC"]  # DEC in degrees
            gleam_flux = hdul[1].data["Total_flux"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Total_flux"]  # Flux error (check column name)

    # Apply flux threshold
    bright_sources = gleam_flux > flux_thresh  # Boolean mask for bright sources

    # Filter sources
    gleam_ra = gleam_ra[bright_sources]
    gleam_dec = gleam_dec[bright_sources]
    gleam_flux = gleam_flux[bright_sources]
    gleam_flux_err = gleam_flux_err[bright_sources]  # Also filter flux errors

    # Convert to SkyCoord
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)

    # Convert input target position to SkyCoord
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    # Find sources within search radius
    sep = target_coord.separation(gleam_coords)
    matched = sep < search_radius

    # Extract matched sources
    gleam_ra = gleam_ra[matched]
    gleam_dec = gleam_dec[matched]
    gleam_flux = gleam_flux[matched]
    gleam_flux_err = gleam_flux_err[matched]

    return gleam_ra, gleam_dec, gleam_flux, gleam_flux_err

def find_sources(catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='gleam'):
    if name == 'gleam':
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAJ2000"]   # RA in degrees
            gleam_dec = hdul[1].data["DEJ2000"]  # DEC in degrees
            gleam_flux = hdul[1].data["int_flux_wide"]  # Wideband peak flux density
            gleam_flux_err = hdul[1].data["err_int_flux_wide"]  # Flux error (check column name)
            # print(gleam_ra, gleam_dec)
    elif name == 'gleamx1':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAdeg"]   # RA in degrees
            gleam_dec = hdul[1].data["DEdeg"]  # DEC in degrees
            gleam_flux = hdul[1].data["Fintwide"]  # Total flux density
            # gleam_flux = hdul[1].data["Fint151"]  # Total flux density
            gleam_flux_err = hdul[1].data["e_Fintwide"]  # Flux error (check column name)
            # gleam_flux_err = hdul[1].data["e_Fint151"]  # Flux error (check column name)
    elif name == 'tgss':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAdeg"]   # RA in degrees
            gleam_dec = hdul[1].data["DEdeg"]  # DEC in degrees
            gleam_flux = hdul[1].data["Stotal"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Stotal"]  # Flux error (check column name)

    elif name == 'nvss':
        with fits.open(catalog) as hdul:
            data = hdul[1].data
    
            # Combine RA and DEC parts into SkyCoord
            ra_str = [f"{ra_h}h{ra_m}m{ra_s}s" for ra_h, ra_m, ra_s in zip(data['RAh'], data['RAm'], data['RAs'])]
            sign = np.where(data['DE-'] == '-', -1, 1)  # Get the sign from 'DE-'
            dec_deg = sign * (data['DEd'] + data['DEm']/60 + data['DEs']/3600)
            
            coords = SkyCoord(ra=ra_str, dec=dec_deg*u.deg, frame='icrs')
    
            gleam_ra = coords.ra.deg
            gleam_dec = coords.dec.deg
    
            # Flux and uncertainty (in Jy)
            gleam_flux = data["S1.4"] / 1000.0  # Convert mJy to Jy
            gleam_flux_err = data["e_S1.4"] / 1000.0  # Convert mJy to Jy
            
            # print(gleam_ra, gleam_dec)

    else:  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RA"]   # RA in degrees
            gleam_dec = hdul[1].data["DEC"]  # DEC in degrees
            gleam_flux = hdul[1].data["Total_flux"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Total_flux"]  # Flux error (check column name)

    # Apply flux threshold
    bright_sources = gleam_flux > flux_thresh  # Boolean mask for bright sources

    # Filter sources
    gleam_ra = gleam_ra[bright_sources]
    gleam_dec = gleam_dec[bright_sources]
    gleam_flux = gleam_flux[bright_sources]
    gleam_flux_err = gleam_flux_err[bright_sources]  # Also filter flux errors

    # Convert to SkyCoord
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)

    # Convert input target position to SkyCoord
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    # Find sources within search radius
    sep = target_coord.separation(gleam_coords)
    matched = sep < search_radius

    # Extract matched sources
    gleam_ra = gleam_ra[matched]
    gleam_dec = gleam_dec[matched]
    gleam_flux = gleam_flux[matched]
    gleam_flux_err = gleam_flux_err[matched]

    return gleam_ra, gleam_dec, gleam_flux, gleam_flux_err

    

def find_sources_peak_int(flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='gleam'):
    if name == 'gleam':
        with fits.open(gleam_catalog) as hdul:
            ra = hdul[1].data["RAJ2000"]   # RA in degrees
            dec = hdul[1].data["DEJ2000"]  # DEC in degrees
            peak_flux = hdul[1].data["peak_flux_wide"]  # Total flux density
            peak_flux_err = hdul[1].data["err_peak_flux_wide"]  # Flux error (check column name)
            int_flux = hdul[1].data["int_flux_wide"]  # Total flux density
            int_flux_err = hdul[1].data["err_int_flux_wide"]  # Total flux density
            local_rms = hdul[1].data["local_rms_wide"]  # Flux error (check column name)
            
    elif name == 'gleamx1':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            ra = hdul[1].data["RAdeg"]   # RA in degrees
            dec = hdul[1].data["DEdeg"]  # DEC in degrees
            peak_flux = hdul[1].data["Fpwide"]  # Total flux density
            peak_flux_err = hdul[1].data["e_Fpwide"]  # Flux error (check column name)
            int_flux = hdul[1].data["Fintwide"]  # Total flux density
            int_flux_err = hdul[1].data["e_Fintwide"]  # Flux error (check column name)
            local_rms = hdul[1].data["lrmswide"]  # Flux error (check column name)

    elif name == 'tgss':  # TGSS or other catalogs
        with fits.open(tgss_catalog) as hdul:
            ra = hdul[1].data["RA"]   # RA in degrees
            dec = hdul[1].data["DEC"]  # DEC in degrees
            peak_flux = hdul[1].data["Peak_flux"]  # Total flux density
            peak_flux_err = hdul[1].data["E_Peak_flux"]  # Flux error (check column name)
            int_flux = hdul[1].data["Total_flux"]  # Total flux density
            int_flux_err = hdul[1].data["E_Total_flux"]  # Total flux density
            local_rms = hdul[1].data["RMS_noise"]  # Flux error (check column name)
            
    else:  # pybdsf
        with fits.open(pybdsf_catalog) as hdul:
            ra = hdul[1].data["RA"]   # RA in degrees
            dec = hdul[1].data["DEC"]  # DEC in degrees
            peak_flux = hdul[1].data["Peak_flux"]  # Total flux density
            peak_flux_err = hdul[1].data["E_Peak_flux"]  # Flux error (check column name)
            int_flux = hdul[1].data["Total_flux"]  # Total flux density
            int_flux_err = hdul[1].data["E_Total_flux"]  # Total flux density
            local_rms = hdul[1].data["Isl_rms"]  # Flux error (check column name)
    
    # Apply flux threshold upon the integrated flux.
    bright_sources = int_flux > flux_thresh  # Boolean mask for bright sources

    # Filter sources
    ra = ra[bright_sources]
    dec = dec[bright_sources]
    peak_flux = peak_flux[bright_sources]
    peak_flux_err = peak_flux_err[bright_sources]  # Also filter flux errors
    int_flux = int_flux[bright_sources]  # Also filter flux errors
    int_flux_err = int_flux_err[bright_sources]  # Also filter flux errors
    local_rms = local_rms[bright_sources]  # Also filter flux errors

    # Convert to SkyCoord
    coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    # Convert input target position to SkyCoord
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    # Find sources within search radius
    sep = target_coord.separation(coords)
    matched = sep < search_radius

    # Extract matched sources
    ra = ra[matched]
    dec = dec[matched]
    peak_flux = peak_flux[matched]
    peak_flux_err = peak_flux_err[matched]
    int_flux = int_flux[matched]
    int_flux_err = int_flux_err[matched]
    local_rms = local_rms[matched]

    return ra, dec, peak_flux, peak_flux_err, int_flux, int_flux_err, local_rms




def find_sources_peak_int_beam(flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='gleam'):
    if name == 'gleam':
        with fits.open(gleam_catalog) as hdul:
            ra = hdul[1].data["RAJ2000"]
            dec = hdul[1].data["DEJ2000"]
            peak_flux = hdul[1].data["peak_flux_wide"]
            peak_flux_err = hdul[1].data["err_peak_flux_wide"]
            int_flux = hdul[1].data["int_flux_wide"]
            int_flux_err = hdul[1].data["err_int_flux_wide"]
            local_rms = hdul[1].data["local_rms_wide"]
            
            # # Extract beam sizes (update column names if needed)
            # beam_maj = hdul[1].data["beam_maj"]  # Major axis in arcsec
            # beam_min = hdul[1].data["beam_min"]  # Minor axis in arcsec

    elif name == 'tgss':  
        with fits.open(tgss_catalog) as hdul:
            ra = hdul[1].data["RA"]
            dec = hdul[1].data["DEC"]
            peak_flux = hdul[1].data["Peak_flux"]
            peak_flux_err = hdul[1].data["E_Peak_flux"]
            int_flux = hdul[1].data["Total_flux"]
            int_flux_err = hdul[1].data["E_Total_flux"]
            local_rms = hdul[1].data["RMS_noise"]
            
            # beam_maj = hdul[1].data["beam_maj"]
            # beam_min = hdul[1].data["beam_min"]

    else:  
        with fits.open(pybdsf_catalog) as hdul:
            ra = hdul[1].data["RA"]
            dec = hdul[1].data["DEC"]
            peak_flux = hdul[1].data["Peak_flux"]
            peak_flux_err = hdul[1].data["E_Peak_flux"]
            int_flux = hdul[1].data["Total_flux"]
            int_flux_err = hdul[1].data["E_Total_flux"]
            local_rms = hdul[1].data["Isl_rms"]
            
            source_maj = hdul[1].data["Maj"]*3600
            source_min = hdul[1].data["Min"]*3600
            beam_maj = hdul[1].data["Maj_img_plane"]*3600
            beam_min = hdul[1].data["Min_img_plane"]*3600
            
            

    # Apply flux threshold
    bright_sources = int_flux > flux_thresh

    # Filter sources
    ra = ra[bright_sources]
    dec = dec[bright_sources]
    peak_flux = peak_flux[bright_sources]
    peak_flux_err = peak_flux_err[bright_sources]
    int_flux = int_flux[bright_sources]
    int_flux_err = int_flux_err[bright_sources]
    local_rms = local_rms[bright_sources]
    beam_maj = beam_maj[bright_sources]
    beam_min = beam_min[bright_sources]

    # Convert to SkyCoord
    coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    # Find sources within search radius
    sep = target_coord.separation(coords)
    matched = sep < search_radius

    # Extract matched sources
    ra = ra[matched]
    dec = dec[matched]
    peak_flux = peak_flux[matched]
    peak_flux_err = peak_flux_err[matched]
    int_flux = int_flux[matched]
    int_flux_err = int_flux_err[matched]
    local_rms = local_rms[matched]
    beam_maj = beam_maj[matched]
    beam_min = beam_min[matched]

    return ra, dec, peak_flux, peak_flux_err, int_flux, int_flux_err, local_rms, beam_maj, beam_min
"""
def find_sources_peak(catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='gleam'):
    if name == 'gleam':
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAJ2000"]   # RA in degrees
            gleam_dec = hdul[1].data["DEJ2000"]  # DEC in degrees
            gleam_flux = hdul[1].data["peak_flux_wide"]  # Wideband peak flux density
            gleam_flux_int = hdul[1].data["int_flux_wide"]  # Wideband peak flux density
            gleam_flux_err = hdul[1].data["err_peak_flux_wide"]  # Flux error (check column name)
            gleam_local_rms = hdul[1].data["local_rms_wide"]  # Flux error (check column name)
    if name == 'gleamx':
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAdeg"]   # RA in degrees
            gleam_dec = hdul[1].data["DEdeg"]  # DEC in degrees
            gleam_flux = hdul[1].data["Fpwide"]  # Wideband peak flux density
            gleam_flux_int = hdul[1].data["Fintwide"]  # Wideband peak flux density
            gleam_flux_err = hdul[1].data["e_Fpwide"]  # Flux error (check column name)
            gleam_local_rms = hdul[1].data["local_rms_wide"]  # Flux error (check column name)
            

    elif name == 'tgss':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RA"]   # RA in degrees
            gleam_dec = hdul[1].data["DEC"]  # DEC in degrees
            gleam_flux = hdul[1].data["Peak_flux"]  # Total flux density
            gleam_flux_int = hdul[1].data["Total_flux"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Peak_flux"]  # Flux error (check column name)
            gleam_local_rms = hdul[1].data["RMS_noise"]  # Flux error (check column name)
            
    else:  # pybdsf
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RA"]   # RA in degrees
            gleam_dec = hdul[1].data["DEC"]  # DEC in degrees
            gleam_flux = hdul[1].data["Peak_flux"]  # Total flux density
            gleam_flux_int = hdul[1].data["Total_flux"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Peak_flux"]  # Flux error (check column name)
            gleam_local_rms = hdul[1].data["Isl_rms"]  # Flux error (check column name)

    # Apply flux threshold upon the integrated flux.
    bright_sources = gleam_flux_int > flux_thresh  # Boolean mask for bright sources

    # Filter sources
    gleam_ra = gleam_ra[bright_sources]
    gleam_dec = gleam_dec[bright_sources]
    gleam_flux = gleam_flux[bright_sources]
    gleam_flux_err = gleam_flux_err[bright_sources]  # Also filter flux errors
    gleam_local_rms = gleam_local_rms[bright_sources]  # Also filter flux errors

    # Convert to SkyCoord
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)

    # Convert input target position to SkyCoord
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    # Find sources within search radius
    sep = target_coord.separation(gleam_coords)
    matched = sep < search_radius

    # Extract matched sources
    gleam_ra = gleam_ra[matched]
    gleam_dec = gleam_dec[matched]
    gleam_flux = gleam_flux[matched]
    gleam_flux_err = gleam_flux_err[matched]
    gleam_local_rms = gleam_local_rms[matched]

    return gleam_ra, gleam_dec, gleam_flux, gleam_flux_err, gleam_local_rms
"""

def find_sources_spind(catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='gleam'):
    if name == 'gleam':
        with fits.open(catalog) as hdul:
            gleam_name = hdul[1].data["Name"]   # RA in degrees
            gleam_ra = hdul[1].data["RAJ2000"]   # RA in degrees
            gleam_dec = hdul[1].data["DEJ2000"]  # DEC in degrees
            gleam_flux = hdul[1].data["int_flux_wide"]  # Wideband peak flux density
            gleam_flux_err = hdul[1].data["err_int_flux_wide"]  # Flux error (check column name)
            gleam_flux_151 = hdul[1].data["int_flux_151"]  # Wideband peak flux density
            gleam_flux_151_err = hdul[1].data["err_int_flux_151"]  # Flux error (check column name)
            alpha = hdul[1].data["alpha"]  # Flux error (check column name)
            alpha_err = hdul[1].data["err_alpha"]  # Flux error (check column name)
            
            # print(gleam_ra, gleam_dec)
    elif name == 'gleamx1':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            hdr = hdul[1].header
            hdr['TNULL384'] = '---'
            hdr['TNULL385'] = '---'
            # hdul.flush()
            gleam_name = hdul[1].data["GLEAM-X"]
            gleam_ra = hdul[1].data["RAdeg"]   # RA in degrees
            gleam_dec = hdul[1].data["DEdeg"]  # DEC in degrees
            gleam_flux = hdul[1].data["Fintwide"]  # Total flux density
            gleam_flux_151 = hdul[1].data["Fint151"]  # Total flux density
            gleam_flux_err = hdul[1].data["e_Fintwide"]  # Flux error (check column name)
            gleam_flux_151_err = hdul[1].data["e_Fint151"]  # Flux error (check column name)
            alpha = hdul[1].data["alpha-SP"] 
            alpha_err = hdul[1].data["e_alpha-SP"]  # Flux error (check column name)

    elif name == 'tgss':  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RAdeg"]   # RA in degrees
            gleam_dec = hdul[1].data["DEdeg"]  # DEC in degrees
            gleam_flux = hdul[1].data["Stotal"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Stotal"]  # Flux error (check column name)

    elif name == 'nvss':
        with fits.open(catalog) as hdul:
            data = hdul[1].data
    
            # Combine RA and DEC parts into SkyCoord
            ra_str = [f"{ra_h}h{ra_m}m{ra_s}s" for ra_h, ra_m, ra_s in zip(data['RAh'], data['RAm'], data['RAs'])]
            sign = np.where(data['DE-'] == '-', -1, 1)  # Get the sign from 'DE-'
            dec_deg = sign * (data['DEd'] + data['DEm']/60 + data['DEs']/3600)
            
            coords = SkyCoord(ra=ra_str, dec=dec_deg*u.deg, frame='icrs')
    
            gleam_ra = coords.ra.deg
            gleam_dec = coords.dec.deg
    
            # Flux and uncertainty (in Jy)
            gleam_flux = data["S1.4"] / 1000.0  # Convert mJy to Jy
            gleam_flux_err = data["e_S1.4"] / 1000.0  # Convert mJy to Jy
            
            # print(gleam_ra, gleam_dec)

    else:  # TGSS or other catalogs
        with fits.open(catalog) as hdul:
            gleam_ra = hdul[1].data["RA"]   # RA in degrees
            gleam_dec = hdul[1].data["DEC"]  # DEC in degrees
            gleam_flux = hdul[1].data["Total_flux"]  # Total flux density
            gleam_flux_err = hdul[1].data["E_Total_flux"]  # Flux error (check column name)

    # Apply flux threshold
    bright_sources = gleam_flux > flux_thresh  # Boolean mask for bright sources

    # Filter sources
    gleam_name = gleam_name[bright_sources]
    gleam_ra = gleam_ra[bright_sources]
    gleam_dec = gleam_dec[bright_sources]
    gleam_flux = gleam_flux[bright_sources]
    gleam_flux_err = gleam_flux_err[bright_sources]  # Also filter flux errors
    gleam_flux_151 = gleam_flux_151[bright_sources]  # Also filter flux errors
    gleam_flux_151_err = gleam_flux_151_err[bright_sources]  # Also filter flux errors
    alpha = alpha[bright_sources]  # Also filter flux errors
    alpha_err = alpha_err[bright_sources]  # Also filter flux errors

    # Convert to SkyCoord
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)

    # Convert input target position to SkyCoord
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    # Find sources within search radius
    sep = target_coord.separation(gleam_coords)
    matched = sep < search_radius

    # Extract matched sources
    gleam_name = gleam_name[matched]
    gleam_ra = gleam_ra[matched]
    gleam_dec = gleam_dec[matched]
    gleam_flux = gleam_flux[matched]
    gleam_flux_err = gleam_flux_err[matched]
    gleam_flux_151 = gleam_flux_151[matched]
    gleam_flux_151_err = gleam_flux_151_err[matched]
    alpha = alpha[matched]
    alpha_err = alpha_err[matched]

    return gleam_ra, gleam_dec, gleam_flux, gleam_flux_err, gleam_flux_151, gleam_flux_151_err, alpha, alpha_err, gleam_name