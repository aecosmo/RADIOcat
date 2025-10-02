import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u


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

def matched_sources(gleam_catalog, pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, match_radius, names):
    
    gleam_ra, gleam_dec, gleam_flux, gleam_flux_err = find_sources(gleam_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[0])
    pybdsf_ra, pybdsf_dec, pybdsf_flux, pybdsf_flux_err = find_sources(pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[1])
    
    
    # Convert to SkyCoord objects
    pybdsf_coords = SkyCoord(ra=pybdsf_ra * u.deg, dec=pybdsf_dec * u.deg)
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)
    
    # Find nearest neighbor in GLEAM catalog for each PyBDSF source
    idx_gleam, d2d, _ = pybdsf_coords.match_to_catalog_sky(gleam_coords)
    
    # Set matching radius (e.g., 45 arcseconds for GLEAM-X resolution)
    matched_mask = d2d < match_radius
    
    # Apply mask to get matched indices
    matched_pybdsf = np.where(matched_mask)[0]  # Indices in PyBDSF catalog
    matched_gleam = idx_gleam[matched_mask]      # Nearest neighbor indices in GLEAM catalog
    matched_sep = d2d[matched_mask].to(u.arcsec) # Separation in arcsec
    return pybdsf_ra[matched_pybdsf], pybdsf_dec[matched_pybdsf], gleam_ra[matched_gleam], gleam_dec[matched_gleam], pybdsf_flux[matched_pybdsf], \
    gleam_flux[matched_gleam],pybdsf_flux_err[matched_pybdsf], gleam_flux_err[matched_gleam]

def matched_sources_with_alpha(gleam_catalog, pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, match_radius, names):
    
    gleam_ra, gleam_dec, gleam_flux, gleam_flux_err, gleam_flux_151, gleam_flux_151_err, alpha, alpha_err, name = find_sources_spind(gleam_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[0])
    pybdsf_ra, pybdsf_dec, pybdsf_flux, pybdsf_flux_err = find_sources(pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[1])
    
    # argsort = np.argsort(-gleam_flux_151)
    # print(gleam_flux_151[argsort][:10])
    
    # Convert to SkyCoord objects
    pybdsf_coords = SkyCoord(ra=pybdsf_ra * u.deg, dec=pybdsf_dec * u.deg)
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)
    
    # Find nearest neighbor in GLEAM catalog for each PyBDSF source
    idx_gleam, d2d, _ = pybdsf_coords.match_to_catalog_sky(gleam_coords)
    
    # Set matching radius (e.g., 45 arcseconds for GLEAM-X resolution)
    matched_mask = d2d < match_radius*2
    
    # Apply mask to get matched indices
    matched_pybdsf = np.where(matched_mask)[0]  # Indices in PyBDSF catalog
    matched_gleam = idx_gleam[matched_mask]      # Nearest neighbor indices in GLEAM catalog
    matched_sep = d2d[matched_mask].to(u.arcsec) # Separation in arcsec
    return pybdsf_ra[matched_pybdsf], pybdsf_dec[matched_pybdsf], gleam_ra[matched_gleam], \
    gleam_dec[matched_gleam], pybdsf_flux[matched_pybdsf], gleam_flux[matched_gleam], \
    pybdsf_flux_err[matched_pybdsf], gleam_flux_err[matched_gleam], gleam_flux_151[matched_gleam], gleam_flux_151_err[matched_gleam],\
    alpha[matched_gleam], alpha_err[matched_gleam], name[matched_gleam]


def matched_sources_tgss(gleam_catalog, pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, match_radius, names):
    
    gleam_ra, gleam_dec, gleam_flux, gleam_flux_err = find_sources(gleam_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[0])
    pybdsf_ra, pybdsf_dec, pybdsf_flux, pybdsf_flux_err = find_sources(pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[1])
    
    
    # Convert to SkyCoord objects
    pybdsf_coords = SkyCoord(ra=pybdsf_ra * u.deg, dec=pybdsf_dec * u.deg)
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)
    
    # Find nearest neighbor in GLEAM catalog for each PyBDSF source
    idx_gleam, d2d, _ = pybdsf_coords.match_to_catalog_sky(gleam_coords)
    
    # Set matching radius (e.g., 45 arcseconds for GLEAM-X resolution)
    matched_mask = d2d < match_radius
    
    # Apply mask to get matched indices
    matched_pybdsf = np.where(matched_mask)[0]  # Indices in PyBDSF catalog
    matched_gleam = idx_gleam[matched_mask]      # Nearest neighbor indices in GLEAM catalog
    matched_sep = d2d[matched_mask].to(u.arcsec) # Separation in arcsec
    return pybdsf_ra[matched_pybdsf], pybdsf_dec[matched_pybdsf], gleam_ra[matched_gleam], gleam_dec[matched_gleam], pybdsf_flux[matched_pybdsf], \
    gleam_flux[matched_gleam], pybdsf_flux_err[matched_pybdsf], gleam_flux_err[matched_gleam]


def matched_sources_nvss(gleam_catalog, pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, match_radius, names):
    
    gleam_ra, gleam_dec, gleam_flux, gleam_flux_err = find_sources(gleam_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[0])
    pybdsf_ra, pybdsf_dec, pybdsf_flux, pybdsf_flux_err = find_sources(pybdsf_catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, names[1])
    
    # argsort = np.argsort(-gleam_flux_151)
    # print(gleam_flux_151[argsort][:10])
    
    # Convert to SkyCoord objects
    pybdsf_coords = SkyCoord(ra=pybdsf_ra * u.deg, dec=pybdsf_dec * u.deg)
    gleam_coords = SkyCoord(ra=gleam_ra * u.deg, dec=gleam_dec * u.deg)
    
    # Find nearest neighbor in GLEAM catalog for each PyBDSF source
    idx_gleam, d2d, _ = pybdsf_coords.match_to_catalog_sky(gleam_coords)
    
    # Set matching radius (e.g., 45 arcseconds for GLEAM-X resolution)
    matched_mask = d2d < match_radius*2
    
    # Apply mask to get matched indices
    matched_pybdsf = np.where(matched_mask)[0]  # Indices in PyBDSF catalog
    matched_gleam = idx_gleam[matched_mask]      # Nearest neighbor indices in GLEAM catalog
    matched_sep = d2d[matched_mask].to(u.arcsec) # Separation in arcsec
    return pybdsf_ra[matched_pybdsf], pybdsf_dec[matched_pybdsf], gleam_ra[matched_gleam], \
    gleam_dec[matched_gleam], pybdsf_flux[matched_pybdsf], gleam_flux[matched_gleam], \
    pybdsf_flux_err[matched_pybdsf], gleam_flux_err[matched_gleam]


