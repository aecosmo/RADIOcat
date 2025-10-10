from astropy.coordinates import SkyCoord
from astropy import units as u
from tqdm import tqdm
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, search_around_sky


def _read_and_filter_catalog(catalog_path, name, target_coord, search_radius):
    """
    Helper function to read, parse, and filter a single catalog.
    This consolidates the redundant code from the original function.
    """
    print(f"--- Reading catalog '{catalog_path}' (format: '{name}') ---")
    table = Table.read(catalog_path, format='fits')
    
    if name in ['chakraborty', 'arnab']:
        ra_col, dec_col = 'RAJ2000', 'DEJ2000'
    elif name == 'pybdsf' or name == 'lofar':
        ra_col, dec_col = 'RA', 'DEC'
    elif name in ['sirothia', 'garn', 'FIRST', 'NVSS']:
        sign = np.where(table['DE-'] == '-', -1, 1)
        final_ra = (table['RAh'].value * u.hour + 
                    table['RAm'].value * u.minute + 
                    table['RAs'].value * u.second)
        final_dec = (sign * table['DEd'].value * u.degree + 
                     table['DEm'].value * u.arcminute + 
                     table['DEs'].value * u.arcsecond)
        coords = SkyCoord(ra=final_ra, dec=final_dec)
        table['RA_deg'] = coords.ra.deg
        table['DEC_deg'] = coords.dec.deg
        ra_col, dec_col = 'RA_deg', 'DEC_deg'
    else:
        raise ValueError(f"Catalog format name '{name}' is not supported.")

    catalog_coords = SkyCoord(ra=table[ra_col], dec=table[dec_col], unit='deg')

    separations = target_coord.separation(catalog_coords)
    regional_mask = separations < search_radius
    
    filtered_table = table[regional_mask]
    print(f"Found {len(filtered_table)} sources within search radius.")
    
    return catalog_coords[regional_mask]


def sky_browse(ref_cat_path, src_cat_path, names, target_ra_hms, target_dec_dms, search_radius, match_radius):
    """
    Counts all sources in a source catalog that are close to each source in a reference catalog.

    This function first performs a cone search on both catalogs to select sources from a
    specific region of the sky. Then, for each source in the reference list, it counts
    how many sources from the source list are within the match_radius.
    
    This function is helpful to know how many counter sources are getting picked in a counter catalog
    for your catalog sources.

    Parameters:
    - ref_cat_path (str): Path to the reference catalog FITS file.
    - src_cat_path (str): Path to the source catalog FITS file.
    - names (list/tuple): A list of two strings with the format names of the 
                          reference and source catalogs (e.g., ['NVSS', 'pybdsf']).
    - target_ra_hms (str): Target RA for the initial cone search (e.g., '16h10m01s').
    - target_dec_dms (str): Target DEC for the initial cone search (e.g., '54d30m36s').
    - search_radius (Quantity): The radius of the sky region to analyze (e.g., 1 * u.deg).
    - match_radius (Quantity): The radius within which to count neighboring sources (e.g., 10 * u.arcsec).

    Returns:
    - astropy.table.Table: A summary table listing each reference source from the filtered region
                           and the number of neighbors it has in the source catalog.
    """

    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))

    ref_coords = _read_and_filter_catalog(ref_cat_path, names[0], target_coord, search_radius)
    src_coords = _read_and_filter_catalog(src_cat_path, names[1], target_coord, search_radius)

    if len(ref_coords) == 0 or len(src_coords) == 0:
        print("\n>>> One of the catalogs is empty after filtering. Cannot count neighbors.")
        return Table()

    print(f"\n--- Counting neighbors within {match_radius} ---")

    idx_ref, idx_src, sep2d, sep3d = search_around_sky(ref_coords, src_coords, match_radius)
    unique_ref_indices, neighbor_counts = np.unique(idx_ref, return_counts=True)
    
    print(f"Total pairings found: {len(idx_ref)}")
    print(f"Total unique reference sources with at least one neighbor: {len(unique_ref_indices)}\n")
    summary_table = Table()
    summary_table['Ref_Source_Index'] = unique_ref_indices
    summary_table['RA_ref'] = ref_coords[unique_ref_indices].ra.to_string(unit=u.hourangle, sep='hms', precision=2)
    summary_table['DEC_ref'] = ref_coords[unique_ref_indices].dec.to_string(unit=u.deg, sep='dms', precision=1)
    summary_table['Num_Neighbors'] = neighbor_counts

    summary_table.sort('Num_Neighbors', reverse=True)

    return summary_table


def cone_search(catalog, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name='chakraborty'):
    ''' Perform a cone search to find sources in src_catalog within search_radius of target coordinates. 
    Parameters: 
    - catalog: Reference catalog with RA and DEC columns.
    - flux_thresh (float): Minimum flux threshold for sources in src_catalog. 
    - target_ra_hms (str): Target RA in HMS format. 
    - target_dec_dms (str): Target DEC in DMS format. 
    - search_radius (Quantity): Search radius as an astropy Quantity (e.g., in arcseconds).     
                                It is your image radius.
    - name (str): Name of the reference source catalog for labeling. 
    Returns: 
    - matched_RA: List of matched RA values. 
    - matched_DEC: List of matched DEC values. 
    - matched_flux: List of matched flux values. 
    - matched_flux_err: List of matched flux errors. 
    '''
    
    print(f"--- Starting search in catalog '{catalog}' with name='{name}' ---")
    catalog_table = Table.read(catalog, format='fits')

    if name == 'chakraborty' or name == 'arnab':
        ra = catalog_table['RAJ2000']
        dec = catalog_table['DEJ2000']
        flux = catalog_table['TotalFlux'] * 1e-3
        flux_err = catalog_table['e_TotalFlux'] * 1e-3
    
    elif name == 'pybdsf':
        ra = catalog_table['RA']
        dec = catalog_table['DEC']
        flux = catalog_table['Total_flux']
        flux_err = catalog_table['E_Total_flux']
    
    elif name == 'lofar':
        ra = catalog_table['RA']
        dec = catalog_table['DEC']
        flux = catalog_table['Total_flux'] 
        flux_err = catalog_table['E_Total_flux'] 
   
    elif name == 'sirothia' or name == 'garn' or name == 'FIRST' or name == 'NVSS': # Grouping similar logic
        ref_RA_h = catalog_table['RAh']
        ref_RA_m = catalog_table['RAm']
        ref_RA_s = catalog_table['RAs']
        ref_DEC_d = catalog_table['DEd']
        ref_DEC_m = catalog_table['DEm']
        ref_DEC_s = catalog_table['DEs']
        sign = np.where(catalog_table['DE-'] == '-', -1, 1)
        
        signed_DEC_d_col = sign * ref_DEC_d

        final_ra = ref_RA_h.value * u.hour + ref_RA_m.value * u.minute + ref_RA_s.value * u.second
        final_dec = signed_DEC_d_col.value * u.degree + ref_DEC_m.value * u.arcminute + ref_DEC_s.value * u.arcsecond
        
        coords = SkyCoord(ra=final_ra, dec=final_dec)
        ra = coords.ra.deg
        dec = coords.dec.deg
        
        
        if name == 'sirothia':
            flux = catalog_table['Stot'] * 1e-3 
            flux_err = catalog_table['rms'] * 1e-6 # uJy to Jy
        elif name == 'garn':
            flux = catalog_table['Sint'] * 1e-3
            flux_err = catalog_table['e_Sint'] * 1e-3
        elif name == 'FIRST':
            flux = catalog_table['Fint'] * 1e-3
            flux_err = catalog_table['Rms'] * 1e-3
        elif name == 'NVSS':
            flux = catalog_table['S1.4'] * 1e-3
            flux_err = catalog_table['e_S1.4'] * 1e-3
    elif name == 'TGSS':
        ra = catalog_table['RAdeg']
        dec = catalog_table['DEdeg']
        flux = catalog_table['Stotal'] * 1e-3
        flux_err = catalog_table['e_Stotal'] * 1e-3
    else:
        raise ValueError(f"Catalog name '{name}' is not supported.")

    print(f"Total sources loaded: {len(flux)}")
    print(f"Flux threshold provided: {flux_thresh} Jy")
    if len(flux) > 0:
        print(f"Brightest source in catalog: {np.max(flux):.4f} Jy")
    
    bright_sources = flux > flux_thresh
    
    num_bright_sources = np.sum(bright_sources)
    print(f"Number of sources found above threshold: {num_bright_sources}")

    if num_bright_sources == 0:
        print("\n>>> No sources passed the flux filter. Check if your flux_thresh is too high.")
        return ra[bright_sources], dec[bright_sources], flux[bright_sources], flux_err[bright_sources]

    ra = ra[bright_sources]
    dec = dec[bright_sources]
    flux = flux[bright_sources]
    flux_err = flux_err[bright_sources]
    
    
    if name == 'pybdsf' or name == 'TGSS' or name == 'lofar':
        ref_coords = SkyCoord(ra=ra, dec=dec)
    else:
        ref_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    
    target_coord = SkyCoord(target_ra_hms, target_dec_dms, unit=(u.hourangle, u.deg))
    
    print(f"Target coordinate: {target_coord.to_string('hmsdms')}")
    print(f"Search radius: {search_radius}")

    sep = target_coord.separation(ref_coords)
    
    print(f"Closest source is at a separation of: {np.min(sep).to(u.arcsec):.2f}")
    
    matched = sep < search_radius
    num_matched_sources = np.sum(matched)
    print(f"Number of sources found within search radius: {num_matched_sources}")
    
    if num_matched_sources == 0:
        print("\n>>> No sources found within the search radius. Check if radius is too small or target coordinates are correct.")

    ra = ra[matched]
    dec = dec[matched]
    flux = flux[matched]
    flux_err = flux_err[matched]
    
    print(f"\n--- Search complete. Returning {len(ra)} sources. ---\n")
    return ra, dec, flux, flux_err


def matched_sources(ref_catfits, src_catfits, flux_thresh, target_ra_hms, target_dec_dms, search_radius, match_radius, names):
    '''
    Match sources from a source catalog to a reference catalog after performing a cone search on both.
    
    Parameters:
    - ref_catfits (str): Filename of the reference FITS catalog.
    - src_catfits (str): Filename of the source FITS catalog.
    - flux_thresh (float): Minimum flux threshold for sources.
    - target_ra_hms (str): Target RA in HMS format for the cone search.
    - target_dec_dms (str): Target DEC in DMS format for the cone search.
    - search_radius (Quantity): Search radius for the initial cone search.
    - match_radius (Quantity): Maximum separation to consider two sources a match.
    - names (list/tuple): A list or tuple of two strings with the names of the 
                          reference and source catalogs (e.g., ['NVSS', 'pybdsf']).
                          
    Returns:
    - Tuple of matched source properties from both catalogs.
    '''
   
    ref_RA, ref_DEC, ref_flux, ref_flux_err = cone_search(
        ref_catfits, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name=names[0]
    )
    src_RA, src_DEC, src_flux, src_flux_err = cone_search(
        src_catfits, flux_thresh, target_ra_hms, target_dec_dms, search_radius, name=names[1]
    )

    if len(src_RA) == 0 or len(ref_RA) == 0:
        print(">>> One of the catalogs returned no sources after filtering. Cannot perform match.")
        return (np.array([]),) * 8
    
    src_coords = SkyCoord(ra=src_RA, dec=src_DEC, unit='deg')
    ref_coords = SkyCoord(ra=ref_RA, dec=ref_DEC, unit='deg')

    idx_ref, d2d, d3d = src_coords.match_to_catalog_sky(ref_coords)

    matched_mask = d2d < match_radius
    
    num_final_matches = np.sum(matched_mask)
    print(f"\n--- Cross-Matching Results ---")
    print(f"Found {num_final_matches} final matches within {match_radius}.")

    # Get the indices for the source catalog that have a valid match.
    matched_src_indices = np.where(matched_mask)[0]
    
    # Get the corresponding indices for the matched reference catalog sources.
    matched_ref_indices = idx_ref[matched_mask]
    # Get the separation distance for the matched pairs.
    matched_sep = d2d[matched_mask].to(u.arcsec)
    
    # Use the indices to retrieve all data for the matched sources
    ref_catRA = ref_RA[matched_ref_indices]
    ref_catDEC = ref_DEC[matched_ref_indices]
    ref_catflux = ref_flux[matched_ref_indices]
    ref_catflux_err = ref_flux_err[matched_ref_indices]

    src_catRA = src_RA[matched_src_indices]
    src_catDEC = src_DEC[matched_src_indices]
    src_catflux = src_flux[matched_src_indices]
    src_catflux_err = src_flux_err[matched_src_indices]
    
    return ref_catRA, ref_catDEC, ref_catflux, ref_catflux_err, src_catRA, src_catDEC, src_catflux, src_catflux_err

def compute_offsets(src_ra, src_dec, ref_ra, ref_dec):
    """
    Compute the positional offsets between a source catalogue and a reference catalogue.

    Parameters
    ----------
    src_ra : array-like
        Right Ascension of sources in the source catalogue (degrees).
    src_dec : array-like
        Declination of sources in the source catalogue (degrees).
    ref_ra : array-like
        Right Ascension of sources in the reference catalogue (degrees).
    ref_dec : array-like
        Declination of sources in the reference catalogue (degrees).

    Returns
    -------
    delta_ra_arcsec : numpy.ndarray
        Offsets in Right Ascension (RA) in arcseconds, accounting for cos(dec).
    delta_dec_arcsec : numpy.ndarray
        Offsets in Declination (DEC) in arcseconds.
    
    Notes
    -----
    This function uses Astropy's `SkyCoord.spherical_offsets_to` method to compute
    angular separations between corresponding coordinates in two catalogues.  
    Offsets are returned in arcseconds for easier interpretation.
    """
    src_coords = SkyCoord(ra=src_ra, dec=src_dec, unit='deg')
    ref_coords = SkyCoord(ra=ref_ra, dec=ref_dec, unit='deg')
    delta_ra_cosdec, delta_dec = ref_coords.spherical_offsets_to(src_coords)
    return delta_ra_cosdec.to(u.arcsec).value, delta_dec.to(u.arcsec).value

#usecase:
#delta_ra1, delta_dec1 = compute_offsets(ref_ra, ref_dec, src_ra, src_dec)

def classify_srcs(catalog, threshold_sigma=3, output_filename=None):
    """
    Classifies sources as point-like or extended and generates a diagnostic plot.

    This function uses the ratio of integrated to peak flux (S_int / S_peak)
    to classify sources. A source is considered extended if the natural log of
    this ratio is significantly greater than zero, based on the propagated errors.

    Parameters:
    ----------
    catalog : dict or pandas.DataFrame
        The input source catalog. Must contain the keys: 'Peak_flux', 'Total_flux',
        'E_Total_flux', 'E_Peak_flux', and 'Isl_rms'.
    threshold_sigma : float, optional
        The number of sigma for the classification threshold. A source is
        extended if ln(S_int/S_peak) > threshold_sigma * propagated_error.
        Default is 3.
    output_filename : str or None, optional
        If a string is provided (e.g., 'source_plot.png'), the plot will be saved.
        If None, the plot will be displayed interactively. Default is None.

    Returns:
    -------
    S_int_by_S_peak : numpy.ndarray
        The calculated ratio of integrated to peak flux for all sources.
    S_peak_by_sigma_L : numpy.ndarray
        The calculated signal-to-noise ratio (peak flux / local RMS).
    mask_extended : numpy.ndarray (bool)
        A boolean mask that is True for extended sources and False for point-like sources.
    """

    S_peak = catalog['Peak_flux']
    S_int = catalog['Total_flux']
    sigma_S = catalog['E_Total_flux']
    sigma_S_peak = catalog['E_Peak_flux']
    sigma_L = catalog['Isl_rms']

    with np.errstate(divide='ignore', invalid='ignore'):
        S_int_by_S_peak = S_int / S_peak
        S_peak_by_sigma_L = S_peak / sigma_L
        propagated_error = np.sqrt((sigma_S / S_int)**2 + (sigma_S_peak / S_peak)**2)
        
        mask_extended = np.log(S_int_by_S_peak) > threshold_sigma * propagated_error
    
    mask_extended = np.nan_to_num(mask_extended, nan=False)
    plt.figure(figsize=(10, 7))

    plt.scatter(S_peak_by_sigma_L[mask_extended], S_int_by_S_peak[mask_extended],
                color='red', s=5, label=f'Extended Sources ({np.sum(mask_extended)})', alpha=0.7)

    plt.scatter(S_peak_by_sigma_L[~mask_extended], S_int_by_S_peak[~mask_extended],
                color='blue', s=5, label=f'Point Sources ({np.sum(~mask_extended)})', alpha=0.7)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$S_{\mathrm{peak}}/\sigma_L$ (Signal-to-Noise)', fontsize=14)
    plt.ylabel(r'$S_{\mathrm{int}}/S_{\mathrm{peak}}$', fontsize=14)
    
    # xticks = [5, 10, 20, 50, 100, 200, 500, 1000]
    # yticks = [0.5, 1, 5, 10, 20]
    # plt.xticks(xticks, [str(x) for x in xticks])
    # plt.yticks(yticks, [str(y) for y in yticks])
    
    plt.axhline(1, color='white', linestyle='--', linewidth=1)
    plt.grid(True, which='both', linestyle=':', linewidth=0.7, alpha=0.5)
    plt.title('Source Classification based on Flux Ratios', fontsize=16)
    plt.legend()
    plt.tight_layout()

    if output_filename:
        plt.savefig(output_filename, dpi=150)
        print(f"Plot saved to {output_filename}")
    else:
        plt.show()

    return S_int_by_S_peak, S_peak_by_sigma_L, mask_extended