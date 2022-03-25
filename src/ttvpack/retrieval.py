# ===WIP===
# Functions to search and retrieve light curves.

download_dir = "./.ttvpack/data/"

import lightkurve as lk
import numpy as np
import astropy.table
import os
import pandas as pd
import eleanor
import astropy
from astropy.coordinates import SkyCoord
from tess_stars2px import tess_stars2px_function_entry


__all__ = ["search_target", "get_target"]


def search_target(tic=None, epic=None, coords=None):
    """ Search MAST and TESSCut for all target data.
    
    This function uses lightkurve and MAST.
    
    Parameters
    ----------
    tic : str or int
        TIC number of the target.
    epic : str or int
        EPIC number of the target.
    coords : tuple or `astropy.coordinates.SkyCoord`
        SkyCoord object or tuple of decimal coordinates of the target.
        
    Returns
    ----------
    result : `lightkurve.SearchResult`
    """
    
    ### How much vetting of the inputs should i do??
    ### error if two or more different inputs are given?
    
    # Find the target coordinates from the given input
    if tic is not None:
        coords = tuple(eleanor.mast.coords_from_tic(tic)[0])
        
    elif coords is not None:
        if isinstance(coords, SkyCoord):
            coords = (coords.ra.deg, coords.dec.deg)
        tic = eleanor.mast.tic_from_coords(coords)[0]
        
    elif epic is None:
        print("Please enter a target identifier.")
        return # TODO: Error!
    

    # search for mast and tesscut products
    if (coords is not None) or (tic is not None):
        mast_search_result = lk.search_lightcurve(f"{coords[0]} {coords[1]}")
        tesscut_search_result = lk.search_tesscut("TIC"+str(tic))
    else:
        mast_search_result = lk.search_lightcurve(str(epic))
        tesscut_search_result = lk.search_tesscut(str(epic))
        
    # Return empty SearchResult if no data
    if ( len(mast_search_result) == 0 and len(tesscut_search_result) == 0 ):
        print("No data could be found for the specified target.")
        return lk.SearchResult(None)
    
    
    # Cast 'mission' column as string due to type mismatch
    mast_search_result.table.replace_column('mission', mast_search_result.table['mission'].astype(str))
    
    
    # Combine search results
    combined_search_table = astropy.table.vstack(
        [mast_search_result.table, tesscut_search_result.table], join_type='outer')
    
    result = lk.SearchResult(combined_search_table)
    
    return result
    
    

def get_target(tic=None, epic=None, coords=None, missions="all", exptime=None, by_number=None):
    """ Download all requested data for a given target.
    
    Parameters
    ----------
    tic : str or int
        TIC number of the target.
    epic : str or int
        EPIC number of the target.
    coords : tuple or `astropy.coordinates.SkyCoord`
        SkyCoord object or tuple of decimal coordinates of the target.
    missions : str, or list of str
        Missions to search for, specifying spacecraft and quarter/campaign/sector. e.g. ["K2 Campaign 05", "TESS Sector 44"]
    by_number? : int or list of int
        As an alternative specification, give the indices of the desired products in the SearchResult object.
    exptime : str
        TODO: figure out how to specify this if at all.
    Returns
    ----------
    lightcurves : lightkurve.LightCurveCollection object
    """
    
    ### Future functionality :: optionally, use everest to process 60 second k2 data
    ### Allow more specific chosing of products. 
    
    search_result = search_target(tic=tic, epic=epic, coords=coords)
    available_missions = set(search_result.mission)
    
    # Select one data product per campaign/sector, choose based on the author and then based on exposure time
    
    author_preference = ["SPOC", "TESScut", "EVEREST", "K2"] 
    product_list = []
    
    for mission in available_missions:
        mission_products = search_result[np.where(search_result.mission == mission)]

        # select based on author preferences
        for author in author_preference:
            if author in mission_products.author:
                mission_products = mission_products[np.where(mission_products.author == author)]
                continue
        # If none of the products have allowed authors
        if not any([(author in mission_products.author) for author in author_preference]):
            continue

        # select further based on exposure time
        min_exptime = np.argmin(mission_products.exptime)
        mission_products = mission_products[min_exptime]

        product_list.append(mission_products)

    
    # Create a SearchResult object with only the selected products
    product_table = astropy.table.vstack(
            [product.table for product in product_list], join_type='outer')
    download_products = lk.SearchResult(product_table)
    
    
    # display the list of selected products
    print("The following products were selected for download:")
    print(download_products.__repr__(), end='\n\n')
    
    
    # Separate the available light curves from the ffis
    eleanor_download_inds    = np.where(download_products.author == "TESScut")
    lightkurve_download_inds = np.where(download_products.author != "TESScut")

    eleanor_downloads    = download_products[eleanor_download_inds]
    lightkurve_downloads = download_products[lightkurve_download_inds]

    
    lightcurves = lk.LightCurveCollection([])
    
    # Download non-ffi data using lightkurve
    ### TODO: change download directories
    if len(lightkurve_downloads) > 0:
        lightkurve_lightcurves = lightkurve_downloads.download_all(download_dir=download_dir+"lightkurve_products/")
        
        for lc in lightkurve_lightcurves.data:
            lightcurves.append(lc)
        

    # Get ffi products using eleanor:
    if len(eleanor_downloads) > 0:

        sectors = [] # sectors of tess ffis to download
        for tess_mission in eleanor_downloads.mission:
            sector_str = tess_mission.removeprefix("TESS Sector ")
            sectors.append(int(sector_str))
            
        # Get the object coords    
        if epic:
            ra = search_result.ra[search_result.ra != 0.][0]  # get first non-zero element
            dec = search_result.dec[search_result.dec != 0.][0]
            coords = (ra, dec)
        elif tic:
            coords = tuple(eleanor.mast.coords_from_tic(tic)[0])
        else:
            if isinstance(coords, SkyCoord):    
                coords = (coords.ra.deg, coords.dec.deg)
            
        eleanor_lightcurves = get_tess_lightcurves(coords, sectors)
        
        for lc in eleanor_lightcurves.data:
            lightcurves.append(lc)
            
            
    # Resolve masked data sets:
    for i in range(len(lightcurves)):
        lightcurves[i] = lightcurves[i].remove_nans()
        
        for key in lightcurves[i].keys():
            if isinstance(lightcurves[i][key], astropy.utils.masked.Masked):
                lightcurves[i][key] = lightcurves[i][key].unmasked
        
    return lightcurves


def get_tess_lightcurves(coords, sectors, **kwargs):
    
    """ Obtain the lightcurves of a TESS target using FFIs.
    
    Parameters
    ----------
    coords : tuple of float
        Coordinates of a target.
    sectors : List of int
        List of sectors to get data from.
    **kwargs : dict, optional
        Extra keyword arguments passed on to eleanor.TargetData()
        
    Returns
    ----------
    lightcurves : lightkurve.LightCurveCollection object containing a light curve for each sector
    """
    
    # Update eleanor's maxsector variable if required
    if max(sectors) > eleanor.maxsector.maxsector:
        eleanor.Update(sector=max(sectors))
    
    # Check the requested sectors against the available sectors on tesscut
    tesspoint_result = tess_stars2px_function_entry("0000", coords[0], coords[1])
    tesspoint_sectors = tesspoint_result[3][tesspoint_result[3] <= eleanor.maxsector.maxsector]
    if -1 in tesspoint_sectors:
        tesspoint_sectors = []
    
    # Ensure the sectors are both avaliable and requested
    intersect = set(tesspoint_sectors).intersection(set(sectors))
    if intersect.issubset(sectors) and (intersect != set(sectors)):
        print(f"Tried to download TESS sectors {sectors}, but only sectors {tesspoint_sectors} are available according to tess-point.")
        print(f"Resolving to download sectors {list(intersect)}.")
        
    sectors = list(intersect)
    sectors.sort()
    
    
    lightcurves = lk.LightCurveCollection([])
    
    # make sure the coords can be matched to a tic id
    if len(eleanor.mast.crossmatch_by_position(coords, 0.01, 'Mast.Tic.Crossmatch')) > 0:
        eleanor_source = eleanor.multi_sectors(coords=coords, sectors=sectors, tc=True)

        for source in eleanor_source:
            data = eleanor.TargetData(source, **kwargs)
            lightcurves.append(data.to_lightkurve())

        ### save products with other data??
    else:
        print("Could not find a TESS target at the given coordinates. Skipping ffi data download.")
    
    return lightcurves


    
### Not used for anythin now ;(
def coords_from_epic(epic):
    """ Retrieve the coordinates of a target given an EPIC ID. 
    First searches the k2 campaign database, and if that fails queries Simbad.
    
    Parameters
    ----------
    epic : str or int
        EPIC number of the target.
    Returns
    ----------
    coords : tuple of float or `None`
        The ra and dec of the target object, or `None` if no coordinates found.
    """
    
    # Download k2 campaign database and search for target coordinates
    
    if os.path.isfile(download_dir+"GO_all_campaigns_to_date.csv"):
        k2_data = pd.read_csv(download_dir+"GO_all_campaigns_to_date.csv")
    else:
        k2_data_url = "https://keplergo.github.io/KeplerScienceWebsite/data/GO_all_campaigns_to_date.csv"
        k2_data = pd.read_csv(k2_data_url)
        k2_data.to_csv(path_or_buf=download_dir+"GO_all_campaigns_to_date.csv")
        
    star_data = k2_data.query(f'`EPIC ID` == {epic}')
    
    ra  = star_data['RA (J2000) [deg]'][star_data.index[0]]
    dec = star_data['Dec (J2000) [deg]'][star_data.index[0]]
    
    if pd.isna(ra):
        return None
    if isinstance(ra, str):
        if (not ra) or ra.isspace():
            return None    
    
    coords = (float(ra), float(dec))
    
    return coords



    
        