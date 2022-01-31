### This script defines functions  from different existing image cutout
### services to retrieve images of WALLABY detections

import numpy
from astropy.table import Table
from astropy.io import fits
import requests
#import time
import pandas as pd
from io import StringIO
from tqdm.auto import tqdm

# PanSTARRS: This uses the convenience functions provided by PanSTARRS 
# here: https://ps1images.stsci.edu/ps1image.html
# PS pixel scale is 0.25 arcseconds/pixel

def ps1_getimages(ra,dec,size=240,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def ps1_geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = ps1_getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


# Recent updated version for bulk image download: 
# https://outerspace.stsci.edu/display/PANSTARRS/PS1+Image+Cutout+Service#PS1ImageCutoutService-BulkImageDownloadPythonScript
ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
def ps1_getimages_bulk(tra, tdec, size=240, filters="grizy", format="fits", imagetypes="stack"):
     
    """Query ps1filenames.py service for multiple positions to get a list of images
    This adds a url column to the table to retrieve the cutout.
     
    tra, tdec = list of positions in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    format = data format (options are "fits", "jpg", or "png")
    imagetypes = list of any of the acceptable image types.  Default is stack;
        other common choices include warp (single-epoch images), stack.wt (weight image),
        stack.mask, stack.exp (exposure time), stack.num (number of exposures),
        warp.wt, and warp.mask.  This parameter can be a list of strings or a
        comma-separated string.
 
    Returns an astropy table with the results
    """
     
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    # if imagetypes is a list, convert to a comma-separated string
    if not isinstance(imagetypes,str):
        imagetypes = ",".join(imagetypes)
    # put the positions in an in-memory file object
    cbuf = StringIO()
    cbuf.write('\n'.join(["{} {}".format(ra, dec) for (ra, dec) in zip(tra,tdec)]))
    cbuf.seek(0)
    # use requests.post to pass in positions as a file
    r = requests.post(ps1filename, data=dict(filters=filters, type=imagetypes),
        files=dict(file=cbuf))
    r.raise_for_status()
    tab = Table.read(r.text, format="ascii")
 
    urlbase = "{}?size={}&format={}".format(fitscut,size,format)
    tab["url"] = ["{}&ra={}&dec={}&red={}".format(urlbase,ra,dec,filename)
            for (filename,ra,dec) in zip(tab["filename"],tab["ra"],tab["dec"])]
    return tab


# SkyMapper data are retrieved using their SimpleImageAccess Service. A list
# of cutout urls can be extracted from the returned table. There seems to be
# a maximum on the size of the requested cutout. So, this could be an issue
# for large sources and will require some form of a mosaic
def skymapper_getcutouts(ra_list, dec_list, size):
    url_list = []
    for (ra,dec) in tqdm(zip(ra_list,dec_list), total=len(ra_list)):
        #print(f"Retrieving for RA = {ra}, Dec = {dec}")
        query_url = "https://api.skymapper.nci.org.au/public/siap/dr2/query?POS={},{}&SIZE={}&FORMAT=image/fits&BAND=u,v,g,r,i,z&INTERSECT=center&RESPONSEFORMAT=CSV".format(ra,dec,size)
        sm_table = pd.read_csv(query_url,delimiter=',',error_bad_lines = False)
        available_bands = list(sm_table['band'].unique())
        target_url = []
        for band in available_bands:
            sm_table_band = sm_table[sm_table['band']==band]
            sm_table_max_exptime = sm_table_band[(sm_table_band['exptime'] == sm_table_band['exptime'].max())]
            sm_table_max_size = sm_table_max_exptime[(sm_table_max_exptime['size'] == sm_table_max_exptime['size'].max())]
            target_url.append(sm_table_max_size['get_fits'].values[0])
        #print(f"Here are the urls retrieved so far = {target_url}")
        url_list.append(target_url)
    return url_list
# WISE, 2MASS, AND GALEX data are retrieved from the hips2fits service
# provided by CDS: https://alasky.u-strasbg.fr/hips-image-services/hips2fits
# These data come with a slightly different initial scaling. However, during
# testing the simple background subtract values from these data are virtually
# identical (<5% difference).
#
# WISE pixel scale is 1.37499998090796 arcseconds/pixel
# unWISE pixel scale is ~2.75 arcseconds/pixels
# 2MASS pixel scale is 1.0000000242 arcseconds/pixel
# GALEX pixel scale is 1.5 arcseconds/pixel

def unwise_cutouts(ra_list, dec_list, size):
    unwise_pix_scale = 2.75
    width = height = int(size*3600/unwise_pix_scale)
    fov = size/60.
    url_list = []
    for (ra,dec) in tqdm(zip(ra_list,dec_list), total=len(ra_list)):
        #print(f"Retrieving unWISE cutouts for RA = {ra}, Dec = {dec}")
        query_url = f'https://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips=CDS%2FP%2FunWISE%2Fcolor-W2-W1W2-W1&width={width}&height={height}&fov={size}&projection=TAN&coordsys=icrs&rotation_angle=0.0&ra={ra}&dec={dec}&format=fits'
        url_list.append(query_url)
    return url_list

def twomass_cutouts(ra_list, dec_list, size):
    twomass_pix_scale = 1.0000000242
    width = height = int(size*3600/twomass_pix_scale)
    fov = size/60.
    url_list = []
    for (ra,dec) in tqdm(zip(ra_list,dec_list), total=len(ra_list)):
        #print(f"Retrieving 2MASS cutouts for RA = {ra}, Dec = {dec}")
        query_url = f'https://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips=CDS%2FP%2F2MASS%2Fcolor&width={width}&height={height}&fov={size}&projection=TAN&coordsys=icrs&rotation_angle=0.0&ra={ra}&dec={dec}&format=fits'
        url_list.append(query_url)
    return url_list

def galex_cutouts(ra_list, dec_list, size):
    galex_pix_scale = 1.5
    width = height = int(size*3600/galex_pix_scale)
    fov = size/60.
    url_list = []
    for (ra,dec) in tqdm(zip(ra_list,dec_list), total=len(ra_list)):
        #print(f"Retrieving GALEX cutouts for RA = {ra}, Dec = {dec}")
        query_url = f'https://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips=CDS%2FGALEXGR6%2FAIS%2Fcolor&width={width}&height={height}&fov={size}&projection=TAN&coordsys=icrs&rotation_angle=0.0&ra={ra}&dec={dec}&format=fits'
        url_list.append(query_url)
    return url_list
