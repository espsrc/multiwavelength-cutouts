### This is the main script that will be run to retrieve cutouts and save
### them in the current directory

from numpy import timedelta64
#import requests
import retrieve_cutouts as ret_cut
#from astropy.io import fits
import pandas as pd
import time
from tqdm.auto import tqdm

t0 = time.time()

wallaby_cat = pd.read_csv('../subsample.csv')

###Below is some simple testing on a single position for the PanSTARRS data.
ra_test = wallaby_cat.ra.to_list()#[151.109056583226,151.414463079159]
dec_test = wallaby_cat.dec.to_list()#[-28.444156643408,-28.442512016763]
size_test = int(680.*4.)


#url_list = ret_cut.ps1_geturl(ra=ra_test, dec=dec_test, size=size_test,filters="grizy",format="fits")
#print(url_list)

#images_list = [fits.open(url) for url in url_list]
#for i, image in enumerate(images_list):
#        fits.writeto('test-{}.fits'.format(i),image[0].data, header=image[0].header, overwrite=True)


ps1_table = ret_cut.ps1_getimages_bulk(ra_test, dec_test,size_test,filters="grizy",imagetypes="stack")
#print(ps1_table['filename'])
#'''
for i, row in enumerate(tqdm(ps1_table, total=len(ps1_table))):
    ra = row['ra']
    dec = row['dec']
    projcell = row['projcell']
    subcell = row['subcell']
    filter = row['filter']

    fname = 'test-{}-.fits'.format(i)

    url = row['url']
    
    #print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
    #r = requests.get(url)
    #open(fname,"wb").write(r.content)
    #image = fits.open(url)
    #fits.writeto('test-{}-.fits'.format(i),image[0].data,header=image[0].header,overwrite=True)
print("{:.1f} s: retrieved FITS files for positions".format(time.time()-t0))
#print(ps1_table['url'])
#'''
print('Retrieving Skymapper Cutouts')
sm_urls = ret_cut.skymapper_getcutouts(ra_test,dec_test,10./60.)
#print(sm_urls)
print('Retrieving unWISE Cutouts')
unwise_urls = ret_cut.unwise_cutouts(ra_test, dec_test, 10./60.)
print('Retrieving 2MASS Cutouts')
twomass_urls = ret_cut.twomass_cutouts(ra_test, dec_test, 10./60.)
print('Retrieving GALEX Cutouts')
galex_urls = ret_cut.galex_cutouts(ra_test, dec_test, 10./60.)

#print(unwise_urls)
#print(galex_urls)
#print(twomass_urls)
#print("Completed Downloading in {} seconds".format(time.time() - t0))
#print(images_list)
#fits.writeto('test.fits',images[0].data)

