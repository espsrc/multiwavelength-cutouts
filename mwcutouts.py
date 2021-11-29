### This is the main script that will be run to retrieve cutouts and save
### them in the current directory

from numpy import timedelta64
#import requests
import retrieve_cutouts as ret_cut
from astropy.io import fits
import time

t0 = time.time()



###Below is some simple testing on a single position for the PanSTARRS data.
ra_test = [112.39242095470075]
dec_test = [-30.002971422550402]
size_test = int(680.*4.)


#url_list = ret_cut.ps1_geturl(ra=ra_test, dec=dec_test, size=size_test,filters="grizy",format="fits")
#print(url_list)

#images_list = [fits.open(url) for url in url_list]
#for i, image in enumerate(images_list):
#        fits.writeto('test-{}.fits'.format(i),image[0].data, header=image[0].header, overwrite=True)


table = ret_cut.ps1_getimages_bulk(ra_test, dec_test,size_test,filters="grizy",imagetypes="stack")
print(table['filename'])
#'''
for i, row in enumerate(table):
    ra = row['ra']
    dec = row['dec']
    projcell = row['projcell']
    subcell = row['subcell']
    filter = row['filter']

    fname = 'test-{}-.fits'.format(i)

    url = row['url']
    print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
    #r = requests.get(url)
    #open(fname,"wb").write(r.content)
    image = fits.open(url)
    #fits.writeto('test-{}-.fits'.format(i),image[0].data,header=image[0].header,overwrite=True)
print("{:.1f} s: retrieved FITS files for positions".format(time.time()-t0))

#'''


print("Completed Downloading in {} seconds".format(time.time() - t0))
#print(images_list)
#fits.writeto('test.fits',images[0].data)