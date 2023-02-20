### This is the main script that will be run to retrieve cutouts and save
### them in the current directory

#Python Packages
import numpy as np
from numpy import timedelta64
import sqlite3
import pandas as pd
import time
from tqdm.auto import tqdm

#Wallaby Packages
import retrieve_cutouts as ret_cut
import db_functions as db_func
import catquery_functions as cqf



db_func.create_wallaby_mw_tables()

print('*'*10+'\n'+'*'*10+"\nREADING IN TABLES\n"+'*'*10+'\n'+'*'*10)


#To test this example for the first time, the db name should be "ngc5044_example2.db"
#and on subsequent runs use "ngc5044_example1.db"
ngc5044_example_conn = sqlite3.connect('ngc5044_example2.db')


wallaby_cat = pd.read_sql('SELECT * from wallcat LIMIT 20',ngc5044_example_conn)#pd.read_csv('../subsample.csv')
#display(wallaby_cat.head())
wallaby_conn = sqlite3.connect('wallaby_mw.db')
wallaby_prev = pd.read_sql('SELECT * from wallaby_mw_image_urls', wallaby_conn)

new_values_idx = db_func.compare_input_dbs(wallaby_prev,wallaby_cat)
print(new_values_idx)
if len(new_values_idx)>0:
	name_test = wallaby_cat.name.iloc[new_values_idx]    #.head(5).str.replace(r"\'",'',regex=True).to_list()
	ra_test = wallaby_cat.ra.iloc[new_values_idx]        #.head(5).to_list()#[151.109056583226,151.414463079159]
	dec_test = wallaby_cat.dec.iloc[new_values_idx]      #.head(5).to_list()#[-28.444156643408,-28.442512016763]
	velo_test = cqf.freq_to_velo(wallaby_cat.freq.iloc[new_values_idx])       #.head(5)).to_list()
else:
	name_test = wallaby_cat.name    #.head(5).str.replace(r"\'",'',regex=True).to_list()
	ra_test = wallaby_cat.ra        #.head(5).to_list()#[151.109056583226,151.414463079159]
	dec_test = wallaby_cat.dec      #.head(5).to_list()#[-28.444156643408,-28.442512016763]
	velo_test = cqf.freq_to_velo(wallaby_cat.freq)       #.head(5)).to_list()
size_test = int(680.*4.)

sixd_cat = pd.read_csv('../6dFGSzDR3_cleaned_galaxiesonly.csv')
gswlc_cat = pd.read_csv('../GSWLC-X2_cleaned.csv')

print('*'*10+'\n'+'*'*10+"\nFETCHING IMAGE URLS\n"+'*'*10+'\n'+'*'*10)

ps1_table = ret_cut.ps1_getimages_bulk(ra_test, dec_test,size_test,filters="grizy",imagetypes="stack")


ps1_urls = ret_cut.ps1_merge_and_concat(ps1_table.to_pandas(),wallaby_cat)
print(ps1_urls)
print('Retrieving Skymapper Cutouts')
sm_urls = ret_cut.skymapper_getcutouts(name_test, ra_test,dec_test,10.)
print(sm_urls)

print('Retrieving unWISE Cutouts')
unwise_urls = ret_cut.unwise_cutouts(name_test, ra_test, dec_test, 10.)
print(unwise_urls)

print('Retrieving 2MASS Cutouts')
twomass_urls = ret_cut.twomass_cutouts(name_test, ra_test, dec_test, 10.)
print(twomass_urls)

print('Retrieving GALEX Cutouts')
galex_urls = ret_cut.galex_cutouts(name_test, ra_test, dec_test, 10.)
print(galex_urls)

print('Retrieving LS-DR10 Cutouts')
ls_urls = ret_cut.ls_cutouts(name_test, ra_test, dec_test, 10.)
print(ls_urls)


merged_multi_df = ret_cut.merge_cutout_df(ps1_urls, sm_urls, unwise_urls, twomass_urls, galex_urls, ls_urls)
print(merged_multi_df.columns)

print('*'*10+'\n'+'*'*10+"\nPERFORMING CROSS-MATCHING\n"+'*'*10+'\n'+'*'*10)
ned_df_list = []
sdss_df_list = []
for (name, ra, dec, velo) in zip(name_test, ra_test, dec_test, velo_test):
	ned_df_list.append(cqf.nedqueryandcheck_df(name, ra, dec, velo))
	sdss_df_list.append(cqf.sdssqueryandcheck(name, ra, dec, velo))
ned_df_total = pd.concat(ned_df_list,ignore_index=True)
sdss_df_total = pd.concat(sdss_df_list, ignore_index=True)
print(ned_df_total)
print(sdss_df_total)
#ned_df_total.to_csv('test.csv',index=False)

sixd_df = cqf.cross_match_6df(sixd_cat, wallaby_cat)

gswlc_df = cqf.cross_match_gswlc(gswlc_cat,wallaby_cat)


print('*'*10+'\n'+'*'*10+"\nWRITING TO DATABASE TABLES\n"+'*'*10+'\n'+'*'*10)
db_func.insert_mw_image_urls(merged_multi_df.astype(str), wallaby_conn)
db_func.insert_ned(ned_df_total, wallaby_conn)
db_func.insert_sdss(sdss_df_total, wallaby_conn)
db_func.insert_6df(sixd_df, wallaby_conn)
db_func.insert_gswlc(gswlc_df,wallaby_conn)




