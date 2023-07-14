### This is the main script that will be run to retrieve cutouts and save
### them in the current directory.
###
### All files will need to be in the same directory.

#Python Packages
import sqlite3
import pandas as pd
import sys,logging
import argparse

#Wallaby Packages
import retrieve_cutouts as ret_cut
import db_functions as db_func
import catquery_functions as cqf

# CLI params and args parser
CLI=argparse.ArgumentParser()

# Input file
CLI.add_argument(
  "--input",  
  nargs=1, 
  type=str,
  default="./ngc5044_example2.db",  
)

# Catalogs
CLI.add_argument(
  "--catalogs",
  nargs="*",
  type=str, 
  default=["skymapper","unwise","twomass","galex","ls"],
)

args = CLI.parse_args()
inputfile = args.input[0]
catalogs = args.catalogs


logging.basicConfig(level=logging.DEBUG)

# UNCOMMENT ON FIRST RUN OF mwcutouts.py
db_func.create_wallaby_mw_tables()
logging.info("READING IN TABLES")

#To test this example for the first time, the db name should be "ngc5044_example2.db"
#and on subsequent runs use "ngc5044_example1.db"
ngc5044_example_conn = sqlite3.connect(inputfile)

wallaby_cat = pd.read_sql('SELECT * from wallcat LIMIT 20',ngc5044_example_conn)#pd.read_csv('../subsample.csv')
wallaby_conn = sqlite3.connect('wallaby_mw.db')
wallaby_prev = pd.read_sql('SELECT * from wallaby_mw_image_urls', wallaby_conn)

new_values_idx = db_func.compare_input_dbs(wallaby_prev,wallaby_cat)
logging.debug(new_values_idx)

if len(new_values_idx)>0:
	name_test = wallaby_cat.name.iloc[new_values_idx]
	ra_test = wallaby_cat.ra.iloc[new_values_idx]      
	dec_test = wallaby_cat.dec.iloc[new_values_idx]
	velo_test = cqf.freq_to_velo(wallaby_cat.freq.iloc[new_values_idx])
else:
	name_test = wallaby_cat.name
	ra_test = wallaby_cat.ra
	dec_test = wallaby_cat.dec
	velo_test = cqf.freq_to_velo(wallaby_cat.freq)
size_test = int(680.*4.)


logging.info("FETCHING IMAGE URLS")

ps1_table = ret_cut.ps1_getimages_bulk(ra_test, dec_test,size_test,filters="grizy",imagetypes="stack")
ps1_urls = ret_cut.ps1_merge_and_concat(ps1_table.to_pandas(),wallaby_cat)
logging.debug(ps1_urls)

if "skymapper" in catalogs:
	logging.info("Retrieving Skymapper Cutouts")
	sm_urls = ret_cut.skymapper_getcutouts(name_test, ra_test,dec_test,10.)
	logging.debug(sm_urls)

if "unwise" in catalogs: 
	logging.info('Retrieving unWISE Cutouts')
	unwise_urls = ret_cut.unwise_cutouts(name_test, ra_test, dec_test, 10.)
	logging.debug(unwise_urls)

if "twomass" in catalogs:
	logging.info('Retrieving 2MASS Cutouts')
	twomass_urls = ret_cut.twomass_cutouts(name_test, ra_test, dec_test, 10.)
	logging.debug(twomass_urls)

if "galex" in catalogs:
	logging.info('Retrieving GALEX Cutouts')
	galex_urls = ret_cut.galex_cutouts(name_test, ra_test, dec_test, 10.)
	logging.debug(galex_urls)

if "ls" in catalogs:
	logging.info('Retrieving LS-DR10 Cutouts')
	ls_urls = ret_cut.ls_cutouts(name_test, ra_test, dec_test, 10.)
	logging.debug(ls_urls)

emptydf = pd.DataFrame({'wallaby_id': [], 'url':[]})

merged_multi_df = ret_cut.merge_cutout_df(ps1_urls, 
					  sm_urls      if "skymapper" in catalogs else emptydf, 
					  unwise_urls  if "unwise" in catalogs else emptydf, 
					  twomass_urls if "twomass" in catalogs else emptydf,
					  galex_urls   if "galex" in catalogs else emptydf,
					  ls_urls      if "ls" in catalogs else emptydf, )
logging.debug(merged_multi_df)

exit()
logging.info('PERFORMING CROSS-MATCHING')

ned_df_list = []
sdss_df_list = []
for (name, ra, dec, velo) in zip(name_test, ra_test, dec_test, velo_test):
	ned_df_list.append(cqf.nedqueryandcheck_df(name, ra, dec, velo))
	sdss_df_list.append(cqf.sdssqueryandcheck(name, ra, dec, velo))
ned_df_total = pd.concat(ned_df_list,ignore_index=True)
sdss_df_total = pd.concat(sdss_df_list, ignore_index=True)
print(ned_df_total)
print(sdss_df_total)


sixd_cat = pd.read_csv('./6dFGSzDR3_cleaned_galaxiesonly.csv')
gswlc_cat = pd.read_csv('./GSWLC-X2_cleaned.csv')

sixd_df = cqf.cross_match_6df(sixd_cat, wallaby_cat)
gswlc_df = cqf.cross_match_gswlc(gswlc_cat,wallaby_cat)

logging.info('WRITING TO DATABASE TABLES')
db_func.insert_mw_image_urls(merged_multi_df.astype(str), wallaby_conn)
db_func.insert_ned(ned_df_total, wallaby_conn)
db_func.insert_sdss(sdss_df_total, wallaby_conn)
db_func.insert_6df(sixd_df, wallaby_conn)
db_func.insert_gswlc(gswlc_df,wallaby_conn)




