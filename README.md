# WALLABY Multiwavelength cutouts

*This repo is still in development.*

To test this on your local machine you'll need the following packages that do not come bundled with your python installation:
- numpy
- pandas
- astropy
- tqdm

Once these are set up, you should be able to replicate the database/table creation and updates with the following steps:

1. Create the "example" databases that will be used in our examples by running the "example_dbs.py" script. You will need to have the "ngc5044_dr1_catalog.csv" file in the same directory to do so. If you have the DB browser installed, you should be able to see these dbs and tables.

2. Once the example databases are created, the primary multiwavelength database and tables can be created and updated. When running "mwcutouts.py" for the first time, ensure that on line 28 the connected database name is "ngc5044_example.db" and that line 21 (with 'db_func.create_wallaby_mw_tables()') is uncommented. Now, run "mwcutouts.py". This will initialize the database and tables with entries, mimicking the first inputs that will be provided to this table. Again, you should be able to view the 5 different tables that are created within this database:
     - NED
     - SDSS
     - 6dF survey
     - GSWLC

3. With the db and tables initialized, we can run a mock "update" of the tables with new entries. First, comment out line 21 (with 'db_func.create_wallaby_mw_tables()')change the name of the database on line 28 to "ngc5044_example1.db". Now, you should set to run "mwcutouts.py" to perform this "update". You should be able to see that only 10 additional entries were added since the first 10 in the "update" are identical to the initialized values.
