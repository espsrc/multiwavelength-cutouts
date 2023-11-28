import psycopg2

# Update the connection information for PostgreSQL
conn = psycopg2.connect(
    database="<db>",
    user="<user>",
    password="<password>",
    host="<host>",
    port="<port>"
)

c = conn.cursor()


def create_wallaby_mw_tables():
    with conn:
        try:
            c.execute("""SELECT * 
                        FROM wallaby_mw_image_urls 
                        """
                        )
        except:
            print("CREATING TABLES")
            c.execute("""CREATE TABLE wallaby_ned (
    wallaby_id VARCHAR,
    delta_velocity_flag INTEGER,
    delta_velocity REAL,
    "Number" INTEGER,
    "Object_Name" VARCHAR,
    "RA" REAL,
    "DEC" REAL,
    "Type" VARCHAR,
    "Velocity" REAL,
    "Redshift" REAL,
    "Redshift_Flag" INTEGER,
    "Magnitude_and_Filter" VARCHAR,
    "Separation" REAL,
    "Refs" VARCHAR,
    "Notes" VARCHAR,
    "Photometry_Points" VARCHAR,
    "Positions" VARCHAR,
    "Redshift_Points" VARCHAR,
    "Diameter_Points" VARCHAR,
    "Associations" VARCHAR
);""")

                c.execute("""CREATE TABLE wallaby_sdss (
    wallaby_id VARCHAR,
    delta_velocity_flag INTEGER,
    delta_velocity REAL,
    objID INTEGER,
    ra REAL,
    dec REAL,
    objtype VARCHAR,
    modelMag_u REAL,
    modelMag_g REAL,
    modelMag_r REAL,
    modelMag_i REAL,
    modelMag_z REAL,
    modelMagErr_u REAL,
    modelMagErr_g REAL,
    modelMagErr_r REAL,
    modelMagErr_i REAL,
    modelMagErr_z REAL,
    z_best REAL,
    zErr REAL
);""")

            c.execute("""CREATE TABLE wallaby_6dF (
    wallaby_id VARCHAR,
    delta_velocity REAL,
    angular_sep REAL,
    id VARCHAR,
    ra_hr INTEGER,
    ra_min INTEGER,
    ra_sec REAL,
    dec_deg INTEGER,
    dec_min INTEGER,
    dec_sec REAL,
    num_6d_meas INTEGER,
    num_fin_meas INTEGER,
    b_J REAL,
    progid INTEGER,
    r_F REAL,
    sg_class INTEGER,
    sum_flags INTEGER,
    v_best INTEGER,
    v_best_err INTEGER,
    v_source INTEGER,
    z_qual INTEGER,
    gal_lat REAL,
    gal_lon REAL,
    A_V REAL,
    weight_fib INTEGER,
    targid INTEGER,
    temp_code VARCHAR,
    z_file VARCHAR,
    specid VARCHAR,
    "RA" VARCHAR,
    "Dec" VARCHAR,
    RA_degrees REAL,
    Dec_degrees REAL
);""")

            c.execute("""CREATE TABLE wallaby_gswlc (
    wallaby_id VARCHAR,
    delta_velocity REAL,
    angular_sep REAL,
    objid REAL,
    glxid REAL,
    plate INTEGER,
    mjd INTEGER,
    fiberid INTEGER,
    ra REAL,
    dec REAL,
    z REAL,
    chi2_sed REAL,
    log_stellar REAL,
    log_stellar_err REAL,
    log_sfr REAL,
    log_sfr_err REAL,
    ext_fuv REAL,
    ext_fuv_err REAL,
    ext_b REAL,
    ext_b_err REAL,
    ext_v REAL,
    ext_v_err REAL,
    f_sed INTEGER,
    uv_survey INTEGER,
    f_uv INTEGER,
    f_midir INTEGER,
    f_mgs INTEGER,
    velo REAL
);""")

            c.execute("""CREATE TABLE wallaby_mw_image_urls (
    wallaby_id VARCHAR,
    panstarrs_url TEXT,
    skymapper_url TEXT,
    unwise_url TEXT,
    twomass_url TEXT,
    galex_url TEXT,
    ls_url TEXT
);""")
        else:
            print("TABLES EXIST, READY FOR NEW ENTRIES")
def insert_mw_image_urls(merged_df,connection):
    merged_df.to_sql('wallaby_mw_image_urls', con=connection, if_exists = 'append',index=False)
    
def insert_ned(merged_df,connection):
    merged_df.rename(columns = {"Object Name": "Object_Name", "No.": "Number",\
                   "Redshift Flag": "Redshift_Flag", "Magnitude and Filter": "Magnitude_and_Filter",\
                   "Photometry Points": "Photometry_Points", "Redshift Points": "Redshift_Points",\
                   "Diameter Points": "Diameter_Points", "References":"Refs"}, inplace=True)
    merged_df = merged_df.astype({'wallaby_id': 'string',
                    'delta_velocity_flag': 'int',
                    'delta_velocity': 'float',
                    'Number': 'int',
                    'Object_Name': 'string',
                    'RA': 'float',
                    'DEC': 'float',
                    'Type': 'string',
                    'Velocity': 'float',
                    'Redshift': 'float',
                    'Redshift_Flag': 'string',
                    'Magnitude_and_Filter': 'string',
                    'Separation': 'float',
                    'Refs': 'string',
                    'Notes': 'string',
                    'Photometry_Points': 'string',
                    'Positions': 'string',
                    'Redshift_Points': 'string',
                    'Diameter_Points': 'string',
                    'Associations': 'string'})
    merged_df.to_sql('wallaby_ned', con=connection, if_exists = 'append',index=False)
    
def insert_sdss(merged_df,connection):
    merged_df = merged_df.astype({'wallaby_id': 'string',
                    'delta_velocity_flag': 'int',
                    'delta_velocity': 'float',
                    'objID': 'int',
                    'ra': 'float',
                    'dec': 'float',
                    'objtype': 'string',
                    'modelMag_u': 'float',
                    'modelMag_g': 'float',
                    'modelMag_r': 'float',
                    'modelMag_i': 'float',
                    'modelMag_z': 'float',
                    'modelMagErr_u': 'float',
                    'modelMagErr_g': 'float',
                    'modelMagErr_r': 'float',
                    'modelMagErr_i': 'float',
                    'modelMagErr_z': 'float',
                    'z_best': 'float',
                    'zErr': 'float'})

    merged_df.to_sql('wallaby_sdss', con=connection, if_exists = 'append',index=False)
    
    
def insert_6df(merged_df, connection):
    merged_df.rename(columns={'RA_deg':'RA_degrees','Dec_deg':'Dec_degrees'},inplace=True)
    merged_df = merged_df.astype({'wallaby_id': 'string',
                    'delta_velocity': 'float',
                    'angular_sep': 'float',
                    'id': 'string',
                    'ra_hr': 'int',
                    'ra_min':'int',
                    'ra_sec':'float',
                    'dec_deg':'int',
                    'dec_min':'int',
                    'dec_sec':'float',
                    'num_6d_meas':'int',
                    'num_fin_meas':'int',
                    'b_J':'float',
                    'progid':'int',
                    'r_F':'float',
                    'sg_class':'int',
                    'sum_flags':'int',
                    'v_best':'int',
                    'v_best_err':'int',
                    'v_source':'int',
                    'z_qual':'int',
                    'gal_lat':'float',
                    'gal_lon':'float',
                    'A_V':'float',
                    'weight_fib':'int',
                    'targid':'int',
                    'temp_code':'string',
                    'z_file':'string',
                    'specid':'string',
                    'RA':'string',
                    'Dec':'string',
                    'RA_degrees':'float',
                    'Dec_degrees':'float'})
    merged_df.to_sql('wallaby_6dF',con = connection, if_exists = 'append', index=False)
    
def insert_gswlc(merged_df, connection):
    merged_df = merged_df.astype({'wallaby_id': 'string',
                    'delta_velocity': 'float',
                    'angular_sep': 'float',
                    'objid': 'float',
                    'glxid': 'float',
                    'plate': 'int',
                    'mjd': 'int',
                    'fiberid': 'int',
                    'ra': 'float',
                    'dec': 'float',
                    'z': 'float',
                    'chi2_sed': 'float',
                    'log_stellar': 'float',
                    'log_stellar_err': 'float',
                    'log_sfr': 'float',
                    'log_sfr_err': 'float',
                    'ext_fuv': 'float',
                    'ext_fuv_err': 'float',
                    'ext_b': 'float',
                    'ext_b_err': 'float',
                    'ext_v': 'float',
                    'ext_v_err': 'float',
                    'f_sed': 'int',
                    'uv_survey': 'int',
                    'f_uv': 'int',
                    'f_midir': 'int',
                    'f_mgs': 'int',
                    'velo': 'float'})
    merged_df.to_sql('wallaby_6dF',con = connection, if_exists = 'append', index=False)
    
    
def compare_input_dbs(db_table_1, db_table_2):
    db1_names = db_table_1.wallaby_id.to_list()
    db2_names = db_table_2.name.to_list()
    
    no_match_idx = [idx for idx, name in enumerate(db2_names) if name not in db1_names]
    return no_match_idx
