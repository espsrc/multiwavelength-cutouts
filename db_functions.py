import sqlite3

conn = sqlite3.connect('wallaby_mw.db')

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
                        wallaby_id text,
                        delta_velocity_flag integer,
                        delta_velocity real,
                        Number integer,
                        Object_Name text,
                        RA real,
                        DEC real,
                        Type text,
                        Velocity real,
                        Redshift real,
                        Redshift_Flag integer,
                        Magnitude_and_Filter text,
                        Separation real,
                        Refs text,
                        Notes text,
                        Photometry_Points text,
                        Positions text,
                        Redshift_Points text,
                        Diameter_Points text,
                        Associations text
                        )""")

            c.execute("""CREATE TABLE wallaby_sdss (
                        wallaby_id text,
                        delta_velocity_flag integer,
                        delta_velocity real,
                        objID integer,
                        ra real,
                        dec real,
                        objtype text,
                        modelMag_u real,
                        modelMag_g real,
                        modelMag_r real,
                        modelMag_i real,
                        modelMag_z real,
                        modelMagErr_u real,
                        modelMagErr_g real,
                        modelMagErr_r real,
                        modelMagErr_i real,
                        modelMagErr_z real,
                        z_best real,
                        zErr real
                        )""")

            c.execute(""" CREATE TABLE wallaby_6dF (
                        wallaby_id text,
                        delta_velocity real,
                        angular_sep real,
                        id text,
                        ra_hr integer,
                        ra_min integer,
                        ra_sec real,
                        dec_deg integer,
                        dec_min integer,
                        dec_sec real,
                        num_6d_meas integer,
                        num_fin_meas integer,
                        b_J real,
                        progid integer,
                        r_F real,
                        sg_class integer,
                        sum_flags integer,
                        v_best integer,
                        v_best_err integer,
                        v_source integer,
                        z_qual integer,
                        gal_lat real,
                        gal_lon real,
                        A_V real, 
                        weight_fib integer,
                        targid integer,
                        temp_code text,
                        z_file text,
                        specid text,
                        RA text,
                        Dec text,
                        RA_degrees real, 
                        Dec_degrees real
                        )""")

            c.execute(""" CREATE TABLE wallaby_gswlc(
                        wallaby_id text,
                        delta_velocity real,
                        angular_sep real,
                        objid real,
                        glxid real,
                        plate integer,
                        mjd integer,
                        fiberid integer,
                        ra real,
                        dec real,
                        z real,
                        chi2_sed real,
                        log_stellar real,
                        log_stellar_err real,
                        log_sfr real,
                        log_sfr_err real,
                        ext_fuv real,
                        ext_fuv_err real,
                        ext_b real,
                        ext_b_err real,
                        ext_v real,
                        ext_v_err real,
                        f_sed integer,
                        uv_survey integer,
                        f_uv integer,
                        f_midir integer,
                        f_mgs integer,
                        velo real
                        )""")

            c.execute(""" CREATE TABLE wallaby_mw_image_urls(
                        wallaby_id text,
                        panstarrs_url text,
                        skymapper_url text,
                        unwise_url text,
                        twomass_url text,
                        galex_url text,
                        ls_url text
                        )""")
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
