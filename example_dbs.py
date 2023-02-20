import sqlite3
import pandas as pd

conn1 = sqlite3.connect('ngc5044_example1.db')

c1 = conn1.cursor()

c1.execute('''CREATE TABLE IF NOT EXISTS wallcat (
		name TEXT,
		id INTEGER,
		x REAL,
		y REAL,
		z REAL,
		x_min REAL,
		x_max REAL,
		y_min REAL,
		y_max REAL,
		z_min REAL,
		z_max REAL,
		n_pix INTEGER,
		f_min REAL,
		f_max REAL,
		f_sum REAL,
		rel REAL,
		rms REAL,
		w20 REAL,
		w50 REAL,
		ell_maj REAL,
		ell_min REAL,
		ell_pa REAL,
		ell3s_maj REAL,
		ell3s_min REAL,
		ell3s_pa REAL,
		kin_pa REAL,
		err_x REAL,
		err_y REAL,
		err_z REAL,
		err_f_sum REAL,
		ra REAL,
		dec REAL,
		freq REAL,
		flag INTEGER,
		wm50 REAL,
		x_peak REAL,
		y_peak REAL,
		z_peak REAL,
		ra_peak REAL ,
		dec_peak REAL,
		freq_peak REAL)''')
		
conn1.commit()
#con2.close()

conn2 = sqlite3.connect('ngc5044_example2.db')

c2 = conn2.cursor()

c2.execute('''CREATE TABLE IF NOT EXISTS wallcat (
		name TEXT,
		id INTEGER,
		x REAL,
		y REAL,
		z REAL,
		x_min REAL,
		x_max REAL,
		y_min REAL,
		y_max REAL,
		z_min REAL,
		z_max REAL,
		n_pix INTEGER,
		f_min REAL,
		f_max REAL,
		f_sum REAL,
		rel REAL,
		rms REAL,
		w20 REAL,
		w50 REAL,
		ell_maj REAL,
		ell_min REAL,
		ell_pa REAL,
		ell3s_maj REAL,
		ell3s_min REAL,
		ell3s_pa REAL,
		kin_pa REAL,
		err_x REAL,
		err_y REAL,
		err_z REAL,
		err_f_sum REAL,
		ra REAL,
		dec REAL,
		freq REAL,
		flag INTEGER,
		wm50 REAL,
		x_peak REAL,
		y_peak REAL,
		z_peak REAL,
		ra_peak REAL ,
		dec_peak REAL,
		freq_peak REAL)''')
		
conn2.commit()
#conn2.close()


def insert_example(df, connection):
	df = df.astype({'name': 'string',
		'id': 'int',
		'x': 'float',
		'y': 'float',
		'z': 'float',
		'x_min': 'float',
		'x_max': 'float',
		'y_min': 'float',
		'y_max': 'float',
		'z_min': 'float',
		'z_max': 'float',
		'n_pix': 'int',
		'f_min': 'float',
		'f_max': 'float',
		'f_sum': 'float',
		'rel': 'float',
		'rms': 'float',
		'w20': 'float',
		'w50': 'float',
		'ell_maj': 'float',
		'ell_min': 'float',
		'ell_pa': 'float',
		'ell3s_maj': 'float',
		'ell3s_min': 'float',
		'ell3s_pa': 'float',
		'kin_pa': 'float',
		'err_x': 'float',
		'err_y': 'float',
		'err_z': 'float',
		'err_f_sum': 'float',
		'ra': 'float',
		'dec': 'float',
		'freq': 'float',
		'flag': 'int',
		'wm50': 'float',
		'x_peak': 'float',
		'y_peak': 'float',
		'z_peak': 'float',
		'ra_peak': 'float' ,
		'dec_peak': 'float',
		'freq_peak': 'float'})
	df.to_sql('wallcat',con = connection, if_exists='append',index=False)
	
wallaby = pd.read_csv('../ngc5044_dr1_catalog.csv')

wallaby_top10 = wallaby.head(10)

insert_example(wallaby,conn1)
insert_example(wallaby_top10,conn2)



