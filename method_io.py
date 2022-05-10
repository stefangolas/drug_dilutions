import os
import sqlite3
import csv
from datetime import datetime as dt
from pace_util import method_local_dir
from pace_util import containing_dirname
from pace_util import ROBOID
from datetime import datetime
import time

from pyhamilton import Plate96
import logging

# This file is for all reading and writing all data related to this specific method that doesn't go to the robot directly.
# That's database init/io, manifest reading, things like that

def ensure_meas_table_exists(db_conn):
    '''
    Definitions of the fields in this table:
    lagoon_number - the number of the lagoon, uniquely identifying the experiment, zero-indexed
    filename - absolute path to the file in which this data is housed
    plate_id - ID field given when measurement was requested, should match ID in data file
    timestamp - time at which the measurement was taken
    well - the location in the plate reader plate where this sample was read, e.g. 'B2'
    measurement_delay_time - the time, in minutes, after the sample was pipetted that the
                            measurement was taken. For migration, we consider this to be 0
                            minutes in the absense of pipetting time values
    reading - the raw measured value from the plate reader
    data_type - 'lum' 'abs' or the spectra values for the fluorescence measurement
    '''
    c = db_conn.cursor()
    c.execute('''CREATE TABLE if not exists measurements
                (lagoon_number, filename, plate_id, timestamp, well, measurement_delay_time, reading, data_type, bacteria_id, barcode)''')
    db_conn.commit()

def db_add_plate_data(plate_data, data_type, plate, vessel_numbers, read_wells, controller_manifest, barcode):
    db_conn = sqlite3.connect(os.path.join(method_local_dir, containing_dirname + '.db'))
    ensure_meas_table_exists(db_conn)
    c = db_conn.cursor()
    timestamp = str(dt.now()) #NOTE: previously plate_data.header.time, though this was not very useful.
    for lagoon_number, read_well in zip(vessel_numbers, read_wells):
        lagoon_num_to_well = {x:x+8*(x//8) for x in range(48)}
        lagoon_id = plate.position_id(lagoon_num_to_well[lagoon_number])
        bacteria_id = controller_manifest[lagoon_id]['id']
        
        filename = plate_data.path
        plate_id = plate_data.header.plate_ids[0]

        well = plate.position_id(read_well)
        measurement_delay_time = 0.0
        reading = plate_data.value_at(*plate.well_coords(read_well))
        data = (lagoon_number, filename, plate_id, timestamp, well, measurement_delay_time, 
                 reading, data_type, bacteria_id, barcode)
        logging.info("attempting to c.execute")
        c.execute('INSERT INTO measurements VALUES (?,?,?,?,?,?,?,?,?,?)', data)
        logging.info("successfully c.execute")
    while True:
        logging.info("line 48 while true loop")
        try:
            logging.info("trying to add to database")
            db_conn.commit() # Unknown why error has been happening. Maybe Dropbox. Repeat until success.
            logging.info("successfully committed to database")
            break
        except (IOError, sqlite3.OperationalError):
            logging.info("except statement line 52")
            time.sleep(1)
    db_conn.close()
    
def db_add_whole_plate(plate_data, data_type): # Specifically when every lagoon well maps to every reader plate well 1:1
    db_add_plate_data(plate_data, data_type, Plate96(''), range(96), range(96))
    
manifest_filename = 'method_local_'+ROBOID+'\\controller_manifest.csv'
def read_manifest(cols_as_tuple=False):
    '''Reads in the current contents of a controller manifest; returns as dict'''
    controller_manifest = {}
    #output_file = 'method_local_'+ROBOID+'\\manifest_history\\controller_manifest_' + str(int(round(time.time()))) + '.csv'
    #timestamp = str(datetime.now())
    while True: # retry in case of simultaneous file access
        try:
            with open(manifest_filename, newline='') as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    if row:
                        if cols_as_tuple:
                            controller_manifest[row[0]] = tuple(row[1:])
                        else:
                            controller_manifest[row[0]]={}
                            controller_manifest[row[0]]['id'] = row[1]
                            controller_manifest[row[0]]['vol'] = row[2]
                            
                    #row.append(timestamp)
                    #write_obj.writerow(row)
            break
        except EnvironmentError:
            time.sleep(30)
            pass
    return controller_manifest