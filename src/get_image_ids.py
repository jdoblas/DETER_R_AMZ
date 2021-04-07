# -*- coding: utf-8 -*-
"""
This script will compute the available S1 images
between date1 and date2
@author: juand
"""
import ee
import pandas
from lib.sar_ee_utils import getS1dataFloat,getS1dataFloatByIngestionDate
from get_config import get_config


def get_image_ids(date1str,date2str,config):    
    options = config['image_selection']

    if options['selection_mode']=='ingestion_date': 
        print(f'Looking for new S1 images ingested between {date1str} and {date2str}')
    else:
        print(f'Looking for new S1 images adquired between {date1str} and {date2str}')
    
    AOI = ee.FeatureCollection(options['area_of_interest'])
    date1=ee.Date(date1str)
    date2=ee.Date(date2str).advance(1,'days')
    if options['selection_mode']=='ingestion_date':
        col=getS1dataFloatByIngestionDate(date1, date2).filterBounds(AOI).sort('system:time_start')    
    else:
        col=getS1dataFloat(date1, date2).filterBounds(AOI).sort('system:time_start')
    if options['include_S1B'] == "False":
         col=col.filterMetadata("platform_number", "equals", "A")
    # Avoids processing of images older than 2 months since date1
    col = col.filterDate(date1.advance(-2,'months'),date2)
    col_ids=pandas.DataFrame(col.aggregate_array('system:index').getInfo())
    if col_ids.size>0:
        print (f'{col_ids.size} images found:')
        print (col_ids.to_string(index=False))
    else:
        print ("No new images found")  

    return col_ids

if __name__ == '__main__':
    ee.Initialize()
    config = get_config('config.ini')
    get_image_ids('2021-04-01','2021-04-01',config)

