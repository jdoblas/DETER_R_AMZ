# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:41:32 2020

@author: juand
"""
import ee,logging,os
#import geemap
from lib.sar_ee_utils import ee_export_vector_silent,compute_pol_area

def get_polygons_from_asset(asset, sar_tmp_mask, config):
    
    options                = config  ['detection'] 
    output_options         = config  ['output']  
    general_scale          = float (options ['general_scale'])
    area_threshold         = float (options ['area_threshold'])
    confirmation_threshold = float (options ['confirmation_threshold'])
    intensity_threshold    = float (options ['intensity_threshold'])
    
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)
    
    img=ee.Image(asset)
    
    vector=img.reduceToVectors(reducer='mean',scale=general_scale,maxPixels=7919954990,eightConnected=True)\
        .map(compute_pol_area)\
        .filterMetadata('area_ha','greater_than',area_threshold)
    
    vector = vector.filterMetadata('intensity','not_greater_than',intensity_threshold).map(lambda ft: ft.set('class','DEGRADACAO')) \
        .merge(vector.filterMetadata('intensity','greater_than',intensity_threshold).map(lambda ft: ft.set('class','DEGRADACAO')))
    
    polygons_CR2=vector.filterMetadata('n_alerts','greater_than',confirmation_threshold)
    polygons_CR1=vector.filterMetadata('n_alerts','not_greater_than',confirmation_threshold)
    
    CR1_size=polygons_CR1.size().getInfo()    
    CR2_size=polygons_CR2.size().getInfo()
    
    print(f"Found {CR1_size} CR1 polygons and {CR2_size} CR2 polygons")
    
    # Ouput polygons
    asset_export_folder = output_options['asset_export_folder']
    output_prefix       = output_options['output_prefix']
    output_asset        = output_prefix + "_" + img_id + '_pol' 
    output_path         = asset_export_folder + '/' + output_asset
    
    if (CR2_size>0):
        try:
            ee.data.deleteAsset(output_path+"_CR2")
        except:
            pass         
        task = ee.batch.Export.table.toAsset(polygons_CR2, output_asset+"_CR2", output_path+"_CR2")
        
        sar_tmp_mask=sar_tmp_mask.merge(polygons_CR2)    
    else:
        task = None
    if (CR1_size>0) and output_options['export_polygons']=='True': 
        print ("Exporting CR1 warning polygons")
        ee_export_vector_silent(polygons_CR1, os.path.join(output_options['local_export_folder'],output_prefix+"_CR1.shp"))
    if (CR2_size>0) and output_options['export_polygons']=='True': 
        print ("Exporting CR2 warning polygons")   
        ee_export_vector_silent(polygons_CR2, os.path.join(output_options['local_export_folder'],output_prefix+"_CR2.shp"))
    return task, sar_tmp_mask


