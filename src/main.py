# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 00:42:27 2020

@author: juand
"""
import ee
#import pandas
import os,sys
import configparser
import logging
import datetime,time
from get_image_ids import get_image_ids
from get_warnings import get_raster_warnings
from get_polygons import get_polygons_from_asset
from lib.sar_ee_utils import execTask #,ee_export_vector_silent

USAGE = f"Usage: python {sys.argv[0]} begindate enddate [<config file>]"

def read_config(args):
    if len(args)==3:
        config_filename=args[2]
    else:
        config_file_name='config.ini'
    config_file= os.path.dirname(os.path.abspath(__file__)) + os.sep + 'config' + os.sep + config_filename
    # Check config file
    if not os.path.exists(config_file):
        raise FileNotFoundError('File not found: ' + config_file + \
                            '. Please configure it first by copying/editing ' +  'config.ini.dist')
        exit(1)
    # Load config
    config = configparser.ConfigParser()
    config.read(config_file)
    return config


def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(USAGE)
    if len(args) < 2:
        raise SystemExit(USAGE)
    else:
        initial_data = args[0]
        end_data = args[1]

    config = read_config(args)
    
    # 0 - initialize EE and logging
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)
    main_start_time = datetime.datetime.now()
    main_start_time_f = main_start_time.strftime("%b %d %Y %H:%M:%S")
    print (f"Initiating run at {main_start_time_f}")

    try:
        ee.Initialize()
    except:
        ee.Authenticate()
        ee.Initialize()
        
    # 1 - get new images
    new_images_id_list = get_image_ids(initial_data, end_data, config)

    # 2 - Initiate main loop
    print ("Starting detection")
    i=1
    sar_tmp_mask=ee.FeatureCollection([])
    for image_id in new_images_id_list.values:
        image_id=image_id[0]
        print ("Processing image "+str(i)+"/"+str(new_images_id_list.count()[0])+": "+image_id)
        img_start_time=datetime.datetime.now()

        # Detection
        raster_task=get_raster_warnings(image_id, sar_tmp_mask, config)                                        
        asset = execTask(raster_task)

        # Vectorize detection raster asset
        print ("Computing warning polygons on " + asset)
        vector_task, sar_tmp_mask = get_polygons_from_asset(asset, sar_tmp_mask, config)
        if vector_task:
            print ("Exporting CR2 polygons as asset")
            execTask(vector_task)
        else:
            print ("No CR2 polygons to export as asset")
        print ("task finished")
        print ("Img done")
        img_end_time=datetime.datetime.now()
        print ("Elapsed time: ",img_end_time-img_start_time)
        print ("-------------------------------------------")
        i=i+1

    # 3 - Update SAR mask

    sar_tmp_mask_size=sar_tmp_mask.size().getInfo()
    if (sar_tmp_mask_size > 0):
        print (f"Finished run. Found {sar_tmp_mask_size} CR2 polygons")
        print ("Exporting results of the run as shapefile")
        ee_export_vector_silent(sar_tmp_mask, os.path.join(local_output_dir,output_prefix+"_CR2_"+initial_data+"_"+end_data+".shp"))      
        task=ee.batch.Export.table.toDrive(collection=sar_tmp_mask,description=output_prefix+"_CR2_"+initial_data+"_"+end_data,fileFormat="SHP",folder=gdrive_exports_path)
        execTask(task)    
        # UPDATE SAR_MASK
        if update_sar_mask:
            print ("Updating SAR mask")
            new_sar_mask=ee.Image(sar_mask_asset).blend(ee.Image(sar_tmp_mask.map(lambda ft:ft.set('desm',1)).reduceToImage(['desm'],'first'))).selfMask()
            update_task=ee.batch.Export.image.toAsset(image=new_sar_mask,description='run_update_sar_mask',assetId=sar_mask_asset+'_tmp',scale=20,maxPixels=10000000000000)
            execTask(update_task)
            try:
                ee.data.deleteAsset(sar_mask_asset+"_old")
            except:
                print("Warning: no old version of sar mask found")
            try:
                ee.data.renameAsset(sar_mask_asset,sar_mask_asset+"_old")
                ee.data.renameAsset(sar_mask_asset+"_tmp",sar_mask_asset)
            except:
                print ("Warning: There was an error updating the CR2 mask. Please verify assets.")
        # HERE FTP RESULTS of the run (will make it manual for now)
    else:
        print ("Not CR2 polygons on this run.")

    # 4 - End of run

    main_end_time=datetime.datetime.now()
    print ("Normal end of job at ",main_end_time.strftime("%b %d %Y %H:%M:%S"))
    print ("Elapsed time: ",main_end_time-main_start_time)
   
    config.detection=config['detection']
    config.masks=config['masks']
    config.output=config['output']
 

"""   
    
    
    
    
    
    
    
    with open(config_file) as f:
        params=yaml.load(f,Loader=yaml.FullLoader)
    detection_threshold=params.get('detection_threshold')
    area_threshold=params.get('area_threshold')
    confirmation_threshold=params.get('confirmation_threshold')
    AOI_shape_asset=params.get('area_of_interest_asset')
    adaptative_threshold=params.get('adaptative_threshold')
    output_assets_path=params.get('gee_assets_export_path')
    output_prefix=params.get('output_prefix')
    gdrive_exports_path=params.get('gdrive_exports_path')
    local_output_dir=params.get('local_export_folder')
    local_log_dir=params.get('local_log_folder')
    yearly_mask_asset=params.get('static_mask_asset')
    optical_mask_asset=params.get('ongoing_optical_mask_asset')
    sar_mask_asset=params.get('ongoing_sar_mask_asset')
    stabilize=params.get('stabilize')
    stabilization_type=params.get('stabilization_type')
    detection_mode=params.get('detection_mode')
    select_by_ingestion=params.get('select_images_by_ingestion_date')
    update_sar_mask=params.get('update_sar_mask')
    filter_type=params.get('filter_type')
    learning_period_years=params.get('learning_period_years')
    shrink=0
    export_footprints=1
    export_polygons=1
    ###### TESTING PARAMETERS (comment out for normal operation) ###########
    # detection_threshold=-1
    # area_threshold=.1
    # shrink=5000
    # export_footprints=0
    # export_polygons=0
    # stabilize=0    
    ######### INITIALIZATION #########
    #     # Initialize logging
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)
    main_start_time=datetime.datetime.now()
    main_start_time_f=main_start_time.strftime("%b %d %Y %H:%M:%S")
    print (f"Initiating run at {main_start_time_f}")
    try:
        ee.Initialize()
    except:
        ee.Authenticate()
        ee.Initialize()
    if not os.path.exists(local_log_dir):
        print ("Error: unable to find local log folder")
        sys.exit()
    if not os.path.exists(local_output_dir):
        print ("Error: unable to find local outputs folder")
        sys.exit()
    if yearly_mask_asset==None:
        print ("Error: you should specify at least a static non-forest mask GEE raster asset")
        sys.exit()
    if sar_mask_asset==None and update_sar_mask:
        print ("Error: you should specify at least a SAR deforestation mask GEE raster asset (can be an empty one)")
        sys.exit()
    if output_prefix==None:
        output_prefix="DETERSAR"
    if learning_period_years:
        learning_period_years=int(learning_period_years)
    if stabilization_type!="Harmonic": 
        stabilization_type="Spatial"
        if learning_period_years==False:
            learning_period_years=3
    else:
        if learning_period_years==False:
            learning_period_years=2
    if detection_mode:
        if (detection_mode.lower()=="mlc"):
            detection_mode="mlc"
            if learning_period_years==False:
                learning_period_years=2
        else:
            detection_mode="alt"
    else:
        detection_mode="alt"
    if (output_assets_path[-1]!="/"): output_assets_path=output_assets_path+"/"
    # If the output asset path don't exists, we create it
    try:
       ee.data.createAsset({'type':'Folder'},output_assets_path[0:-1])
    except:
        pass
    if filter_type==None:
        filter_type="QM5RL"
    ###### Find new images ######
    if select_by_ingestion==False: 
        print(f'Looking for new S1A images acquired between {initial_data} and {end_data}')
    else:
        print(f'Looking for new S1A images ingested between {initial_data} and {end_data}')
    new_images_id_list=get_new_images(initial_data,end_data,\
            AOI_shape_asset,os.path.join(local_log_dir,"new_images_"+initial_data+"_"+end_data+".txt"),include_S1B,select_by_ingestion)
    if new_images_id_list.size>0:
        print (f'{new_images_id_list.size} images found:')
        print (new_images_id_list.to_string(index=False))
    else:
        print ("No new images found")
    ###### MAIN LOOP ############
    print ("Starting detection")
    i=1
    sar_tmp_mask=ee.FeatureCollection([])
    for image_id in new_images_id_list.values:
        IMAGE_ID=image_id[0]
        print ("Processing image "+str(i)+"/"+str(new_images_id_list.count()[0])+": "+IMAGE_ID)
        img_start_time=datetime.datetime.now()
        assetId="RASTER_DETECTION_"+IMAGE_ID
        raster_task=get_raster_warnings(img_id=IMAGE_ID,\
                                        AOI_shape=AOI_shape_asset,\
                                        output_name=assetId,\
                                        output_asset_path=output_assets_path,\
                                        local_output_dir=local_log_dir,\
                                        detection_threshold=detection_threshold,\
                                        stabilize=int(stabilize),\
                                        stabilization_type=stabilization_type,\
                                        detection_mode=detection_mode,\
                                        shrink=shrink,\
                                        learningPeriodYears=learning_period_years,\
                                        yearly_mask=yearly_mask_asset,
                                        optical_mask=optical_mask_asset,
                                        sar_mask=sar_mask_asset,
                                        sar_tmp_mask=sar_tmp_mask,\
                                        export_footprints=export_footprints,\
                                        include_S1B=include_S1B,\
                                        adaptative_threshold=adaptative_threshold,\
                                        filter_type=filter_type,\
                                        output_prefix=output_prefix)
        execTask(raster_task)
        print ("Computing warning polygons")
        vector_task, sar_tmp_mask = get_polygons_from_asset(output_assets_path+assetId,\
                                        output_name=IMAGE_ID,\
                                        output_path=output_assets_path,\
                                        local_output_dir=local_log_dir,\
                                        area_threshold=area_threshold,\
                                        sar_tmp_mask=sar_tmp_mask,\
                                        export_polygons=export_polygons,\
                                        confirmation_threshold=confirmation_threshold,\
                                        output_prefix=output_prefix)
        if vector_task:
            print ("Exporting CR2 polygons as asset")
            execTask(vector_task)
        else:
            print ("No CR2 polygons to export as asset")
        print ("task finished")
        print ("Img done")
        img_end_time=datetime.datetime.now()
        print ("Elapsed time: ",img_end_time-img_start_time)
        print ("-------------------------------------------")
        i=i+1
    if (sar_tmp_mask.size().getInfo()>0):
        sar_tmp_mask_size=sar_tmp_mask.size().getInfo()
        print (f"Finished run. Found {sar_tmp_mask_size} CR2 polygons")
        print ("Exporting results of the run as shapefile")
        ee_export_vector_silent(sar_tmp_mask, os.path.join(local_output_dir,output_prefix+"_CR2_"+initial_data+"_"+end_data+".shp"))      
        task=ee.batch.Export.table.toDrive(collection=sar_tmp_mask,description=output_prefix+"_CR2_"+initial_data+"_"+end_data,fileFormat="SHP",folder=gdrive_exports_path)
        execTask(task)    
        # UPDATE SAR_MASK
        if update_sar_mask:
            print ("Updating SAR mask")
            new_sar_mask=ee.Image(sar_mask_asset).blend(ee.Image(sar_tmp_mask.map(lambda ft:ft.set('desm',1)).reduceToImage(['desm'],'first'))).selfMask()
            update_task=ee.batch.Export.image.toAsset(image=new_sar_mask,description='run_update_sar_mask',assetId=sar_mask_asset+'_tmp',scale=20,maxPixels=10000000000000)
            execTask(update_task)
            try:
                ee.data.deleteAsset(sar_mask_asset+"_old")
            except:
                print("Warning: no old version of sar mask found")
            try:
                ee.data.renameAsset(sar_mask_asset,sar_mask_asset+"_old")
                ee.data.renameAsset(sar_mask_asset+"_tmp",sar_mask_asset)
            except:
                print ("Warning: There was an error updating the CR2 mask. Please verify assets.")
        # HERE FTP RESULTS of the run (will make it manual for now)
    else:
        print ("Not CR2 polygons on this run.")
    # END OF RUN
    main_end_time=datetime.datetime.now()
    print ("Normal end of job at ",main_end_time.strftime("%b %d %Y %H:%M:%S"))
    print ("Elapsed time: ",main_end_time-main_start_time)
"""
if __name__ == '__main__':
    main()
