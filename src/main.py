# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 00:42:27 2020

@author: juand
"""
import datetime
import logging
import os
import sys
import ee
from get_config import get_config
from get_image_ids import get_image_ids
from get_polygons import get_polygons_from_asset
from get_warnings import get_raster_warnings
from lib.sar_ee_utils import execTask, ee_export_vector_silent

USAGE = f"Usage: python {sys.argv[0]} begindate enddate [<config file>]"


def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(USAGE)
    if len(args) < 2:
        raise SystemExit(USAGE)
    else:
        initial_data = args[0]
        end_data = args[1]
    if len(args) == 3:
        config_file = args[2]
    else:
        config_file = 'config.ini'
    config = get_config(config_file)

    # 0 - initialize EE and logging
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)
    main_start_time = datetime.datetime.now()
    main_start_time_f = main_start_time.strftime("%b %d %Y %H:%M:%S")
    print(f"Initiating run at {main_start_time_f}")

    try:
        ee.Initialize()
    except:
        ee.Authenticate()
        ee.Initialize()

    # 1 - get new images
    new_images_id_list = get_image_ids(initial_data, end_data, config)

    # 2 - Initiate main loop
    print("Starting detection")
    i = 1
    sar_tmp_mask = ee.FeatureCollection([])
    for image_id in new_images_id_list.values:
        image_id = image_id[0]
        print("Processing image " + str(i) + "/" + str(new_images_id_list.count()[0]) + ": " + image_id)
        img_start_time = datetime.datetime.now()

        # Detection
        raster_task = get_raster_warnings(image_id, sar_tmp_mask, config)
        asset = execTask(raster_task)

        # Vectorize detection raster asset
        print("Computing warning polygons on " + asset)
        vector_task, sar_tmp_mask = get_polygons_from_asset(asset, sar_tmp_mask, config)
        if vector_task:
            print("Exporting CR2 polygons as asset")
            execTask(vector_task)
        else:
            print("No CR2 polygons to export as asset")
        print("task finished")
        print("Img done")
        img_end_time = datetime.datetime.now()
        print("Elapsed time: ", img_end_time - img_start_time)
        print("-------------------------------------------")
        i = i + 1

    # 3 - Export results
    output_options = config['output']
    sar_tmp_mask_size = sar_tmp_mask.size().getInfo()
    if sar_tmp_mask_size > 0:
        print(f"Finished run. Found {sar_tmp_mask_size} CR2 polygons")
        print("Exporting results of the run as shapefile")
        ee_export_vector_silent(ee.FeatureCollection(sar_tmp_mask),
                                os.path.join(output_options['local_export_folder'],
                                             output_options[
                                                 'output_prefix'] + "_CR2_" + initial_data + "_" + end_data + ".shp"))
        # write trigger file
        with open(os.path.join(output_options['local_export_folder'], 'trigger.txt'), 'w') as fp:
            pass
        # Export to drive
        task = ee.batch.Export.table.toDrive(collection=sar_tmp_mask,
                                             description=output_options[
                                                             'output_prefix'] + "_CR2_" + initial_data + "_" + end_data,
                                             fileFormat="SHP", folder=output_options['gdrive_export_folder'])
        execTask(task)

    # 4 - Update SAR mask
    sar_tmp_mask_size_desm = sar_tmp_mask.filterMetadata('class', 'equals', 'DESMATAMENTO').size().getInfo()

    if sar_tmp_mask_size_desm > 0 and config['masks']['update_sar_mask'] == "True":
        sar_mask_asset = config['masks']['sar_mask_asset']
        print("Updating SAR mask")
        new_sar_mask = ee.Image(sar_mask_asset).blend(
            ee.Image(sar_tmp_mask
                     .filterMetadata('class', 'equals', 'DESMATAMENTO')
                     .map(lambda ft: ft.set('desm', 1))
                     .reduceToImage(['desm'], 'first'))) \
            .selfMask()
        update_task = ee.batch.Export.image.toAsset(image=new_sar_mask, description='run_update_sar_mask',
                                                    assetId=sar_mask_asset + '_tmp', scale=20,
                                                    maxPixels=10000000000000)
        execTask(update_task)
        try:
            ee.data.deleteAsset(sar_mask_asset + "_old")
        except:
            print("Warning: no old version of sar mask found")
        try:
            ee.data.renameAsset(sar_mask_asset, sar_mask_asset + "_old")
            ee.data.renameAsset(sar_mask_asset + "_tmp", sar_mask_asset)
        except:
            print("Warning: There was an error updating the CR2 mask. Please verify assets.")


    # 5 - End of run

    main_end_time = datetime.datetime.now()
    print("Normal end of job at ", main_end_time.strftime("%b %d %Y %H:%M:%S"))
    print("Elapsed time: ", main_end_time - main_start_time)


if __name__ == '__main__':
    main()
