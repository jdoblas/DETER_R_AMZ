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
        initial_date = args[0]
        end_date =     args[1]
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
    new_images_id_list = get_image_ids(initial_date, end_date, config)

    # 2 - Initiate main loop
    print("Starting detection")
    i = 1
    detected_pols = ee.FeatureCollection([])
    for image_id in new_images_id_list.values:
        image_id = image_id[0]
        print("Processing image " + str(i) + "/" + str(new_images_id_list.count()[0]) + ": " + image_id)
        img_start_time = datetime.datetime.now()

        # Detection
        raster_task = get_raster_warnings(image_id, detected_pols, config)
        asset = execTask(raster_task)

        # Vectorize detection raster asset
        print("Computing warning polygons on " + asset)
        detected_pols = get_polygons_from_asset(asset, detected_pols, config)
        print("Img done")
        img_end_time = datetime.datetime.now()
        print("Elapsed time: ", img_end_time - img_start_time)
        print("-------------------------------------------")
        i = i + 1

    # 3 - Export results
    output_options = config['output']
    detected_pols_size = detected_pols.size().getInfo()

    if detected_pols_size > 0:

        confirmation_threshold = float(config['detection']['confirmation_threshold'])

        polygons_CR2 = detected_pols.filterMetadata('n_alerts', 'greater_than', confirmation_threshold)
        polygons_CR1 = detected_pols.filterMetadata('n_alerts', 'not_greater_than', confirmation_threshold)

        CR1_size = polygons_CR1.size().getInfo()
        CR2_size = polygons_CR2.size().getInfo()

        print(f"Finished run. Found {CR1_size} CR1 polygons and {CR2_size} CR2 polygons")
        print("Exporting results of the run as shapefile")
        ee_export_vector_silent(ee.FeatureCollection(polygons_CR2),
                                os.path.join(output_options['local_export_folder'],
                                             output_options[
                                                 'output_prefix'] + "_CR2_" + initial_date + "_" + end_date + ".shp"))
        ee_export_vector_silent(ee.FeatureCollection(polygons_CR1),
                                os.path.join(output_options['local_export_folder'],
                                             output_options[
                                                 'output_prefix'] + "_CR1_" + initial_date + "_" + end_date + ".shp"))
        # write trigger file
        with open(os.path.join(output_options['local_export_folder'], 'trigger.txt'), 'w') as fp:
            description = output_options['output_prefix'] + "_CR2_" + initial_date + "_" + end_date + ".shp"
            fp.write(description)
            fp.close()
        # Export to drive
        task = ee.batch.Export.table.toDrive(collection=polygons_CR2,
                                             description=output_options[
                                                             'output_prefix'] + "_CR2_" + initial_date + "_" + end_date,
                                             fileFormat="SHP", folder=output_options['gdrive_export_folder'])
        execTask(task)

        # 4 - Update SAR mask
        if CR2_size > 0 and config['masks']['update_sar_mask'] == "True":
            complementary_sar_col = config['masks']['complementary_sar_col']
            print("Updating SAR mask")
            polygons_CR2_CR_raster = ee.Image(polygons_CR2
                                              .filterMetadata('class', 'equals', 'CLEAR_CUT')
                                              .map(lambda ft: ft.set('desm', 1))
                                              .reduceToImage(['desm'], 'first'))\
                                            .set('system:time_start', ee.Date(initial_date).millis())\
                                            .set('system:time_end',   ee.Date(end_date).millis())
            try:
                ee.data.deleteAsset(complementary_sar_col + '/' + output_options['output_prefix'] + "_CR2_" + initial_date + "_" + end_date)
            except:
                pass
            update_task = ee.batch.Export.image.toAsset(image=polygons_CR2_CR_raster,
                                                        description='run_update_sar_mask',
                                                        assetId=complementary_sar_col + '/' + output_options['output_prefix'] + "_CR2_" + initial_date + "_" + end_date,
                                                        scale=int(config['detection']['general_scale']),
                                                        region=polygons_CR2.geometry().buffer(1000).bounds(),
                                                        maxPixels=10000000000000)
            execTask(update_task)

    # 5 - End of run

    main_end_time = datetime.datetime.now()
    print("Normal end of job at ", main_end_time.strftime("%b %d %Y %H:%M:%S"))
    print("Elapsed time: ", main_end_time - main_start_time)


if __name__ == '__main__':
    main()
