# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 17:13:46 2020

@author: juand
"""
import ee
import os, datetime
import pandas
import math
import logging
from lib.detrend import harmonic_detrending
from lib.sar_ee_utils import toDB, getS1dataFloat, refinedLeeFilter, toGamma0natural, \
    makeMedianFilter, QueganYuFilter, ee_export_vector_silent, \
    computeLogisticThreshold
# from lib.slope_correction import slope_correction
from get_masks import get_forest_mask, get_deforestation_mask
from timeit import default_timer as timer

def get_raster_warnings(img_id, detected_pols, config):
    # 0 - Initialize
    options = config['detection']
    output_options = config['output']
    general_scale = float(options['general_scale'])
    global_AOI = ee.FeatureCollection(config['image_selection']['area_of_interest'])
    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)

    # 1 - Get image from EE
    full_img_id = 'COPERNICUS/S1_GRD_FLOAT/' + img_id
    img = ee.Image(full_img_id)
    AOI = img.clip(global_AOI.geometry()).geometry()

    # 2 - Exports clipped footprint (optional)
    if output_options['export_footprints'] == 'True':
        print("Exporting footprint")
        ee_export_vector_silent(ee.FeatureCollection(AOI),
                                os.path.join(output_options['local_export_folder'],
                                             output_options['output_prefix'] + "_IMG_" + img_id + ".kml"))

    # 3 - Gets S1 collection    

    # Define relevant dates
    date3 = img.date().advance(1, 'days')
    date2 = date3.advance(-1 * int(options['detection_months']), 'months')
    date1 = date2.advance(-1 * int(options['learning_period_years']), 'years')
    # Define S1 colection from EE
    colS1 = getS1dataFloat(date1, date3) \
        .filterBounds(AOI) \
        .sort('system:time_start') \
        .map(lambda img: img.updateMask(img.gt(.00002))) \
        .map(toGamma0natural) \
        .select(["VHg0", "VVg0", "LIA"])

    if config['image_selection']['allow_orbit_overlap'] == 'False':
        colS1 = colS1.filterMetadata("relativeOrbitNumber_start", "equals", img.get('relativeOrbitNumber_start'))
    if config['image_selection']['include_S1B'] == 'False':
        colS1 = colS1.filterMetadata("platform_number", "equals", "A")
    if options['harmonic_detrend'] == 'True':
            colS1 = harmonic_detrending(colS1, 'VHg0')

        # First status messages
    print("Detecting warnings")
    print("Learning period start: " + date1.format('YYYY-MM-dd').getInfo())
    print("Detection start: " + date2.format('YYYY-MM-dd').getInfo())
    print("Detection end: " + date3.format('YYYY-MM-dd').getInfo())
    print("# of collected images: " + str(colS1.size().getInfo()))

    # 4 - Filter collection 

    median5 = makeMedianFilter(5)
    colS1_f = QueganYuFilter(colS1.select(0), median5).map(refinedLeeFilter)

    # add date band and filters out extreme Local Incidence Angles to avoid misdetections

    colS1_f2 = colS1_f.combine(colS1.select("LIA")) \
        .map(lambda img: img.updateMask(img.select("LIA").gt(28)).updateMask(img.select("LIA").lt(50))).select(0)

    # 5 - Detection

    print("Initializing detection")

    # Define collections
    forestMask = get_forest_mask(detected_pols, AOI, config)

    learnCol = colS1.select(0).filterDate(date1, date2)
    detectionCol = colS1_f2.select(0).map(toDB).filterDate(date2, date3) \
        .map(lambda img: img.updateMask(forestMask).addBands(
         ee.Image(img.date().difference("2020-01-01", "days")).int16().rename("julian_day")))

    # Compute detection threshold image
    pureDeforestationMask = get_deforestation_mask(AOI, config)
    if options['threshold_mode'] == 'logistic':
        ALT_threshold = computeLogisticThreshold(pureDeforestationMask, general_scale, float(options['threshold_min']),
                                                 float(options['threshold_max']))
    else:
        ALT_threshold = float(options['threshold_max'])
    mean_dif_mean_p1 = ee.Image(1.31)
    sd_dif_mean_p1 = ee.Image(0.35)
    detection_threshold_db = toDB(learnCol.median()) \
        .subtract(mean_dif_mean_p1.add(ALT_threshold.multiply(sd_dif_mean_p1)))

    # Detect warning areas
    detectionCol_db_below_threshold = detectionCol.select(0).map(lambda img: img.subtract(detection_threshold_db))
    detectionCol_warning = detectionCol_db_below_threshold.map(lambda img: img.lt(0))
    warning_count = detectionCol_warning.sum()
    warning_mask = warning_count.gte(1)

    # Compute 1st warning day
    nimg = detectionCol_warning.toArray().arrayProject([0]).arrayArgmax().arrayFlatten([["n"]])
    dd = ee.Image(detectionCol.select(['julian_day']).toArray().arrayGet(nimg.addBands(0)))

    # Morphological postprocessing of the warning raster
    contract_pixels = int(config['post_processing']['contract_pixels'])
    opening_pixels = int(config['post_processing']['contract_pixels'])

    detectionCol_warning_postprocessed = warning_mask \
        .unmask(0, False) \
        .focal_min(contract_pixels, "square").focal_max(contract_pixels, "square") \
        .focal_max(opening_pixels, "circle").focal_min(opening_pixels, "circle") \
        .selfMask() \
        .rename("alert") \
        .addBands(warning_count.rename("n_alerts").updateMask(warning_mask)) \
        .addBands(dd.rename("daydetec").updateMask(warning_mask)) \
        .addBands(detectionCol_db_below_threshold.min().multiply(-10).rename("intensity").updateMask(warning_mask))

    # 6 - Export of the results

    asset_export_folder = output_options['asset_export_folder']
    output_prefix = output_options['output_prefix']
    output_asset = output_prefix + "_" + img_id + '_raster'
    output_path = asset_export_folder + '/' + output_asset
    try:
        ee.data.deleteAsset(output_path)
    except:
        pass

    local_task = ee.batch.Export.image.toAsset(detectionCol_warning_postprocessed, \
                                               scale=general_scale, \
                                               assetId=output_path, \
                                               description=output_asset, \
                                               region=AOI.bounds().getInfo()['coordinates'], maxPixels=1587244950)
    return local_task


if __name__ == "__main__":
    import ee, configparser
    from lib.sar_ee_utils import execTask

    ee.Initialize()
    sar_tmp_mask = ee.FeatureCollection([])
    config_filename = 'config.ini'
    config_file = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'config' + os.sep + config_filename
    config = configparser.ConfigParser()
    config.read(config_file)
    test_img = 'S1A_IW_GRDH_1SDV_20190103T092406_20190103T092431_025311_02CCE7_CE4E'
    task = get_raster_warnings(test_img, sar_tmp_mask, config)
    execTask(task)

