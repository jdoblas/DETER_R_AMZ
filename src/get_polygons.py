# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:41:32 2020

@author: juand
"""
import ee, logging, os
# import geemap
from lib.sar_ee_utils import ee_export_vector_silent, compute_pol_area


def get_polygons_from_asset(asset, detected_pols, config):
    options = config['detection']
    output_options = config['output']
    general_scale = float(options['general_scale'])
    area_threshold = float(options['area_threshold'])
    confirmation_threshold = float(options['confirmation_threshold'])
    intensity_threshold = float(options['intensity_threshold'])

    logging.getLogger('googleapicliet.discovery_cache').setLevel(logging.ERROR)

    img = ee.Image(asset)

    mask = img.select(0) \
                .addBands(ee.Image.pixelArea()) \
                .reduceConnectedComponents(ee.Reducer.sum()) \
                .lt(area_threshold*10000) \
                .unmask(0, False) \
                .Not()
    reducer = ee.Reducer.median().combine(ee.Reducer.mode(),"",True)
    vector = img.updateMask(mask)\
        .reduceToVectors(reducer=reducer, scale=general_scale, maxPixels=7919954990, eightConnected=True) \
        .map(compute_pol_area) \
        .filterMetadata('area_ha', 'greater_than', area_threshold)\
        .select(['area_ha', 'n_alerts_median', 'intensity_median', 'label', 'system:index', 'daydetec_mode'],
                ['area_ha', 'n_alerts', 'intensity', 'label', 'system:index', 'daydetec'])

    #print (vector.first().propertyNames().getInfo())

    vector = vector.filterMetadata('intensity', 'not_greater_than', intensity_threshold).map(
        lambda ft: ft.set('class', 'DEGRADATION')) \
        .merge(vector.filterMetadata('intensity', 'greater_than', intensity_threshold).map(
        lambda ft: ft.set('class', 'CLEAR_CUT')))

    detected_pols = detected_pols.merge(vector)

    polygons_CR2 = vector.filterMetadata('n_alerts', 'greater_than', confirmation_threshold)
    polygons_CR1 = vector.filterMetadata('n_alerts', 'not_greater_than', confirmation_threshold)

    CR1_size = polygons_CR1.size().getInfo()
    CR2_size = polygons_CR2.size().getInfo()

    print(f"Found {CR1_size} CR1 polygons and {CR2_size} CR2 polygons")

    # Output polygons
    asset_name = asset.split("/")[-1]

    if (CR1_size > 0) and output_options['export_img_polygons'] == 'True':
        print("Exporting CR1 warning polygons")
        ee_export_vector_silent(polygons_CR1,
                                os.path.join(output_options['local_export_folder'], asset_name + "_CR1.shp"))
    if (CR2_size > 0) and output_options['export_img_polygons'] == 'True':
        print("Exporting CR2 warning polygons")
        ee_export_vector_silent(polygons_CR2,
                                os.path.join(output_options['local_export_folder'], asset_name + "_CR2.shp"))
    return detected_pols


if __name__ == "__main__":
    import ee, configparser

    ee.Initialize()
    sar_tmp_mask_test = ee.FeatureCollection([])
    config_filename = 'config.ini'
    config_file = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'config' + os.sep + config_filename
    config = configparser.ConfigParser()
    config.read(config_file)
    test_asset = 'users/detersaree/DETER_R_2021_OUTPUT/DETER_R_AMZ_S1A_IW_GRDH_1SDV_20210413T094535_20210413T094604_037430_04697D_4087_raster'
    sar_tmp_mask_test = get_polygons_from_asset(test_asset, sar_tmp_mask_test, config)
