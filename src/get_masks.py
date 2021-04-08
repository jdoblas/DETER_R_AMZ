# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:38:21 2020

@author: juand
"""
import ee
from get_config import get_config

def get_forest_mask(detected_pols, AOI, config):
    """
    """
    options = config['masks']
    static_mask_asset = options['static']
    updated_deforestation_mask = options['updated_deforestation']
    complementary_sar_mask_col = options['complementary_sar_col']

    clear_cut_detected_pols = detected_pols.filterMetadata('intensity', 'greater_than', float(config['detection']['intensity_threshold'])) \
                                        .filterMetadata('n_alerts', 'greater_than', float(config['detection']['confirmation_threshold']))\
                                        .filterBounds(AOI)

    if clear_cut_detected_pols.size().getInfo() > 0:
        ongoing_sar_mask = ee.Image(clear_cut_detected_pols
                                    .map(lambda ft: ft.set('desm', 1))
                                    .reduceToImage(['desm'], 'first'))
    else:
        ongoing_sar_mask = ee.Image(0).selfMask()

    deforestationMask = ee.Image(static_mask_asset) \
        .blend(ee.Image(updated_deforestation_mask).toByte()) \
        .blend(ee.Image(ongoing_sar_mask)) \
        .blend(ee.ImageCollection(complementary_sar_mask_col)
               .map(lambda img: img.unmask(0, False))
               .reduce(ee.Reducer.anyNonZero())
               .selfMask()) \
        .unmask(0, False)
    forestMask = deforestationMask.Not().clip(AOI)

    return forestMask


def get_deforestation_mask(AOI, config):
    """
    """
    options = config['masks']
    historic_mask_asset = options['historic_deforestation']
    updated_deforestation_mask = options['updated_deforestation']

    deforestationMask = ee.Image(historic_mask_asset) \
        .blend(ee.Image(updated_deforestation_mask)) \
        .unmask(0, False) \
        .clip(AOI)
    return deforestationMask

if __name__ == '__main__':
    ee.Initialize()
    config = get_config('config.ini')
    featCol = ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[-51.74386812876962, -3.5601460002191527],
                 [-51.750734583847745, -3.5717963731411193],
                 [-51.73837496470712, -3.573166995537336],
                 [-51.74386812876962, -3.5601460002191527]]),
            {
              "system:index": "0", "intensity": 11, "n_alerts": 2
            })]);
    geom = ee.Geometry.Polygon(
        [[[-51.837789638785566, -3.530825394959526],
          [-51.837789638785566, -3.642872000212739],
          [-51.7351361353676, -3.642872000212739],
          [-51.7351361353676, -3.530825394959526]]], None, False);
    mask = get_forest_mask(featCol, geom, config)
    params = {'min':0, 'max':1, 'palette':['brown','green']}
    print (mask.getThumbURL(params))

