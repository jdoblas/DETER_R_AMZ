# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:38:21 2020

@author: juand
"""
import ee


def get_forest_mask(sar_tmp_mask, AOI, config):
    """
    """
    options = config['masks']
    static_mask_asset = options['static']
    updated_deforestation_mask = options['updated_deforestation']
    complementary_sar_mask = options['complementary_sar']
    ongoing_sar_mask = ee.Image(sar_tmp_mask \
                                .filterMetadata('class', 'equals', 'DESMATAMENTO')
                                .filterBounds(AOI)
                                .map(lambda ft: ft.set('desm', 1))
                                .reduceToImage(['desm'], 'first')).unmask(0, False)

    deforestationMask = ee.Image(static_mask_asset) \
        .blend(ee.Image(updated_deforestation_mask)) \
        .blend(ee.Image(complementary_sar_mask)) \
        .blend(ee.Image(ongoing_sar_mask)) \
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
