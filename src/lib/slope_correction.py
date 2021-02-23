# -*- coding: utf-8 -*-
"""
Radiometric slope correction of Sentinel-1 data on Google Earth Engine
by Andreas Vollrath
taken from https://github.com/ESA-PhiLab/radiometric-slope-correction
"""
import ee
import numpy as np

def slope_correction(image):
    '''This function applies the slope correction on a image of Sentinel-1 data
       
       :param collection: ee.Image of Sentinel-1
       :param elevation: name of ee.Image of DEM
       :param model: model to be applied (volume/surface)
       :param buffer: buffer in meters for layover/shadow amsk
        
        :returns: ee.Image
    '''
    elevation_str = 'USGS/SRTMGL1_003'
    model     = 'volume'
    buffer    = 0
    
    def _volumetric_model_SCF(theta_iRad, alpha_rRad):
        '''Code for calculation of volumetric model SCF
        
        :param theta_iRad: ee.Image of incidence angle in radians
        :param alpha_rRad: ee.Image of slope steepness in range
        
        :returns: ee.Image
        '''
        
        # create a 90 degree image in radians
        ninetyRad = ee.Image.constant(90).multiply(np.pi/180)
        
        # model
        nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan()
        denominator = (ninetyRad.subtract(theta_iRad)).tan()
        return nominator.divide(denominator) 
    
    
    def _surface_model_SCF(theta_iRad, alpha_rRad, alpha_azRad):
        '''Code for calculation of direct model SCF
        
        :param theta_iRad: ee.Image of incidence angle in radians
        :param alpha_rRad: ee.Image of slope steepness in range
        :param alpha_azRad: ee.Image of slope steepness in azimuth
        
        :returns: ee.Image
        '''
        
        # create a 90 degree image in radians
        ninetyRad = ee.Image.constant(90).multiply(np.pi/180)
        
        # model  
        nominator = (ninetyRad.subtract(theta_iRad)).cos()
        denominator = (alpha_azRad.cos()
          .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos()))

        return nominator.divide(denominator)


    def _erode(image, distance):
      '''Buffer function for raster

      :param image: ee.Image that shoudl be buffered
      :param distance: distance of buffer in meters
        
      :returns: ee.Image
      '''
      
      d = (image.Not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()))
    
      return image.updateMask(d.gt(distance))
    
    
    def _masking(alpha_rRad, theta_iRad, buffer):
        '''Masking of layover and shadow
        
        
        :param alpha_rRad: ee.Image of slope steepness in range
        :param theta_iRad: ee.Image of incidence angle in radians
        :param buffer: buffer in meters
        
        :returns: ee.Image
        '''
        # layover, where slope > radar viewing angle 
        layover = alpha_rRad.lt(theta_iRad).rename('layover')

        # shadow 
        ninetyRad = ee.Image.constant(90).multiply(np.pi/180)
        shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow')
        
        # add buffer to layover and shadow
        if buffer > 0:
            layover = _erode(layover, buffer)   
            shadow = _erode(shadow, buffer)  

        # combine layover and shadow
        no_data_mask = layover.And(shadow).rename('no_data_mask')
        
        return layover.addBands(shadow).addBands(no_data_mask)
                        
        
    def _correct(image):
        '''This function applies the slope correction and adds layover and shadow masks        
        '''        
        elevation = ee.Image(elevation_str)
        # get the image geometry and projection
        geom = image.geometry()
        proj = image.select(1).projection()
        
        # calculate the look direction
        heading = (ee.Terrain.aspect(image.select('angle'))
                                     .reduceRegion(ee.Reducer.mean(), geom, 1000)
                                     .get('aspect'))
                   

        # Sigma0 to Power of input image (not necessary, as we will use linear scaled images)
        #sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0))
        sigma0Pow = image.select(['VV', 'VH'])
        
        # the numbering follows the article chapters
        # 2.1.1 Radar geometry 
        theta_iRad = image.select('angle').multiply(np.pi/180)
        phi_iRad = ee.Image.constant(heading).multiply(np.pi/180)
        
        # 2.1.2 Terrain geometry
        alpha_sRad = ee.Terrain.slope(elevation).select('slope').multiply(np.pi/180).setDefaultProjection(proj).clip(geom)
        phi_sRad = ee.Terrain.aspect(elevation).select('aspect').multiply(np.pi/180).setDefaultProjection(proj).clip(geom)
        
        # we get the height, for export 
        height = elevation.setDefaultProjection(proj).clip(geom)
        
        # 2.1.3 Model geometry
        #reduce to 3 angle
        phi_rRad = phi_iRad.subtract(phi_sRad)

        # slope steepness in range (eq. 2)
        alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()

        # slope steepness in azimuth (eq 3)
        alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()

        # local incidence angle (eq. 4)
        theta_liaRad = (alpha_azRad.cos().multiply((theta_iRad.subtract(alpha_rRad)).cos())).acos()
        theta_liaDeg = theta_liaRad.multiply(180/np.pi)

        # 2.2 
        # Gamma_nought
        gamma0 = sigma0Pow.divide(theta_iRad.cos())
        #gamma0dB = ee.Image.constant(10).multiply(gamma0.log10()).select(['VV', 'VH'], ['VV_gamma0', 'VH_gamma0'])
        #ratio_gamma = (gamma0dB.select('VV_gamma0')
        #                .subtract(gamma0dB.select('VH_gamma0'))
        #                .rename('ratio_gamma0'))

        if model == 'volume':
            scf = _volumetric_model_SCF(theta_iRad, alpha_rRad)

        if model == 'surface':
            scf = _surface_model_SCF(theta_iRad, alpha_rRad, alpha_azRad)

        # apply model for Gamm0_f
        gamma0_flat = gamma0.divide(scf)
        #gamma0_flatDB = (ee.Image.constant(10)
        #                 .multiply(gamma0_flat.log10())
        #                 .select(['VV', 'VH'],['VV_gamma0flat', 'VH_gamma0flat'])
        #                )

        #masks = _masking(alpha_rRad, theta_iRad, buffer)

        # calculate the ratio for RGB vis
        #ratio_flat = (gamma0_flatDB.select('VV_gamma0flat')
        #                .subtract(gamma0_flatDB.select('VH_gamma0flat'))
        #                .rename('ratio_gamma0flat')
        #             )

        return (gamma0_flat
              .addBands(image.select('angle'))
              .copyProperties(image)
              .set('system:time_start', image.get('system:time_start'))
              )
    
    
    # run and return correction
    return _correct(image)