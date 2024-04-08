#!/usr/bin/env python

"""
The initial version of this file is derived from the UMRider Package developed at NCMRWF by Arulalan and team.

Only some of the general use library functions of UMRide is available in this file. 

Seperate Package for UMRider is available at https://github.com/NCMRWF/UMRider

__author__ = 'arulalant'
__version__ = 'v3.1.0'
__long_name__ = 'NCUM Parallel Rider'

References:
1. Iris. v1.8.1 03-Jun-2015. Met Office. UK. https://github.com/SciTools/iris/archive/v1.8.1.tar.gz

2. Saji M. (2014), "Utility to convert UM fieldsfile output to NCEP GRIB1 format:
                    A User Guide", NMRF/TR/01/2014, April 2014, pp. 51, available at
                    http://www.ncmrwf.gov.in/umfld2grib.pdf

Disclaimers (if any!)
This is just test code as of now and is meant for a specific purpose only!

Copyright: ESSO-NCMRWF, MoES, 2015-2016, 2016-2017.

Author : Arulalan.T
previous Update : 11-Mar-2016
latest Update : 29-Aug-2016
"""



def getCubeAttr(tmpCube):
    """
    This module returns basic coordinate & attribute info about any Iris data cube.
    :param tmpCube: a temporary Iris cube containing a single geophysical field/parameter
    :return: stdNm: CF-compliant Standard name of the field/parameter
    :return: fcstTm: forecast time period for eg: 00, 06, 12 etc -- units as in hours
    :return: refTm: reference time -- units as date  in Gregorian
    :return: lat as scalar array (1D) units as degree (from 90S to 90N)
    :return: lon as scalar array (1D) units as degree (from 0E to 360E)
    Original by MNRS
    """
    stdNm = tmpCube.standard_name
    stdNm = stdNm if stdNm else tmpCube.long_name
    stash = str(tmpCube.attributes['STASH'])
    fcstTm = tmpCube.coord('forecast_period')
    refTm = tmpCube.coord('forecast_reference_time')
    lat = tmpCube.coord('latitude')
    lon = tmpCube.coord('longitude')

    return stdNm, stash, fcstTm, refTm, lat, lon

# create a class #2 to initiate mp daemon processes
class _NoDaemonProcess(mp.Process):
    # make 'daemon' attribute always return False
    # A class created by AAT
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)
# end of class #2

# create a class #3 to set-up worker-pools
class _MyPool(mppool.Pool):
    # We sub-class multiprocessing.pool. Pool instead of multiprocessing.Pool
    # because the latter is only a wrapper function, not a proper class.
    ### http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
    ### refer the above link to invoke child processes
    # A class created by AAT
    Process = _NoDaemonProcess
# end of class #3

def _createDepthBelowLandSurfaceCoords1Lev(cube):
    # Dr. Saji / UM_Model_DOC suggested that UM produce Root zone soil model
    # level number is equivalent to 0 to 2m. (i.e. from 1 to 4 layer no)
    
    global __soilFirstSecondFixedSurfaceUnit__ 

    if __soilFirstSecondFixedSurfaceUnit__ == 'cm':
        # So we kept here unit as 'cm'. But points are muliplied by
        # 100 with its  corresponding cm values. Why because, that 
        # 100 will be factorized (divied) in grib_message by setting 
        # scaleFactorOfFirstFixedSurface as 2 and 
        # scaleFactorOfFirstFixedSurface as 2. So we must follow 
        # Lets create new coords with 0, 2m infomation.   
        depth_below_land_surface = iris.coords.DimCoord(numpy.array([15000]), 
                         bounds=numpy.array([[0, 30000]]), units=Unit('cm'),
                                       long_name='depth_below_land_surface')
    elif __soilFirstSecondFixedSurfaceUnit__ == 'mm':  
        # We kept here unit as 'mm'. But points are muliplied by
        # 1000 with its  corresponding cm values. Why because, that 
        # 1000 will be factorized (divied) in grib_message by setting 
        # scaleFactorOfFirstFixedSurface as 3 and 
        # scaleFactorOfSecondFixedSurface as 3. So we must follow 
        # this procedure to get correct results.

        # Lets create new coords with 0, 3m infomation.   
        depth_below_land_surface = iris.coords.DimCoord(numpy.array([1500000]), 
                          bounds=numpy.array([[0, 3000000]]), units=Unit('mm'), 
                                        long_name='depth_below_land_surface')    
    else:
        return 
    # end of if __soilFirstSecondFixedSurfaceUnit__ == 'cm':
        
    # add the above created new coords to the cube 
    cube.add_aux_coord(depth_below_land_surface)    
# end of def _createDepthBelowLandSurfaceCoords1Lev():

def _updateDepthBelowLandSurfaceCoords4Levs(depth_below_land_surface):
    # Dr. Saji / UM_Model_DOC suggested that UM produce soil model
    # level number is equivalent to 10cm, 35cm, 1m & 2m. 
    
    global __soilFirstSecondFixedSurfaceUnit__ 

    if __soilFirstSecondFixedSurfaceUnit__ == 'cm':
        # So we kept here unit as 'cm'. But points are muliplied by
        # 100 with its  corresponding cm values. Why because, that 
        # 100 will be factorized (divied) in grib_message by setting 
        # scaleFactorOfFirstFixedSurface as 2 and 
        # scaleFactorOfFirstFixedSurface as 2. So that in grib2 will
        # be able to read as 0.1 m, 0.35m, 1m & 2m. Iris will convert 
        # cm to m while saving into grib2 file. So we must follow 
        # this procedure to get correct results. (tested and suggested for 
        # WRF-Noah model input)
        
        # 1000 cm -> 10 m -> 10 m / 100 (scaleFactorOfFirstFixedSurface = 2) -> 0.1 m
        # 3500 cm -> 35 m -> 35 m / 100 (scaleFactorOfSecondFixedSurface = 2) -> 0.35 m
        # 10000 cm -> 100 m -> 100 m / 100 (scaleFactorOfSecondFixedSurface = 2) -> 1.0 m
        # 30000 cm -> 300 m -> 300 m / 100 (scaleFactorOfSecondFixedSurface = 2) -> 3.0 m
        
        depth_below_land_surface.points = numpy.array([500, 2250, 6750, 20000])
        # we must set the bounds in vertical depths, since we required
        # to mention the four different layers depth properly.
        depth_below_land_surface.bounds = numpy.array([[0, 1000], 
                                   [1000, 3500], [3500,10000],[10000,30000]])
        depth_below_land_surface.units = Unit('cm')
    elif __soilFirstSecondFixedSurfaceUnit__ == 'mm':
        # Here we kept unit as 'mm'. But points are muliplied by
        # 1000 with its  corresponding mm values. Why because, that 
        # 1000 will be factorized (divied) in grib_message by setting 
        # scaleFactorOfFirstFixedSurface as 3 and 
        # scaleFactorOfSecondFixedSurface as 3. So that in grib2 will
        # be able to read as 0.1m, 0.35m, 1m & 3m. Iris will convert 
        # mm to m while saving into grib2 file. So we must follow 
        # this procedure to get correct results.
        # Moreover IMD-MFI model required to be scaling range of 100000. So we 
        # must follow this procedure only (i.e. mm to m conversion and not cm to m conversion) 
        
        # 100000 mm -> 100 m -> 100 m / 1000 (scaleFactorOfFirstFixedSurface = 3) -> 0.1 m
        # 350000 mm -> 350 m -> 350 m / 1000 (scaleFactorOfSecondFixedSurface = 3) -> 0.35 m
        # 1000000 mm -> 1000 m -> 1000 m / 1000 (scaleFactorOfSecondFixedSurface = 3) -> 1.0 m
        # 3000000 mm -> 3000 m -> 3000 m / 1000 (scaleFactorOfSecondFixedSurface = 3) -> 3.0 m
        
        depth_below_land_surface.points = numpy.array([50000, 225000, 675000, 2000000])
        # we must set the bounds in vertical depths, since we required
        # to mention the four different layers depth properly.
        depth_below_land_surface.bounds = numpy.array([[0, 100000], 
                           [100000, 350000], [350000,1000000],[1000000,3000000]])
        depth_below_land_surface.units = Unit('mm')
    # end of if __soilFirstSecondFixedSurfaceUnit__ == 'cm':
    
    depth_below_land_surface.long_name = 'depth_below_land_surface'    
    depth_below_land_surface.standard_name = 'depth'
# end of def _updateDepthBelowLandSurfaceCoords4Levs():

def _convert2VolumetricMoisture(cube, levels=[100.0, 250.0, 650.0, 2000.0]):
    #### Lets convert moisture_content_of_soil_layer into 
    ##  volumetric_moisture_of_soil_layer by divide each layer 
    ## with its layer depth in mm.
    ## Unit also gets changed from Kg/m2 into m3/m3. How ?
    ## voulumetric_soil_moisture = moisture_content_of_soil_layer / (density_of_water x depth_of_soil_layer)
    ## density_of_water is 1000 Kg/m3
    ## depth_of_soil_layer of first layer from 0 to 10 cm = 10/100 m
    ## depth_of_soil_layer of first layer from 10 to 35 cm = 25/100 m
    ## depth_of_soil_layer of first layer from 35 to 100 cm = 65/100 m
    ## depth_of_soil_layer of first layer from 100 to 300 cm = 200/100 m
    
    ## So if we apply the above values of density 1000 Kg/m3 
    ## & depth in meter in the denominator of voulumetric_soil_moisture
    ## equavation, we will endup with just divide first layer 
    ## by 100, second layer by 250, third layer by 650 and 
    ## fourth layer by 2000.
    
    ## By this way, we converted moisture_content_of_soil_layer 
    ## from Kg/m2 into volumetric_soil_moisture_of_layer m3/m3.
    
    ## Reference : "Comparison of the Met Office soil moisture
    ## analyses with SMOS retrievals (2010-2011)", MARCH 2013.
    
    ## Link : http://www.researchgate.net/publication/257940913
    
    print "before volumetirc", cube.data.min(), cube.data.max()
    if isinstance(levels, (list, tuple)):   
        # This block of code for 4 different layers 
        for idx, dval in enumerate(levels):
            cube.data[idx] /= dval
    elif isinstance(levels, float):
        # this block of code for single layer 
        cube.data /= levels
    # end of for idx, denominator in enumerate([...]):
    print "after volumetirc", cube.data.min(), cube.data.max()
    # WRF-WPS requires minimum vlaue as 0.005. If it is < 0.005 then 
    # Noah thorws segmentation error due to low value of soil moisture. 
    # Reference : look at the lines from 1219 t0 1260 in the below link
    # http://www.cisl.ucar.edu/staff/huangwei/WRFV3/dyn_em/module_initialize_real.F.html
    # though the above Noah code replace the <0.005 grid values with 0.005,
    # but it does only for the first time step (say analysis 00hr), and then 
    # model will blow up for the next time step (say forecast 06hr).
    # And either we should do mask grid points < 0.005 or replace with 0.0051.
    # Here we are replacing with 0.0051 since soil moisture masking will not 
    # make proper sense!. so replace the values less than 0.005 with 0.0051.
    cube.data[numpy.ma.logical_and(cube.data > 0.0, cube.data < 0.005)] = 0.0051
    
    # update the units as m3 / m3
    cube.units = Unit('m3 m-3')
    # make sure that standard_name as None, so that it will  
    # not messup with units while writing as grib2. 
    cube.standard_name = None
    # set long name as volumetric_moisture_of_soil_layer, 
    # though its not standard cf name, I made it as 
    # understandable long name which points into volumetirc 
    # grib2 param code in _grib_cf_map.py.
    cube.long_name = 'volumetric_moisture_of_soil_layer'        
# end of def _convert2VolumetricMoisture(cube):

def _convert2WEASD(cube):
    # http://www.nrcs.usda.gov/wps/portal/nrcs/detail/or/snow/?cid=nrcs142p2_046155
    if cube.standard_name == 'snowfall_amount':
        cube.standard_name = 'liquid_water_content_of_surface_snow'
        # convert the snow amount data to water equivalent by divide by 10 or 
        # multiply by 0.1.
        # reference link : look above        
        cube.data *= 0.1 
    # end of if cube.standard_name == 'snowfall_amount':
# end of def _convert2WEASD(cube):

