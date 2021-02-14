#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 12:05:33 2020

@author: coralie
"""

import numpy as np
import iris
import iris.coord_categorisation
import iris
from iris.coords import AuxCoord, DimCoord
from osgeo import gdal
import sys


def read_geotiff(f, proj=None, rasterband=1, **cubekwargs):

    ''' Read a geotiff file, generate lat/lon coordinates from the metadata, and output a cube.
    Note it assumes the GeoTransform outputs lat/lon and not metres, 
    if not input your own proj object to execute the transformation '''

    ds = gdal.Open(f)
    band = ds.GetRasterBand(rasterband)
    data = band.ReadAsArray()

    # get coordinates
    ny, nx = data.shape
    gt = ds.GetGeoTransform()
    
    # Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)   Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
    
    # In case of north up images, the GT(2) and GT(4) coefficients are zero, and the GT(1) is pixel width, and
    # GT(5) is pixel height.  The (GT(0),GT(3)) position is the top left corner of the top left pixel of the raster.

    # so lon uses the value of nx to create an array with nx number of cells. each cell is then multiplied by the pixel 
    # width and then the x position of the top left pixel in the top left corner 

    lon = np.arange(nx)*gt[1]+gt[0]
    
    # lon does the same by using the value for ny and uses pixel height and the y position of the top left pixel

    lat = np.arange(ny)*gt[5]+gt[3]

    # if proj is provided it means that 'lat' and 'lon' are actually metres
    # that need converting, producing 2D lat/lon arrays
    if proj:
        (mx, my) = np.meshgrid(lon, lat)
        (lon, lat) = proj(mx, my, inverse=True)

        cube = iris.cube.Cube(data, **cubekwargs)

        longitude = AuxCoord(lon, standard_name='longitude', units='degrees')
        latitude = AuxCoord(lat, standard_name='latitude', units='degrees')
        cube.add_aux_coord(longitude, data_dims=[0,1])
        cube.add_aux_coord(latitude, data_dims=[0,1])

    else: # lat/lon are 1D
 
        longitude = DimCoord(lon, standard_name='longitude',units='degrees')
        latitude = DimCoord(lat, standard_name='latitude',units='degrees')
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(latitude,0),(longitude,1)],**cubekwargs)
    
    return cube


#filepath = '/home/coralie/bash_project/Africa/'
#latmin, latmax, lonmin, lonmax = 5, 6, 15, 16

filepath = sys.argv[1]
filename = sys.argv[2]
lat, lon = float(sys.argv[3]), float(sys.argv[4])

# filepath = '/home/coralie/bash_project/Africa/'
# filename = 'Hansen_GFC-2019-v1.7_lossyear_10N_000E.tif'
# lat, lon = 6, 6
latmin, latmax, lonmin, lonmax = lat - 0.5, lat + 0.5, lon - 0.5, lon + 0.5

StandardNomenclature = str(lonmin) + '-' + str(lonmax) + ' lat' + str(latmin) + '-' + str(latmax)

ForestLoss = read_geotiff(filepath + filename, long_name='Year of Gross Forest Cover Loss Event', var_name='yearLoss')

constraint_lon = iris.Constraint(longitude = lambda cell: lonmin <= cell < lonmax)    
constraint_lat = iris.Constraint(latitude = lambda cell: latmin <= cell < latmax )

ForestLoss_constrained = ForestLoss.extract(constraint_lon & constraint_lat)      

iris.save(ForestLoss_constrained, filepath + 'ForestLoss2000-19 lon' + StandardNomenclature + '.nc')
