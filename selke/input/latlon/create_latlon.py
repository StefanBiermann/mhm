#!/usr/bin/env python
#############################################
# calculate latitude and longitude coordinates
# given header file specifications
#
# created by Matthias Zink 04.04.2012
# modified Stephan Thober, Nov 2013 - mHM confirmative
# modified Matthias Zink , Feb 2014 - BugFix for wrong latitudes, added command line arguments
#          Stephan Thober, Jun 2015 - added option for performing no coordinate transformation,
#                                     updated import to new ufz module, cs can be float,
#                                     create 1 dimensional xx array now with linspace, which is cleaner and solves a bug in the creation of the y direction
#
#############################################


def header_to_latlon(headerfile, coord_sys):
    # This function returns the latitude longitude given a ASCII header file

    # check input files
    if not os.path.isfile(headerfile):
      sys.exit("No "+headerfile+" file found here, are you in the right directory?")

    # read header information
    header_info = np.loadtxt( headerfile, dtype='|S20')
    ncols       = np.int(header_info[0,1])
    nrows       = np.int(header_info[1,1])
    xllcorner   = np.float(header_info[2,1])
    yllcorner   = np.float(header_info[3,1])
    cs          = np.float(header_info[4,1])
    missVal     = header_info[5,1]

    # create x and y grid
    xx          = np.linspace( xllcorner + cs/2,                xllcorner + cs/2 + (ncols-1)*cs, ncols)
    yy          = np.linspace( yllcorner + cs/2 + (nrows-1)*cs, yllcorner + cs/2,                nrows)
    xx, yy      = np.meshgrid(xx,yy)

    #
    # determine latitude and longitude of the Aimgrid
    if coord_sys == '':
        lons = xx
        lats = yy
    else:
        projAim    = Proj(init=coord_sys)
        lons, lats = projAim(xx, yy, inverse=True)

    return lons, lats, xx, yy, missVal


def latlon_to_nc(fhandle, lons, lats, xx, yy, missVal, suffix):
    # This function writes the latitude longitude to file
    varAtt  = ([['axis'     , 'X']])
    varName = 'xc' + suffix[0]
    dims    = xx.shape[1]
    var     = xx[0,:] #np.arange(ncols)+1
    writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
    #
    varAtt  = ([['axis'     , 'Y']])
    varName = 'yc' + suffix[0]
    dims    = yy.shape[0]
    var     = yy[:,0] #np.arange(nrows)+1
    writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
    # Variables
    # lon
    varAtt  = ([['units'        , 'degrees_east'],
               ['long_name'     , 'longitude' + suffix[1]],
               ['missing_value' , missVal]])
    #
    varName = 'lon' + suffix[0]
    var     = lons
    dims    = ['yc' + suffix[0], 'xc' + suffix[0]]
    writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, vartype='f8', comp=True)
    # lat
    varAtt  = ([['units'        , 'degrees_north'],
               ['long_name'     , 'latitude' + suffix[1]],
               ['missing_value' , missVal]])
    varName = 'lat' + suffix[0]
    var     = lats
    dims    = ['yc' + suffix[0],'xc' + suffix[0]]
    writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, vartype='f8', comp=True)

    
# COORDINATE SYSTEM
#   if empty string is given, no coordinate transformation will be transformed
#   code 31463 defines GK (DHDN3 Zone 3)
#   [website: http://www.spatialreference.org/ref/epsg/31463/]
#   equal +proj=tmerc +ellps=bessel +lon_0=12 +x_0=3,500,000 +y_0=0
# coord_sys = 'epsg:31467' # coordinate system for test basin 2
coord_sys = 'epsg:3035' # coordinate system for test basin 1

# HEADER FILE
#   specifies the grid properties
#   for example, use a copy of the header.txt
#   and adapt cellsize, ncols, nrows to your hydrologic resolution
headerfile_l11 = ''
headerfile_l1  = ''
headerfile_l0  = ''

# OUTPUT FILE
#   path to the output file, latlon.nc is hard-coded in mHM
outfile = 'latlon.nc'

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Cretaes latitude and longitude grids for <coord_sys> with the domain defined in  <headerfile> into NetCDF file <outfile>.")

# usage example with command line arguments
# -------------------------------------------------------------------------
#
# py create_latlon.py -c 'epsg:31463' -f header.txt -o latlon.nc
#
# -------------------------------------------------------------------------

parser.add_option('-c', '--coord_sys', action='store', dest='coord_sys', type='string',
                  default=coord_sys, metavar='Property',
                  help='Coordinate system specifier according to http://www.spatialreference.org. (default: epsg:31467), give empty string for regular latlon grid')
parser.add_option('-f', '--header_l0', action='store', dest='headerfile_l0', type='string',
                  default=headerfile_l0, metavar='Header file for level 0',
                  help='Header file containing information about e.g. number of rows and columns. (default: header.txt).')
parser.add_option('-g', '--header_l1', action='store', dest='headerfile_l1', type='string',
                  default=headerfile_l1, metavar='Header file for level 1',
                  help='Header file containing information about e.g. number of rows and columns at level 1. (default: header.txt).')
parser.add_option('-e', '--header_l11', action='store', dest='headerfile_l11', type='string',
                  default=headerfile_l11, metavar='Header file for level 11',
                  help='Header file containing information about e.g. number of rows and columns at level 1. (default: header.txt).')
parser.add_option('-o', '--outfile', action='store', dest='outfile', type='string',
                  default=outfile, metavar='NetCDF file',
                  help='Name of NetCDF file. (default: latlon.nc).')

(opts, args) = parser.parse_args()

headerfile_l11 = opts.headerfile_l11
headerfile_l1  = opts.headerfile_l1
headerfile_l0  = opts.headerfile_l0
outfile        = opts.outfile   
coord_sys      = opts.coord_sys

# check whether any headerfiles are given
if headerfile_l0 == '' and headerfile_l1 == '' and headerfile_l11 == '':
    raise ValueError('***ERROR: no headerfile specified, use -h switch for more information')

#############################################

import numpy as np                       # array manipulation
import netCDF4 as nc                     # netCDF interphase
import time, os, sys                     # call current time for timestamp
from pyproj import Proj
from ufz    import writenetcdf      # from ufz

# create netCDF Dataset
fName   = outfile
fhandle = nc.Dataset(fName, 'w', format='NETCDF4')
# add file attributes
FiAtt  = ([['description', 'lat lon file'],
           ['projection', coord_sys],
           ['history','Created ' + time.ctime(time.time()) ]])
writenetcdf(fhandle, fileattributes=FiAtt)

if headerfile_l0 != '':
    # get lat lon for level 0 header file
    lons, lats, xx, yy, missVal = header_to_latlon(headerfile_l0, coord_sys)
    suffix = ['_l0', ' at level 0']
    # write lat and lon for level 1 to file
    latlon_to_nc(fhandle, lons, lats, xx, yy, missVal, suffix)
#
if headerfile_l1 != '':
    # get lat lon for level 1 header file
    lons, lats, xx, yy, missVal = header_to_latlon(headerfile_l1, coord_sys)
    suffix = ['', ' at level 1']
    # write lat and lon for level 1 to file
    latlon_to_nc(fhandle, lons, lats, xx, yy, missVal, suffix)
#
if headerfile_l11 != '':
    # get lat lon for level 1 header file
    lons, lats, xx, yy, missVal = header_to_latlon(headerfile_l11, coord_sys)
    suffix = ['_l11', ' at level 11']
    # write lat and lon for level 1 to file
    latlon_to_nc(fhandle, lons, lats, xx, yy, missVal, suffix)

#
# close netcdf dataset
fhandle.close()
print outfile + " created."
