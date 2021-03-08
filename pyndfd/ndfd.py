# Copyright (c) 2015 Marty Sullivan
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''

	NDFD Forecast Retrieval Routines

	Author: 	Marty J. Sullivan and Brandon N. Benton
	Purpose:	Routines that will cache NDFD forecast variables locally
			to allow for easy and fast forecast analysis by lat/lon

'''

###########
#         #
# IMPORTS #
#         #
###########

from datetime import datetime, timedelta
from getpass import getuser
from math import isnan, sqrt
from .ndfd_defs import ndfdDefs
import numpy as np
from numpy.ma.core import MaskedConstant as NAN
from os import makedirs, path
from pyproj import Geod, Proj
from shutil import rmtree
import shutil
from sys import stderr
from tempfile import gettempdir
from urllib.request import urlretrieve, urlopen
import json
import cfgrib
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
i
#############
#           #
# CONSTANTS #
#           #
#############

DEFS = ndfdDefs()
G = Geod(ellps='clrk66')

CACHE_SERVER_BUFFER_MIN = 15

NDFD_LOCAL_SERVER = None
RTMA_LOCAL_SERVER = None
NDFD_REMOTE_SERVER = 'http://tgftp.nws.noaa.gov/SL.us008001/ST.opnl/DF.gr2/'
RTMA_REMOTE_SERVER = NDFD_REMOTE_SERVER

########################
#                      #
# FUNCTION DEFINITIONS #
#                      #
########################

'''
  Function:	setLocalCacheServer
  Purpose:	Set a server to use instead of weather.noaa.gov
  Params:
	uri:	String denoting the server URI to use
'''

def setLocalCacheServerNDFD(uri):
    global NDFD_LOCAL_SERVER 
    NDFD_LOCAL_SERVER = uri


'''
  Function: 	getLatestForecastTime
  Purpose:  	For caching purposes, compare this time to cached time to see if
		the cached variable needs to be updated
'''

def getLatestForecastTime():
    latestTime = datetime.utcnow()
    if latestTime.minute <= CACHE_SERVER_BUFFER_MIN:
        latestTime = (datetime.utcnow() - timedelta(hours=1))
    return latestTime.replace(minute=0, second=0, microsecond=0)


'''
  Function:	getElevationVariable
  Purpose:	Cache the static elevation variable if not already 
                cached and return the path of the cached file
  Params:
	area:	The NDFD grid area to retrieve elevation for
  
  Notes:
	- Cannot be retrieved from weather.noaa.gov, must use a local cache server
	  using the format in const NDFD_STATIC
	- Puerto Rico terrian info not currently available. 
	- Terrain data for NDFD will be updated sometime in 2015
'''

def getElevationVariable(area):
    if area == 'puertori':
        raise ValueError('Elevation currently not available for Puerto Rico. Set elev=False')
    if NDFD_LOCAL_SERVER == None:
        raise RuntimeError('Local cache server must provide elevation data. Specify cache server with ndfd.setLocalCacheServer(uri)')
    if not path.isdir(NDFD_TMP):
        makedirs(NDFD_TMP)
    remoteVar = NDFD_LOCAL_SERVER + NDFD_STATIC.format(area) + NDFD_VAR.format('elev')
    localDir = NDFD_TMP + NDFD_STATIC.format(area)
    localVar = localDir + NDFD_VAR.format('elev')
    if not path.isdir(localDir):
        makedirs(localDir)
    if not path.isfile(localVar):
        urlretrieve(remoteVar, localVar)
    if not path.isfile(localVar):
        raise RuntimeError('Cannot retrieve NDFD variables at this time. Try again in a moment.')
    return localVar

'''
  Function:	getSmallestGrid
  Purpose:	Use the provided lat, lon coordinates to find the smallest
		NDFD area that contains those coordinates. Return the name of the area.
  Params:
	lat:	Latitude 
	lon:	Longitude
'''

def getSmallestGrid(lat, lon):
    smallest = 'neast'
    minDist = G.inv(lon, lat, DEFS['grids'][smallest]['lonC'], DEFS['grids'][smallest]['latC'])[-1]
    for area in DEFS['grids'].keys():
        if area == 'conus' or area == 'nhemi' or area == 'npacocn':
            continue
        curArea = DEFS['grids'][area]
        smallArea = DEFS['grids'][smallest]
        dist = G.inv(lon, lat, curArea['lonC'], curArea['latC'])[-1]
        if dist < minDist:
            minDist = dist
            smallest = area

    return smallest

'''
  Function:	getNearestXrGridPoint
  Purpose:	Find the nearest grid point to the provided coordinates. 
                Return the indexes to the numpy array as well as the lat/lon 
                and grid coordinates of the grid point.
  Params:
	ds:		    xarray data for grib file
	lat:		Latitude
	lon:		Longitude
'''

def getNearestXrGridPoint(ds,lat,lon):
    lats = ds.variables['latitude']
    lons = ds.variables['longitude']
    lat_arr = np.array(lats.data)-np.float64(lat)
    lon_arr = np.array(lons.data)-np.float64(lon)
    '''
    for i in range(lat_arr.shape[0]):
        for j in range(lat_arr.shape[1]):
            lat_arr[i,j]-=np.float64(lat)
            lon_arr[i,j]-=np.float64(lon)
    '''
    dist = (lat_arr)**2+(lon_arr)**2
    idy, idx = np.where(dist == dist.min())
    gLat = lats.data[idy, idx]
    gLon = lons.data[idy, idx]
    return idx, idy, gLat, gLon      

'''
  Function:	validateArguments
  Purpose:	Validate the arguments passed into an analysis function to make sure
		they will work with each other.
  Params:
	var:		The NDFD variable being requested
	area:		The NDFD grid area being requested
	timeStep:	The time step to be used in the returned analysis
	minTime:	The minimum forecast time to analyze
	maxTime:	The maximum forecast time to analyze
  Notes:
	- maxTime is not currently being evaluated
'''

def validateArguments(var, area, timeStep, minTime, maxTime):
    if timeStep < 1:
        raise ValueError('timeStep must be >= 1')
    if minTime != None and minTime < getLatestForecastTime():
        raise ValueError('minTime is before the current forecast time.')
    if maxTime != None and maxTime < getLatestForecastTime():
        raise ValueError('maxTime is before the current forecast time.')
 
    try:
        areaVP = DEFS['vars'][area]
    except IndexError:
        raise ValueError('Invalid Area.')

    validVar = False
    for vp in areaVP:
        if var in areaVP[vp]:
            validVar = True
            break
    if not validVar:
        raise ValueError('Variable not available in area: ' + area)

'''
  Function:	getLocationData
  Purpose:	To get raw data at a specific location for 
                current time and ndfd forecast period
'''
    
def getLocationData(ds_path, var, lat, lon, timeStep=1, elev=False, area=None):
    
    if area == None:
        area = getSmallestGrid(lat, lon)
    validateArguments(var, area, timeStep, minTime, maxTime)

    analysis = {}
    analysis['forecastTime'] = getLatestForecastTime()
    analysis['variables'] = { }
    
    ds = xr.open_dataset(ds_path,engine='cfgrib')
    times = ds.valid_time.values.astype('datetime64[s]').tolist()
    x,y,gLat,gLon = getNearestXrGridPoint(ds,lat,lon)
    
    tmp=ds.isel(y=y,x=x)

    for variable in tmp.data_vars:
        varData = {}
        for i in range(tmp.dims['step']):
            time = times[i]
            value = tmp.variables[variable][i,0,0].values
            varData[time.strftime("%Y-%m-%d:%H")] = float(value)
    
        analysis['variables'][variable] = varData
    '''
    for g in varRTMA:
        
        ds=xr.open_dataset(g,engine='cfgrib')
        x,y,gLat,gLon=getNearestXrGridPoint(ds,lat,lon)
        
        tmp = ds.isel(y=y,x=x)

        for variable in tmp.data_vars:
            
            time = tmp.valid_time.values.astype('datetime64[s]').tolist()
            value = tmp.variables[variable][0,0].values
            #value=tmp[0,0].values
            analysis['variables'][variable] = {time.strftime("%Y-%m-%d:%H"): float(value),**analysis['variables'][variable]}
    
    '''

    return analysis

'''
  Function:	plotData
  Purpose:      Plot data from getLocationData
'''

def plotData(data,var):

    xtmp = [mdates.date2num(datetime.strptime(x,'%Y-%m-%d:%H')) for x in data[var].keys()]
    x = sorted(xtmp)
    ytmp = [y for y in data[var].values()]
    y = []

    #hack for removing weird temp values
    for i in range(len(ytmp)):
        if ytmp[i]<0:
            ytmp[i] = (ytmp[i-1]+ytmp[i+1])/2

    for i in range(len(ytmp)):
        idx = xtmp.index(x[i])
        y.append(ytmp[idx])
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d:%H'))

    for i in range(len(y)):
        if y[i] == ytmp[0]: color = 'red'
        else: color = 'blue'

        plt.plot_date(x[i], y[i], color = color)
    
    plt.plot(x,y)
    plt.ylim([min(y)-np.std(y), max(y)+np.std(y)])
    plt.xticks(rotation=90)
    fig.autofmt_xdate()

    fig_name = NDFD_TMP + path.sep + '%s.png' %(var)
    fig.savefig(fig_name)

    print('Saved %s' %(fig_name))
