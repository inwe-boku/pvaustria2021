import urllib.parse
import requests
from io import StringIO
import pandas as pd
import math
import numpy
from scipy.spatial import cKDTree
import geopandas as gpd
import yaml
import click
import os
import pathlib
from osgeo import gdal
import struct
import osmnx as ox
import itertools

def nearestbld(bds, buildings):
    nA = numpy.array(list(zip(bds['geometry'].centroid.x, bds['geometry'].centroid.y)) )
    nB = numpy.array(list(zip(buildings['geometry'].centroid.x, buildings['geometry'].centroid.y)) )
    btree = cKDTree(nA)
    dist, idx = btree.query(nB,k=1)
    buildings['nearest'] = idx
    return(buildings)

def haversine(lat1, lon1, lat2, lon2):  # returns the distance via the haversine formula (km)
    R = 6371
    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.cos(math.radians(lat1)) * \
                 math.cos(math.radians(lat2)) * math.sin(dLon / 2) * math.sin(dLon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = R * c
    return(d)

def km2deg(bbox):
  lat1 = (bbox[0]+bbox[2])/2
  lon1 = (bbox[1]+bbox[3])/2

  lat2 = lat1 + 1
  lon2 = lon1 + 0
  dlat = haversine(lat1, lon1, lat2, lon2)

  lat2 = lat1 + 0
  lon2 = lon1 + 1
  dlon = haversine(lat1, lon1, lat2, lon2)

  return(1/dlat, 1/dlon)

def nlatlon(resolution, bbox):
    reslat, reslon = km2deg(bbox)
    return(reslat*resolution, reslon*resolution)

def writeGEO(data, path, dataname):
    #data.to_file(filename = os.path.join(path, 'geojson', dataname+'.geojson'), driver="GeoJSON")
    #data.to_file(filename = os.path.join(path, 'shape', dataname+'.shp'), driver = 'ESRI Shapefile')
    data.to_file(filename = os.path.join(path, 'data.gpkg'), layer = dataname, driver = 'GPKG')
    return(0)

def getPVJRC(url, timeres):
    r = requests.get(url)
    data = StringIO(r.text)
    if timeres == 'hourly':
        dateparse = lambda x: pd.datetime.strptime(x, '%Y%m%d:%H%M')
        df = pd.read_csv(data, header = 1, skiprows=9, skipfooter=10, parse_dates=['Date'], date_parser=dateparse, engine='python')
        df.assign(Date=df.Date.dt.round('H'))
        df=df.rename(columns = {'EPV':'Wh'})
    elif timeres == 'monthly':
        dateparse = lambda x: pd.datetime.strptime(x, '%m')
        df = pd.read_csv(data, sep='\t\t', header = 1, parse_dates=['Month'], date_parser=dateparse, skiprows=8, skipfooter=12, engine = 'python')
        df=df.rename(columns = {'Month':'Date'})
        df=df.rename(columns = {'E_m':'Wh'})
    return(df)

def getshapes():
    buildings = gpd.read_file("shapefile.shp")

def getyml(yamlfile):
    with open(yamlfile, 'r') as stream:
        try:
            yml = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return(yml)

def pointraster(resolution, area):
    bbox=[area.bounds.miny.min(), area.bounds.minx.min(), area.bounds.maxy.max(), area.bounds.maxx.max()]
    nlat, nlon = nlatlon(resolution, bbox)
    lats = list(numpy.arange(bbox[0], bbox[2], nlat))
    lons = list(numpy.arange(bbox[1], bbox[3], nlon))
    llats = len(lats)
    llons = len(lons)

    if dbg:
        print('resolution: {:.2f} km;\nlatitude degrees: {:.5f} - longitude degrees: {:.5f}\ngridsize: {:d}x{:d}'.format(resolution, nlat, nlon, llats, llons))
        print('bounding box:',bbox)

    lats = numpy.repeat(lats,llons)
    lons = numpy.tile(lons,llats)

    print('Sampling points in bounding box: ', llats*llons)

    df = pd.DataFrame({'Latitude':lats, 'Longitude':lons})
    points = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    points.crs = 'epsg:4326'

    points = gpd.sjoin(points, area, how='left')
    points = points.dropna()
    points = points.reset_index(drop=True)
    points = points[['geometry', 'GKZ', 'KG_NR']]
    points['type'] = 'none'

    if dbg:
        print(points.head())

    print('Sampling points in area: ', len(points))
    return(points)

def checkrasterint(points, resolution, rasterfile, maxval, checktype):
    ds = gdal.Open(rasterfile)
    gt=ds.GetGeoTransform()
    rb=ds.GetRasterBand(1)
    bbox=[points.bounds.miny.min(), points.bounds.minx.min(), points.bounds.maxy.max(), points.bounds.maxx.max()]
    nlat, nlon = nlatlon(resolution, bbox)
    values = [0,1/5,-1/5,1/3,-1/3]
    tups = list(itertools.product(values, values))
    for idx, row in points.iterrows():
        print('\tsampling ',checktype,' {:d} of {:d}'.format(idx+1, len(points)), end='\r')
        if (row.type != 'none'):
            continue
        intlist = []
        for tup in tups:
            py = int((row.geometry.y + tup[0]*nlat - gt[3]) / gt[5]) #y pixel
            px = int((row.geometry.x + tup[1]*nlon - gt[0]) / gt[1]) #x pixel

            intval = int(struct.unpack('h' , rb.ReadRaster(px,py,1,1,buf_type=gdal.GDT_UInt16))[0])
            intlist.append(intval)
        intval = int(sum(intlist)/len(intlist))
        if (intval >= maxval):
            points.at[idx, 'type'] = checktype
    return(points)

def getfreelandpoints(points, resolution, landusefile, luvalues):
    ds = gdal.Open(landusefile)
    gt=ds.GetGeoTransform()
    rb=ds.GetRasterBand(1)
    bbox=[points.bounds.miny.min(), points.bounds.minx.min(), points.bounds.maxy.max(), points.bounds.maxx.max()]
    nlat, nlon = nlatlon(resolution, bbox)
    values = [0,1/5,-1/5,1/3,-1/3]
    tups = list(itertools.product(values, values))
    if dbg:
        print(gt)

    for idx, row in points.iterrows():
        print('\tsampling LU {:d} of {:d}'.format(idx+1, len(points)), end='\r')
        if (row.type != 'none'):
            continue
        intlist = []
        for tup in tups:
            py = int((row.geometry.y + tup[0]*nlat - gt[3]) / gt[5]) #y pixel
            px = int((row.geometry.x + tup[1]*nlon - gt[0]) / gt[1]) #x pixel

            intval = int(struct.unpack('h' , rb.ReadRaster(px,py,1,1,buf_type=gdal.GDT_UInt16))[0])
            if (intval in luvalues):
                intlist.append(1)
            else:
                intlist.append(0)
        value = int(sum(intlist) / len(intlist) * 100)
        if (value > 40):
            points.at[idx, 'type'] = 'free'
    return(points)

def getosmbuildings(area, fbuildings):
    buildings = gpd.GeoDataFrame()
    bdsample = gpd.GeoDataFrame(columns=['geometry', 'KG_NR', 'fparea'])
    for idx, row in area.iterrows(): # get building area per municipality
        print('\tOSM fetch {:d} of {:d}'.format(idx+1, len(area)))#, end='\r')
        if dbg:
            print("\t\t", row['KG'], " - ", row['BL'])
        if fbuildings.empty == False:
            rgdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(row[['geometry']]), crs = 'epsg:4326')
            fbuildings.to_crs('epsg:3395')
            rgdf.to_crs('epsg:3395')
            bd = gpd.clip(fbuildings, rgdf)
            #bd = gpd.sjoin(fbuildings, rgdf, how="inner", op='intersects')
            #bd = bd.dropna()
            bd = bd.reset_index(drop=True)
            bd.crs = rgdf.crs
        else:
            bd = ox.footprints_from_polygon(row.geometry)
        if (len(bd)) == 0:
            continue
        bd = bd[['geometry']]
        bd = bd.to_crs('epsg:6933')
        bd['fparea'] = bd['geometry'].area.astype(int)
        bd = bd.to_crs('epsg:3395')
        bd['KG_NR'] = row['KG_NR']
        ssize = min(len(bd), round(3*len(bd)**0.40)) #int(len(bd)/math.sqrt(len(bd))
        print('\t\tbuildings {:d}, sample {:d}, total {:d}'.format(len(bd), ssize, len(fbuildings)))
        bds = bd.sample(ssize)
        bd = nearestbld(bds, bd)
        bdj = bd.groupby([bd.nearest]).sum()
        bds['index'] = bds.reset_index().index
        bds = bds.set_index('index').join(bdj,rsuffix='r')
        bds = bds[['geometry', 'KG_NR', 'fparear']]
        bds.rename(columns={'fparear': 'fparea'}, inplace=True)
        bdsample = pd.concat([bdsample, bds])
        buildings = pd.concat([buildings, bd])

    buildings['KG_NR']=pd.to_numeric(buildings.KG_NR)
    buildings['fparea']=pd.to_numeric(buildings.fparea)
    bdsample['KG_NR']=pd.to_numeric(bdsample.KG_NR)
    bdsample['fparea']=pd.to_numeric(bdsample.fparea)
    bdsample.crs = {'init' :'epsg:3395'}
    bdsample = bdsample.reset_index(drop=True)
    return(buildings, bdsample)

def calcfreepv(pointspv, urloptions, urlbase, seasons, store, timeres):
    count = 0
    ln = len(pointspv[pointspv['type'] == 'free'])

    for idx, row in pointspv.iterrows():
        if row.type != 'free':
            continue
        count += 1
        print('\tPV data {:d} of {:d}'.format(count, ln), end='\r')
        if dbg:
            print('')
        lonlat = {'lat': '{:.4f}'.format(row['geometry'].y),
                  'lon': '{:.4f}'.format(row['geometry'].x)}
        if dbg:
            print(lonlat)
        urloptions['mountingplace'] = 'free'
        url = urlbase + urllib.parse.urlencode({**lonlat, **urloptions}).strip("'")
        dataname = 'C'+'x'.join("{!s}".format(str(v).replace('.', '_')) for (k,v) in lonlat.items())
        if dbg:
            print(dataname)
        try:
            pvdata = store[dataname]
            if dbg:
                print('existing, loaded from file')
        except:
            if dbg:
                print('not here, downloading from JRC')
            pvdata = getPVJRC(url, timeres)
            store[dataname] = pvdata

        pvdata['year'] = pvdata['Date'].dt.year
        pvdata['month'] = pvdata['Date'].dt.month
        pvdata['day'] = pvdata['Date'].dt.day

        df = pvdata.groupby([pvdata.month]).sum()
        nryears = len(pvdata['Date'].dt.year.unique())
        for s in seasons:
            if timeres == 'daily':
                kWhpkWp = int(df.loc[seasons[s]].Wh.mean()/1000*12/nryears)
            elif timeres == 'monthly':
                kWhpkWp = int(df.loc[seasons[s]].Wh.mean()*12)
            if dbg:
                print(s, ': ', kWhpkWp, 'kWh/a per kWp')
            fname1 = 'kWhp_'+ s
            fname2 = 'kWhm2_'+ s
            if s == 'all':
                fname1 = 'kWhp'
                fname2 = 'kWhm2'
            pointspv.at[idx, fname1] = kWhpkWp
            pointspv.at[idx, fname2] = int(kWhpkWp/cfg['landuse']['m2kWp'])
    return(pointspv)

def calcbldpv(bdsample, area, urloptions, urlbase, seasons, store, timeres):
    count = 0
    ln = len(bdsample)
    bdsample = bdsample.to_crs('epsg:4326')

    for idx, row in bdsample.iterrows():
        count += 1
        print('\tPV data {:d} of {:d}'.format(count, ln), end='\r')
        if dbg:
            print()
            print()
        lonlat = {'lat': '{:.4f}'.format(row['geometry'].centroid.y),
                  'lon': '{:.4f}'.format(row['geometry'].centroid.x)}
        if dbg:
            print(lonlat)
        urloptions['mountingplace'] = 'building'
        url = urlbase + urllib.parse.urlencode({**lonlat, **urloptions}).strip("'")

        dataname = 'C'+'x'.join("{!s}".format(str(v).replace('.', '_')) for (k,v) in lonlat.items())
        if dbg:
            print(dataname)
        try:
            pvdata = store[dataname]
            if dbg:
                print('existing, loaded from file')
            try:
                df=df.rename(columns = {'E_m':'Wh'})
            except:
                if dbg:
                    print('renamed again')

        except:
            if dbg:
                print('not here, downloading from JRC')
            pvdata = getPVJRC(url, timeres)
            print(url)
            print(pvdata)
            store[dataname] = pvdata

        pvdata['year'] = pvdata['Date'].dt.year
        pvdata['month'] = pvdata['Date'].dt.month
        pvdata['day'] = pvdata['Date'].dt.day
        df = pvdata.groupby([pvdata.month]).sum()
        nryears = len(pvdata['Date'].dt.year.unique())
        for s in seasons:
            if timeres == 'daily':
                kWhpkWp = int(df.loc[seasons[s]].Wh.mean()/1000*12/nryears)
            elif timeres == 'monthly':
                kWhpkWp = int(df.loc[seasons[s]].Wh.mean()*12)
            if dbg:
                print(s, ': ', kWhpkWp, 'kWh/a per kWp')
            fname1 = 'kWhp_' + s
            fname2 = 'kWhm2_' + s
            fname3 = 'kWh_' + s
            if s == 'all':
                fname1 = 'kWhp'
                fname2 = 'kWhm2'
                fname3 = 'kWh'
            bdsample.at[idx, fname1] = kWhpkWp
            bdsample.at[idx, fname2] = int(kWhpkWp/cfg['buildings']['m2kWp'])
            bdsample.at[idx, fname3] = int(kWhpkWp/cfg['buildings']['m2kWp']*row['fparea']*cfg['buildings']['rPV'])
    bds = bdsample.filter(regex='KG_NR|fparea|kWh$|kWh_')
    dfs = bds.groupby([bds.KG_NR]).sum()
    bds = bdsample.filter(regex='KG_NR|kWhp|kWhm2')
    dfm = bds.groupby([bds.KG_NR]).mean()
    areabuildingspv = area
    areabuildingspv['KG_NR']=pd.to_numeric(areabuildingspv.KG_NR)
    areabuildingspv = areabuildingspv.set_index('KG_NR').join(dfs, rsuffix='rs')
    areabuildingspv = areabuildingspv.join(dfm, rsuffix='rm')
    return(bdsample, areabuildingspv)

@click.command()
@click.option('-resolution', '-r', help = 'output resolution in km x km', default = 0.25, type=float)
@click.option('-areaname', '-a', help = 'specify the name of the area', default = 'nockberge', type = str)
@click.option('-inpdir', '-i', help = 'input directory', default = 'Data/input', type = str)
@click.option('-outdir', '-o', help = 'output directory', default = 'Data/output', type = str)
@click.option('-pvddir', '-p', help = 'raw pv data directory', default = 'Data/pvdata', type = str)
@click.option('-landusefile', '-l', help = 'land use file (optional)', default = '', type = str)
@click.option('-elevationfile', '-e', help = 'elevation file (optional)', default = '', type = str)
@click.option('-osmbld', '-b', help = 'buildings from OSM (optional)', default = False, type = bool)
@click.option('-configfile', '-c', help = 'configuration file', default = 'retour.yml', type = str)
@click.option('-jrcconfig', '-j', help = 'jrc download config file', default = 'pvjrc-monthly.yml', type = str)
@click.option('-debug', '-d', help = 'debug output', default = False, type=bool)
def main(resolution, areaname, inpdir, outdir, pvddir, landusefile, elevationfile, osmbld, configfile, jrcconfig, debug):
    pd.options.mode.chained_assignment = None
    config = getyml(os.path.join(os.path.dirname(os.path.realpath(__file__)), configfile))
    jrccfg = getyml(os.path.join(os.path.dirname(os.path.realpath(__file__)), jrcconfig))
    areafile = os.path.join(inpdir, areaname, 'area.shp')
    area = gpd.read_file(areafile)

    global dbg
    global cfg
    global jrc
    dbg = debug
    cfg = config
    jrc = jrccfg
    timeres = str(jrc['timeres'])
    outdirarea = os.path.join(outdir, areaname)

    if dbg:
        print(areafile)
    pathlib.Path(os.path.join(outdir, areaname, 'geojson')).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.join(outdir, areaname, 'shape')).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.join(pvddir)).mkdir(parents=True, exist_ok=True)
    #writeGEO(area, outdirarea, 'dbg_area')

    if dbg:
        print(jrc['options'])

    points = pointraster(resolution, area)
    writeGEO(points, outdirarea, 'dbg_points')

    if landusefile:
        if elevationfile:

            try:
                maxalt = int(cfg['landuse']['maxalt'])
            except:
                maxalt = 0
            if (maxalt > 0):
                print('')
                print('Processing Points Altitude: <', maxalt, 'masl')
                points = checkrasterint(points, resolution, elevationfile, maxalt, 'altitude')

            try:
                maxslope = int(cfg['landuse']['maxslope'])
            except:
                maxslope = 0
            if (maxslope > 0):
                print('')
                print('Processing Points Slope: <', maxslope, ' degrees')
                opts = gdal.DEMProcessingOptions(scale=111120)
                slopefile = '/tmp/slope.tif'
                gdal.DEMProcessing(slopefile, elevationfile, 'slope', options=opts)
                points = checkrasterint(points, resolution, slopefile, maxslope, 'slope')

            writeGEO(points, outdirarea, 'dbg_pointsaltslp')
        print('')
        print('Processing Points LU')
        pointslu = getfreelandpoints(points, resolution, landusefile, cfg['landuse']['free'])
        writeGEO(pointslu, outdirarea, 'dbg_pointslu')
        print('')
        store = pd.HDFStore(os.path.join(pvddir, 'free_'+
                                         str(jrc['options']['raddatabase'])+'_'+
                                         str(jrc['options']['startyear'])+'_'+
                                         str(jrc['options']['endyear'])+'_'+
                                         timeres+'.hdf5'))
        freepvpts = calcfreepv(pointslu, jrc['options'], jrc['url'], cfg['seasons'], store, timeres)
        store.close()
        writeGEO(freepvpts, outdirarea, 'freepvpts')


    if osmbld:
        print('')
        print('Processing buildings')
        bldfile = os.path.join(inpdir, areaname, 'buildings.shp')
        try:
            buildings = gpd.read_file(bldfile)
            fromfile = True
        except:
            buildings = gpd.GeoDataFrame()
            fromfile = False
        buildings, bdsample = getosmbuildings(area, buildings)
        print('')
        writeGEO(buildings, outdirarea, 'dbg_buildings')
        store = pd.HDFStore(os.path.join(outdirarea, 'building_'+
                                         str(jrc['options']['raddatabase'])+'_'+
                                         str(jrc['options']['startyear'])+'_'+
                                         str(jrc['options']['endyear'])+'_'+
                                         timeres+'.hdf5'))
        bdsample, areabuildingspv = calcbldpv(bdsample, area, jrc['options'], jrc['url'], cfg['seasons'], store, timeres)
        store.close()
        writeGEO(bdsample, outdirarea, 'dbg_bdsample')
        writeGEO(areabuildingspv, outdirarea, 'areabuildingspv')
        print('')
    print('done')

if __name__ == "__main__":
    main()
