import geopandas as gpd
import pandas as pd

#world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world = gpd.read_file('Data/world/TM_WORLD_BORDERS-0.3.shp')
cnames = ['Austria','Sweden','Kenya', 'Japan', 'Brazil', 'Canada', 'Lesotho']
corsiz = [83879, 450295, 580367, 377973, 8514215, 9984670, 30355]
epsgs = [6933, 3395, 6931, 6932, 3035]
proj  = {'Sph   Bonne': '+proj=bonne +lat_1=60 +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs',
         'World Bonne': '+proj=bonne +lat_1=60 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
         'Sph Cassini': '+proj=cass +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs',
         'Wld Cassini': '+proj=cass +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
         'Sph Craster': '+proj=crast +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs',
         'NSIDC EASE ': '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +R=6371228 +units=m +no_defs',
         'NSIDC EASE2': '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs', #epsg 6933
         'LAEA Europe': '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',
         'S Mollweide': '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs',
         'Sinusoidal ': '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs'}

for c in cnames:
    i = cnames.index(c)
    print(c,'\r\t\t\t','wiki:  ','{:8d}'.format(int(corsiz[i])))
    carea = world[world['NAME'] == c]
    for e in epsgs:
        carea = carea.to_crs(epsg=e)
        area = '{:8d}'.format(int(pd.to_numeric(carea['geometry'].area)/10**6))
        perc = abs(1 - (corsiz[i]/float(area)))*100
        print('\t\t','epsg: ',e,':\t',area,'\t off:','{:.2f}'.format(perc))
    for k in proj:
        carea = carea.to_crs(proj[k])
        perc = abs(1 - (corsiz[i]/float(area)))*100
        print('\t\t',k,':\t',area,'\t off:','{:.2f}'.format(perc))
    print('')
