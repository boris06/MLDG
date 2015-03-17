from numpy import *
from netCDF4 import Dataset
from shapely.geometry import *
from shapely.ops import *
import fiona
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pysal.cg as cg
from datetime import datetime, timedelta
import orange, orngAssoc, string
import re
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as ptch
from matplotlib.collections import PatchCollection
import ConfigParser
import logging
import time
import igraph
import sys
import fileinput

def load_meshmask(meshfile):

    f = Dataset(str(meshfile).replace("/","\\"), 'r', format='NETCDF3_CLASSIC')    

    nav_lon = f.variables['nav_lon'][:,:]
    nav_lat = f.variables['nav_lat'][:,:]
    tmask = f.variables['tmask'][:,:]

    navLonTran = nav_lon.transpose()
    navLatTran = nav_lat.transpose()
    tmaskTran = tmask[0,0,:,:].transpose()
    im = navLonTran.shape[0]
    jm = navLonTran.shape[1]                 
    
    f.close()

    F_lon_min = navLonTran.min()
    F_lon_max = navLonTran.max()
    F_lat_min = navLatTran.min()
    F_lat_max = navLatTran.max()

    # correction #
    F_lon_min = navLonTran.min() - 1./16./2.
    F_lon_max = navLonTran.max() + 1./16./2.
    F_lat_min = navLatTran.min() - 1./16./2.
    F_lat_max = navLatTran.max() + 1./16./2.

    return (navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max)

def write_config(configfile,meshfile):

    print meshfile    

    (navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    config = ConfigParser.RawConfigParser()

    config.add_section('GeoFrame')
    config.set('GeoFrame', 'F_lon_min', F_lon_min)
    config.set('GeoFrame', 'F_lon_max', F_lon_max)
    config.set('GeoFrame', 'F_lat_min', F_lat_min)
    config.set('GeoFrame', 'F_lat_max', F_lat_max)

    with open(configfile, 'wb') as configfile:
        config.write(configfile)

    return    

def generate_mask(meshfile):

    (navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    dlon = zeros((im, jm))
    dlat = zeros((im, jm))

    for i in range(im-1):
        for j in range(jm-1):
            dlon[i,j] = navLonTran[i+1,j] - navLonTran[i,j]
            dlat[i,j] = navLatTran[i,j+1] - navLatTran[i,j]

    dlon[im-1,:] = dlon[im-2,:]            
    dlat[im-1,:] = dlat[im-2,:]

    dlon[:,jm-1] = dlon[:,jm-2]   
    dlat[:,jm-1] = dlat[:,jm-2]

    dlon[im-1,jm-1] = dlon[im-2,jm-2]    
    dlat[im-1,jm-1] = dlat[im-2,jm-2]       
 
    cells = []
    for i in range(im):
        for j in range(jm):
            if (tmaskTran[i,j] == 1):
                cells.append(Polygon([(navLonTran[i,j]-dlon[i,j]/2,navLatTran[i,j]-dlat[i,j]/2), \
                                     (navLonTran[i,j]-dlon[i,j]/2,navLatTran[i,j]+dlat[i,j]/2), \
                                     (navLonTran[i,j]+dlon[i,j]/2,navLatTran[i,j]+dlat[i,j]/2), \
                                     (navLonTran[i,j]+dlon[i,j]/2,navLatTran[i,j]-dlat[i,j]/2), \
                                     (navLonTran[i,j]-dlon[i,j]/2,navLatTran[i,j]-dlat[i,j]/2)]))
    mask = cascaded_union(cells)

    return mask

def show_mask(mask,configfile,figfile):

    #(navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    config = ConfigParser.ConfigParser()
    config.read(configfile)
    F_lon_min = config.get('GeoFrame', 'F_lon_min')
    F_lon_max = config.get('GeoFrame', 'F_lon_max')
    F_lat_min = config.get('GeoFrame', 'F_lat_min')
    F_lat_max = config.get('GeoFrame', 'F_lat_max')
      
    map = Basemap(llcrnrlon=F_lon_min,llcrnrlat=F_lat_min,urcrnrlon=F_lon_max,urcrnrlat=F_lat_max,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='white',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    
    nmask = len(mask.geoms)    
    for i in range(nmask):
        if (mask.geoms[i].geom_type == 'Polygon'):
            lon_wmc = mask.geoms[i].centroid.x
            lat_wmc = mask.geoms[i].centroid.y
            x = []
            y = []
            polygon = list(mask.geoms[i].exterior.coords)
            lp = len(polygon)
            for j in range(lp):
                x.append(polygon[j][0])
                y.append(polygon[j][1])
            (xp,yp)=map(x,y)
            plt.fill(xp,yp,'blue',alpha=0.5)
            (xt,yt)=map(lon_wmc,lat_wmc)
            plt.text(xt,yt,`i+1`,fontsize=24,fontweight='normal',ha='center',va='center',color='white')
            nisl = len(mask.geoms[i].interiors)
            if (nisl > 0):
                for k in range(nisl):
                    island = list(mask.geoms[i].interiors[k].coords)
                    x = []
                    y = []
                    li = len(island)
                    for j in range(li):
                        x.append(island[j][0])
                        y.append(island[j][1])
                    (xp,yp)=map(x,y)
                    plt.fill(xp,yp,'red',alpha=0.5)
        else:
            print mask.geoms[i].geom_type
            
    plt.savefig(figfile,dpi=plt.gcf().dpi,format='pdf')            
    plt.savefig(str(figfile).replace(".pdf",".png"),dpi=plt.gcf().dpi,format='png')            
                
    plt.show()

    return

def save_mask(mask,shapefile):

    # Define a polygon feature geometry with one attribute
##    schema = {
##        'geometry': 'Polygon',
##        'properties': {'id': 'int'},
##    }
##
    schema = {
        'geometry':  mask.geom_type,
        'properties': {'id': 'int'},
    }

    # Write a new Shapefile
    with fiona.open(shapefile, 'w', 'ESRI Shapefile', schema) as c:
        ## If there are multiple geometries, put the "for" loop here
        c.write({
            'geometry': mapping(mask),
            'properties': {'id': 1},
        })

    return

def load_mask(shapefile):

    with fiona.collection(shapefile, "r") as input:
        schema = input.schema.copy()
        #geom = []
##        for f in input:
##            mask = shape(f['geometry'])
        f = input[0]
        mask = shape(f['geometry'])

    return mask    

def generate_vertices(nx,ny,mask,include_region,configfile):

    #(navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    config = ConfigParser.RawConfigParser()
    config.read(configfile)
    F_lon_min = config.get('GeoFrame', 'F_lon_min')
    F_lon_max = config.get('GeoFrame', 'F_lon_max')
    F_lat_min = config.get('GeoFrame', 'F_lat_min')
    F_lat_max = config.get('GeoFrame', 'F_lat_max')

    F_lon_min = float(F_lon_min)    
    F_lon_max = float(F_lon_max)    
    F_lat_min = float(F_lat_min)    
    F_lat_max = float(F_lat_max)    

    lons = linspace(F_lon_min,F_lon_max,nx+1)    
    lats = linspace(F_lat_min,F_lat_max,ny+1)
    ilons = len(lons)
    ilats = len(lats)

    nmask = len(mask.geoms)    

    verts = []
    name = []
    long_name = []
    lon_wmc = []
    lat_wmc = []
    count = 0
    for i in range(ilons-1):
        for j in range(ilats-1):
            vert = Polygon([(lons[i],lats[j]), \
                            (lons[i],lats[j+1]), \
                            (lons[i+1],lats[j+1]), \
                            (lons[i+1],lats[j]), \
                            (lons[i],lats[j])])
            for k in range(nmask):
                if ((k+1) in include_region):
                    intersect = vert.intersection(mask.geoms[k])
                    if (intersect.geom_type == 'Polygon'):
                        verts.append(intersect.convex_hull)
                        count = count + 1
                        name.append("A%03d" % count)
                        long_name.append("A%03d" % count)
                        lon_wmc.append(float("%9.4f" % intersect.centroid.x))
                        lat_wmc.append(float("%9.4f" % intersect.centroid.y))
                    elif (intersect.geom_type == 'MultiPolygon'):
                        for kk in range(len(intersect.geoms)):
                            if (intersect.geoms[kk].geom_type == 'Polygon'):
                                geom = intersect.geoms[kk]
                                verts.append(geom.convex_hull)
                                count = count + 1
                                name.append("A%03d" % count)
                                long_name.append("A%03d" % count)
                                lon_wmc.append(float("%9.4f" % geom.centroid.x))
                                lat_wmc.append(float("%9.4f" % geom.centroid.y))
                    else:
                        pass
                        #print "%d %d %s" % (i,j,intersect.geom_type)


    return (verts,name,long_name,lon_wmc,lat_wmc)

def generate_arbitrary_vertices(polygons,name,long_name,mask,include_region,configfile):

    #(navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    config = ConfigParser.ConfigParser()
    config.read(configfile)
    F_lon_min = config.get('GeoFrame', 'F_lon_min')
    F_lon_max = config.get('GeoFrame', 'F_lon_max')
    F_lat_min = config.get('GeoFrame', 'F_lat_min')
    F_lat_max = config.get('GeoFrame', 'F_lat_max')    

    nmask = len(mask.geoms)    

    verts = []
    lon_wmc = []
    lat_wmc = []
    count = 0
    for vert in polygons:
        for k in range(nmask):
            if ((k+1) in include_region):
                intersect = vert.intersection(mask.geoms[k])
                if (intersect.geom_type == 'Polygon'):
                    verts.append(intersect.convex_hull)
                    count = count + 1
                    lon_wmc.append(float("%9.4f" % intersect.centroid.x))
                    lat_wmc.append(float("%9.4f" % intersect.centroid.y))
                elif (intersect.geom_type == 'MultiPolygon'):
                    for kk in range(len(intersect.geoms)):
                        if (intersect.geoms[kk].geom_type == 'Polygon'):
                            geom = intersect.geoms[kk]
                            verts.append(geom.convex_hull)
                            count = count + 1
                            lon_wmc.append(float("%9.4f" % geom.centroid.x))
                            lat_wmc.append(float("%9.4f" % geom.centroid.y))
                else:
                    print intersect.geom_type

    return (verts,name,long_name,lon_wmc,lat_wmc)

def map_model_cells_to_vertices(meshfile,verts,name,mapfile):

    (navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    result = chararray((im, jm), max(len(nam) for nam in name))

    (pl,dict1) = initialize_PointLocator(verts,name)    

    fout = open(mapfile, 'w')
    
    for j in range(jm):
        for i in range(im):
            try:
                x1 = float('nan')
                y1 = float('nan')
                key1 =""
                x1 = pl.contains_point((navLonTran[i,j],navLatTran[i,j]))[0].centroid[0]
                y1 = pl.contains_point((navLonTran[i,j],navLatTran[i,j]))[0].centroid[1]
                key1 = "%15.5f %15.5f" % (round(x1,6),round(y1,6))
                result[i,j] = dict1[key1]
            except:
                result[i,j] = '*'
            fout.write('%s\t' % result[i,j])
        fout.write('\n')
    
    return    

def aggregate_vertices(verts,name,old_names,new_names,new_long_name):

    verts_aggr = []
    lon_wmc_aggr = []
    lat_wmc_aggr = []
    name_aggr = []
    long_name_aggr = []
    for i in range(len(old_names)):
        ind = []
        for j in range(len(name)):
            if (name[j] in old_names[i]):
                ind.append(j)
        verts_ind = []
        for index in ind:
            verts_ind.append(verts[index])
        verts_aggr.append(cascaded_union(verts_ind))
        name_aggr.append(new_names[i])
        long_name_aggr.append(new_long_name[i])
        lon_wmc_aggr.append(float("%9.4f" % verts_aggr[i].centroid.x))
        lat_wmc_aggr.append(float("%9.4f" % verts_aggr[i].centroid.y))

    return (verts_aggr,name_aggr,long_name_aggr,lon_wmc_aggr,lat_wmc_aggr)    

def save_vertices(verts,name,long_name,lon_wmc,lat_wmc,shapefile):

    vert_geom_type ='Polygon'    
    schema = {
        'geometry':  vert_geom_type,
        'properties': {'id': 'int',
                       'name': 'str',
                       'long_name': 'str',
                       'lon_wmc': 'float:9.4',
                       'lat_wmc': 'float:9.4'},
    }

    # Write a Shapefile
    with fiona.open(shapefile, 'w', 'ESRI Shapefile', schema) as c:
        ## If there are multiple geometries, put the "for" loop here
        for i in range(len(verts)):
            c.write({
            'geometry': mapping(verts[i]),
            'properties': {'id': i+1,
                           'name': name[i],
                           'long_name': long_name[i],
                           'lon_wmc': lon_wmc[i],
                           'lat_wmc': lat_wmc[i]},
            })

    return

def load_vertices(shapefile):

    with fiona.collection(shapefile, "r") as input:
        schema = input.schema.copy()
        verts = []
        name = []
        long_name = []
        lon_wmc = []
        lat_wmc = []
        for f in input:
            verts.append(shape(f['geometry']))
            name.append(f['properties']['name'])
            long_name.append(f['properties']['long_name'])
            lon_wmc.append(f['properties']['lon_wmc'])
            lat_wmc.append(f['properties']['lat_wmc'])

    return (verts,name,long_name,lon_wmc,lat_wmc)    


def show_vertices(verts,name,long_name,lon_wmc,lat_wmc,configfile):
    
    #(navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)

    config = ConfigParser.ConfigParser()
    config.read(configfile)
    F_lon_min = config.get('GeoFrame', 'F_lon_min')
    F_lon_max = config.get('GeoFrame', 'F_lon_max')
    F_lat_min = config.get('GeoFrame', 'F_lat_min')
    F_lat_max = config.get('GeoFrame', 'F_lat_max')    
    
    map = Basemap(llcrnrlon=F_lon_min,llcrnrlat=F_lat_min,urcrnrlon=F_lon_max,urcrnrlat=F_lat_max,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='white',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    
    nvert = len(verts)
    for i in range(nvert):
        if (verts[i].geom_type == 'Polygon'):
            x = []
            y = []
            polygon = list(verts[i].exterior.coords)
            lp = len(polygon)
            for j in range(lp):
                x.append(polygon[j][0])
                y.append(polygon[j][1])
            (xp,yp)=map(x,y)
            plt.fill(xp,yp,'blue',alpha=0.5)
            (xt,yt)=map(lon_wmc[i],lat_wmc[i])
            plt.text(xt,yt,name[i],fontsize=10,fontweight='normal',ha='center',va='center',color='black')
        else:
            print "%d %s" % (i, verts[i].geom_type)
            
    #plt.savefig(figfile,dpi=plt.gcf().dpi,format='pdf')            
    #plt.savefig(str(figfile).replace(".pdf",".png"),dpi=plt.gcf().dpi,format='png')            
                
    plt.show()        
                
    return

def show_vertices_mask(verts,name,lon_wmc,lat_wmc,mask,configfile):

    #(navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)
    
    config = ConfigParser.ConfigParser()
    config.read(configfile)
    F_lon_min = config.get('GeoFrame', 'F_lon_min')
    F_lon_max = config.get('GeoFrame', 'F_lon_max')
    F_lat_min = config.get('GeoFrame', 'F_lat_min')
    F_lat_max = config.get('GeoFrame', 'F_lat_max')    
    
    map = Basemap(llcrnrlon=F_lon_min,llcrnrlat=F_lat_min,urcrnrlon=F_lon_max,urcrnrlat=F_lat_max,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='white',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    
    nmask = len(mask.geoms)    
    for i in range(nmask):
        if (mask.geoms[i].geom_type == 'Polygon'):
            lon_wmc_mask = mask.geoms[i].centroid.x
            lat_wmc_mask = mask.geoms[i].centroid.y
            x = []
            y = []
            polygon = list(mask.geoms[i].exterior.coords)
            lp = len(polygon)
            for j in range(lp):
                x.append(polygon[j][0])
                y.append(polygon[j][1])
            (xp,yp)=map(x,y)
            plt.fill(xp,yp,'blue',alpha=0.5)
            (xt,yt)=map(lon_wmc_mask,lat_wmc_mask)
            plt.text(xt,yt,`i+1`,fontsize=24,fontweight='normal',ha='center',va='center',color='white')
            nisl = len(mask.geoms[i].interiors)
            if (nisl > 0):
                for k in range(nisl):
                    island = list(mask.geoms[i].interiors[k].coords)
                    x = []
                    y = []
                    li = len(island)
                    for j in range(li):
                        x.append(island[j][0])
                        y.append(island[j][1])
                    (xp,yp)=map(x,y)
                    plt.fill(xp,yp,'red',alpha=0.5)
        else:
            print mask.geoms[i].geom_type

    nvert = len(verts)
    for i in range(nvert):
        if (verts[i].geom_type == 'Polygon'):
            x = []
            y = []
            polygon = list(verts[i].exterior.coords)
            lp = len(polygon)
            for j in range(lp):
                x.append(polygon[j][0])
                y.append(polygon[j][1])
            (xp,yp)=map(x,y)
            plt.fill(xp,yp,'green',alpha=0.5)
            (xt,yt)=map(lon_wmc[i],lat_wmc[i])
            plt.text(xt,yt,name[i],fontsize=10,fontweight='normal',ha='center',va='center',color='black')
        else:
            print verts[i].geom_type
            
    #plt.savefig(figfile,dpi=plt.gcf().dpi,format='pdf')            
    #plt.savefig(str(figfile).replace(".pdf",".png"),dpi=plt.gcf().dpi,format='png')            
                
    plt.show()        
                
    return

def load_errfile(errfile):
        
    lines = open(errfile, "r").readlines()

    header = lines[0].split()
    data = lines[1:]
    lh = len(header)
    nerr = len(data)

    print nerr

    lon_err1 = zeros(nerr,dtype=float)
    lat_err1 = zeros(nerr,dtype=float)
    flag_err1 = zeros(nerr,dtype=int)
    lon_err2 = zeros(nerr,dtype=float)
    lat_err2 = zeros(nerr,dtype=float)
    flag_err2 = zeros(nerr,dtype=int)
    
    for i in range(nerr):     
        line1 = data[i].rpartition('\n')
        line2 = line1[0].split('\t')
        lon_err1[i] = float(line2[0])
        lat_err1[i] = float(line2[1])
        flag_err1[i] = int(line2[2])
        lon_err2[i] = float(line2[3])
        lat_err2[i] = float(line2[4])
        flag_err2[i] = int(line2[5])
            
    return (nerr,lon_err1,lat_err1,flag_err1,lon_err2,lat_err2,flag_err2)

def show_vert_mask_fsm_err(verts,name,lon_wmc,lat_wmc,mask,configfile,errfile):

    (nerr,lon_err1,lat_err1,flag_err1,lon_err2,lat_err2,flag_err2) = load_errfile(errfile)    

    #(navLonTran,navLatTran,tmaskTran,im,jm,F_lon_min,F_lon_max,F_lat_min,F_lat_max) = load_meshmask(meshfile)
    
    config = ConfigParser.ConfigParser()
    config.read(configfile)
    F_lon_min = config.get('GeoFrame', 'F_lon_min')
    F_lon_max = config.get('GeoFrame', 'F_lon_max')
    F_lat_min = config.get('GeoFrame', 'F_lat_min')
    F_lat_max = config.get('GeoFrame', 'F_lat_max')    
    
    map = Basemap(llcrnrlon=F_lon_min,llcrnrlat=F_lat_min,urcrnrlon=F_lon_max,urcrnrlat=F_lat_max,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='white',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every five degrees.
    map.drawmeridians(arange(0,360,5),labels=[0,0,0,1],fontsize=10)
    map.drawparallels(arange(-90,90,5),labels=[1,0,0,0],fontsize=10)
    
    nmask = len(mask.geoms)    
    for i in range(nmask):
        if (mask.geoms[i].geom_type == 'Polygon'):
            lon_wmc_mask = mask.geoms[i].centroid.x
            lat_wmc_mask = mask.geoms[i].centroid.y
            x = []
            y = []
            polygon = list(mask.geoms[i].exterior.coords)
            lp = len(polygon)
            for j in range(lp):
                x.append(polygon[j][0])
                y.append(polygon[j][1])
            (xp,yp)=map(x,y)
            plt.fill(xp,yp,'blue',alpha=0.5)
            (xt,yt)=map(lon_wmc_mask,lat_wmc_mask)
            plt.text(xt,yt,`i+1`,fontsize=24,fontweight='normal',ha='center',va='center',color='white')
            nisl = len(mask.geoms[i].interiors)
            if (nisl > 0):
                for k in range(nisl):
                    island = list(mask.geoms[i].interiors[k].coords)
                    x = []
                    y = []
                    li = len(island)
                    for j in range(li):
                        x.append(island[j][0])
                        y.append(island[j][1])
                    (xp,yp)=map(x,y)
                    plt.fill(xp,yp,'red',alpha=0.5)
        else:
            print mask.geoms[i].geom_type

    nvert = len(verts)
    for i in range(nvert):
        if (verts[i].geom_type == 'Polygon'):
            x = []
            y = []
            polygon = list(verts[i].exterior.coords)
            lp = len(polygon)
            for j in range(lp):
                x.append(polygon[j][0])
                y.append(polygon[j][1])
            (xp,yp)=map(x,y)
            plt.fill(xp,yp,'green',alpha=0.5)
            (xt,yt)=map(lon_wmc[i],lat_wmc[i])
            plt.text(xt,yt,name[i],fontsize=10,fontweight='normal',ha='center',va='center',color='black')
        else:
            print verts[i].geom_type

    (xe1,ye1)=map(lon_err1,lat_err1)
    map.scatter(xe1,ye1,10,color='red',zorder=1000)

    (xe2,ye2)=map(lon_err2,lat_err2)
    map.scatter(xe2,ye2,10,color='yellow',zorder=1000)

        
    plt.show()        
                
    return
        
def show_traj_info(fname):

    f = Dataset(str(fname).replace("/","\\"), 'r', format='NETCDF3_CLASSIC')    

    traj_lon = f.variables['traj_lon'][:,:]

    trajLonShape = traj_lon.shape
    ltraj = trajLonShape[0]
    ntraj = trajLonShape[1]

    f.close()

    del traj_lon

    return (ltraj,ntraj)    

def show_traj(fname,nplot,F_lon_min,F_lon_max,F_lat_min,F_lat_max):

    f = Dataset(str(fname).replace("/","\\"), 'r', format='NETCDF3_CLASSIC')

    traj_lon = f.variables['traj_lon'][:,:]
    
    trajLonShape = traj_lon.shape
    ltraj = trajLonShape[0]
    ntraj = trajLonShape[1]
    traj_lon[nonzero(traj_lon>1.e19)]='NaN'

    
    traj_lat = f.variables['traj_lat'][:,:]
    trajLatShape = traj_lat.shape 
    traj_lat[nonzero(traj_lat>1.e19)]='NaN'

    f.close()

    map = Basemap(llcrnrlon=F_lon_min,llcrnrlat=F_lat_min,urcrnrlon=F_lon_max,urcrnrlat=F_lat_max,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    #map.fillcontinents(color='coral',lake_color='aqua')
    map.fillcontinents(color='white',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every five degrees.
    map.drawmeridians(arange(0,360,5),labels=[0,0,0,1],fontsize=10)
    map.drawparallels(arange(-90,90,5),labels=[1,0,0,0],fontsize=10)

    (xp,yp) = map(traj_lon[0,:],traj_lat[0,:])
    map.scatter(xp,yp,10,marker='o',color='k',zorder=100)

    color = ['m','c','r','g','b','y']
    lc = size(color)
    count = 0
    if (nplot > ntraj):
        nstep = 1
    else:
        nstep = ntraj/nplot   
    for itraj in range(0, ntraj, nstep):
        ic = mod(count,lc)
        (xp,yp) = map(traj_lon[:,itraj],traj_lat[:,itraj])
        map.plot(xp, yp, color=color[ic], linewidth=1.5,zorder=1)
        count = count + 1

    del traj_lon, traj_lat    
    
    plt.show()    

def generate_examples(verts,name,trajfile,base_year,base_month,base_day,yearmin,monthmin,yearmax,monthmax,depthmin,depthmax, \
                      outfile,delta_t,samp_int):

    (pl,dict1) = initialize_PointLocator(verts,name)

    #base_date = datetime(1999,1,2)
    base_date = datetime(base_year,base_month,base_day)

    yyyymm_min = yearmin * 100 + monthmin
    yyyymm_max = yearmax * 100 + monthmax        

    f = Dataset(str(trajfile).replace("/","\\"), 'r', format='NETCDF3_CLASSIC')

    init_t = f.variables['init_t'][:]    

    traj_lon = f.variables['traj_lon'][:,:]

    trajLonShape = traj_lon.shape
    ltraj = trajLonShape[0]
    ntraj = trajLonShape[1]

    print ltraj
    print ntraj

    f.close()

    fileout=open(outfile,'w')
    fileout.write("year\tmonth\tarea_this_epoch\tarea_next_epoch\n")
    fileout.write("d\td\td\td\n")
    fileout.write("\t\t\tc\n")    

    errfile = outfile.replace('.txt','.err.txt')
    fileerr=open(errfile,'w')
    fileerr.write("lon_err1\tlat_err1\tflag_err1\tlon_err2\tlat_err2\tflag_err2\n")

    num_out = 0
    num_err = 0
    for itraj in range(ntraj):
        f = Dataset(str(trajfile).replace("/","\\"), 'r', format='NETCDF3_CLASSIC')
        
        traj_lon = f.variables['traj_lon'][:,itraj]
        traj_lat = f.variables['traj_lat'][:,itraj]
        traj_depth = f.variables['traj_depth'][:,itraj]

        f.close()

        traj_year = zeros((ltraj),dtype=int)    
        traj_month = zeros((ltraj),dtype=int)    
        for traj_day in range(0,ltraj-delta_t,samp_int):
            if ((traj_lon[traj_day]<1.e19) and
                (traj_lon[traj_day+delta_t]<1.e19)):
                time_traj = base_date + ((int(init_t[itraj]))+traj_day)*timedelta(days=1)
                traj_year[traj_day] = time_traj.year
                traj_month[traj_day] = time_traj.month
                year_this_epoch = time_traj.year
                month_this_epoch = time_traj.month
                depth_this_epoch = traj_depth[traj_day]
                yyyymm = year_this_epoch * 100 + month_this_epoch
                if ((yyyymm >= yyyymm_min) & (yyyymm <= yyyymm_max) &
                    (depth_this_epoch < depthmin) & (depth_this_epoch > depthmax)):
                    flag_err1 = 0
                    flag_err = 0
                    x1 = float('nan')
                    y1 = float('nan')
                    key1 =""
                    try:
                        x1 = pl.contains_point((traj_lon[traj_day],traj_lat[traj_day]))[0].centroid[0]
                        y1 = pl.contains_point((traj_lon[traj_day],traj_lat[traj_day]))[0].centroid[1]
                        flag_err1 = 1
                        key1 = "%15.5f %15.5f" % (round(x1,6),round(y1,6))
                        flag_err1 = 2
                        area_this_epoch = dict1[key1]
                        flag_err1 = 3                        
                    except:
                        flag_err = 1
                        #fileerr.write("1\t%15.9f\t%15.9f\t%d\n" %  (traj_lon[traj_day,itraj],traj_lat[traj_day,itraj],flag_err1))
                    flag_err2 = 0
                    x2 = float('nan')
                    y2 = float('nan')
                    key2 =""
                    try:
                        x2 = pl.contains_point((traj_lon[traj_day+delta_t],traj_lat[traj_day+delta_t]))[0].centroid[0]
                        y2 = pl.contains_point((traj_lon[traj_day+delta_t],traj_lat[traj_day+delta_t]))[0].centroid[1]
                        flag_err2 = 1
                        key2 = "%15.5f %15.5f" % (round(x2,6),round(y2,6))
                        flag_err2 = 2
                        area_next_epoch = dict1[key2]
                        flag_err2 = 3                        
                    except:
                        flag_err = 1
                        #fileerr.write("2\t%15.9f\t%15.9f\t%d\n" %  (traj_lon[traj_day+delta_t,itraj],traj_lat[traj_day+delta_t,itraj],flag_err2))
                    if (flag_err == 0):
                        num_out = num_out + 1
                        fileout.write("%d\t%d\t%s\t%s\n" %
                                      (year_this_epoch,month_this_epoch,area_this_epoch,area_next_epoch))
                    else:
                        num_err = num_err + 1
                        fileerr.write("%15.9f\t%15.9f\t%d\t%15.9f\t%15.9f\t%d\n" %
                                      (traj_lon[traj_day],traj_lat[traj_day],flag_err1,
                                       traj_lon[traj_day+delta_t],traj_lat[traj_day+delta_t],flag_err2))
                    
##                    nex_out = nex_out + 1
##                    #print nex_out
##                    try:
##                        key1 = "%15.9f %15.9f" % pl.contains_point((traj_lon[traj_day],traj_lat[traj_day]))[0].centroid
##                        area_this_epoch = dict1[key1]
##                        x = pl.contains_point((traj_lon[traj_day+samp_int],traj_lat[traj_day+samp_int]))[0].centroid[0]
##                        y = pl.contains_point((traj_lon[traj_day+samp_int],traj_lat[traj_day+samp_int]))[0].centroid[1]
##                        key2 = "%15.9f %15.9f" % (x,y)
##                        area_next_epoch = dict1[key2]
##                        fileout.write("%d\t%d\t%s\t%s\n" %
##                                      (year_this_epoch,month_this_epoch,area_this_epoch,area_next_epoch))
##                    except:
##                        print 'error: ' + `key1` + ' ' + `key2`

        del traj_lon
        del traj_lat
        del traj_depth
            
    fileout.close()
    fileerr.close()
   
    return (num_out,num_err)

def generate_examples_fast(verts,name,trajfile,base_year,base_month,base_day,yearmin,monthmin,yearmax,monthmax,depthmin,depthmax, \
                      outfile,delta_t,samp_int):

    (pl,dict1) = initialize_PointLocator(verts,name)

    #base_date = datetime(1999,1,2)
    base_date = datetime(base_year,base_month,base_day)

    yyyymm_min = yearmin * 100 + monthmin
    yyyymm_max = yearmax * 100 + monthmax        

    f = Dataset(str(trajfile).replace("/","\\"), 'r', format='NETCDF3_CLASSIC')

    init_t = f.variables['init_t'][:]    

    traj_lon = f.variables['traj_lon'][:,:]

    trajLonShape = traj_lon.shape
    traj_lat = f.variables['traj_lat'][:,:]
    traj_depth = f.variables['traj_depth'][:,:]
    ltraj = trajLonShape[0]
    ntraj = trajLonShape[1]

    print ltraj
    print ntraj

    f.close()

    fileout=open(outfile,'w')
    fileout.write("year\tmonth\tarea_this_epoch\tarea_next_epoch\n")
    fileout.write("d\td\td\td\n")
    fileout.write("\t\t\tc\n")

    errfile = outfile.replace('.txt','.err.txt')
    fileerr=open(errfile,'w')
    fileerr.write("lon_err1\tlat_err1\tflag_err1\tlon_err2\tlat_err2\tflag_err2\n")

    num_out = 0
    num_err = 0
    for itraj in range(ntraj):
        traj_year = zeros((ltraj),dtype=int)    
        traj_month = zeros((ltraj),dtype=int)    
        for traj_day in range(0,ltraj-delta_t,samp_int):
            if ((traj_lon[traj_day,itraj]<1.e19) and
                (traj_lon[traj_day+delta_t,itraj]<1.e19)):
                time_traj = base_date + ((int(init_t[itraj]))+traj_day)*timedelta(days=1)
                traj_year[traj_day] = time_traj.year
                traj_month[traj_day] = time_traj.month
                year_this_epoch = time_traj.year
                month_this_epoch = time_traj.month
                depth_this_epoch = traj_depth[traj_day,itraj]
                yyyymm = year_this_epoch * 100 + month_this_epoch
                if ((yyyymm >= yyyymm_min) & (yyyymm <= yyyymm_max) &
                    (depth_this_epoch < depthmin) & (depth_this_epoch > depthmax)):
                    flag_err1 = 0
                    flag_err = 0
                    x1 = float('nan')
                    y1 = float('nan')
                    key1 =""
                    try:
                        x1 = pl.contains_point((traj_lon[traj_day,itraj],traj_lat[traj_day,itraj]))[0].centroid[0]
                        y1 = pl.contains_point((traj_lon[traj_day,itraj],traj_lat[traj_day,itraj]))[0].centroid[1]
                        flag_err1 = 1
                        key1 = "%15.5f %15.5f" % (round(x1,6),round(y1,6))
                        flag_err1 = 2
                        area_this_epoch = dict1[key1]
                        flag_err1 = 3                        
                    except:
                        flag_err = 1
                        #fileerr.write("1\t%15.9f\t%15.9f\t%d\n" %  (traj_lon[traj_day,itraj],traj_lat[traj_day,itraj],flag_err1))
                    flag_err2 = 0
                    x2 = float('nan')
                    y2 = float('nan')
                    key2 =""
                    try:
                        x2 = pl.contains_point((traj_lon[traj_day+delta_t,itraj],traj_lat[traj_day+delta_t,itraj]))[0].centroid[0]
                        y2 = pl.contains_point((traj_lon[traj_day+delta_t,itraj],traj_lat[traj_day+delta_t,itraj]))[0].centroid[1]
                        flag_err2 = 1
                        key2 = "%15.5f %15.5f" % (round(x2,6),round(y2,6))
                        flag_err2 = 2
                        area_next_epoch = dict1[key2]
                        flag_err2 = 3                        
                    except:
                        flag_err = 1
                        #fileerr.write("2\t%15.9f\t%15.9f\t%d\n" %  (traj_lon[traj_day+delta_t,itraj],traj_lat[traj_day+delta_t,itraj],flag_err2))
                    if (flag_err == 0):
                        num_out = num_out + 1
                        fileout.write("%d\t%d\t%s\t%s\n" %
                                      (year_this_epoch,month_this_epoch,area_this_epoch,area_next_epoch))
                    else:
                        num_err = num_err + 1
                        fileerr.write("%15.9f\t%15.9f\t%d\t%15.9f\t%15.9f\t%d\n" %
                                      (traj_lon[traj_day,itraj],traj_lat[traj_day,itraj],flag_err1,
                                       traj_lon[traj_day+delta_t,itraj],traj_lat[traj_day+delta_t,itraj],flag_err2))
   
    fileout.close()
    fileerr.close()
   
    return (num_out,num_err)

def initialize_PointLocator(verts,name):

    vertpl = []
    dict1 ={}
    for i in range(len(verts)):
        vertpl.append(cg.asShape(verts[i]))
        key = "%15.5f %15.5f" % (round(verts[i].centroid.x,6),round(verts[i].centroid.y,6))
        #dict1[key] = "A%03d" % (i+1)
        dict1[key] = name[i]

    pl = cg.PolygonLocator(vertpl)    

    return (pl,dict1)

def generate_assoc_rules(fvert,fname,minSupport,minConf,maxIS,fout):

    fname1 = "%s" % fname
    
    data = orange.ExampleTable(fname1)

    rules = orange.AssociationRulesInducer(data, support=minSupport, confidence=minConf, maxItemSets=maxIS, classification_rules=True)

    print "%i rules with support higher than or equal to %9.6f found." % (len(rules), minSupport)

    (verts,area_name,long_name,area_lon_wmc,area_lat_wmc) = load_vertices(fvert)

    narea = len(verts)    

    #(narea,area_name,area_lon_wmc,area_lat_wmc,area_poly)=read_areas(fvert)

    fname2 = "%s" % fout

    fid=open(fname2,'w')
    fid.write("supp\tconf\tcove\tyear\tmonth\tarea_this_epoch\tlonm_this_epoch\tlatm_this_epoch\tarea_next_epoch\tlonm_next_epoch\tlatm_next_epoch\n")

    nr = 0

    patternYMA = re.compile("year=\d+.month=\d+.area_this_epoch=[a-z]+", re.IGNORECASE)
    patternYA = re.compile("year=\d+.area_this_epoch=[a-z]+", re.IGNORECASE)
    patternMA = re.compile("month=\d+.area_this_epoch=[a-z]+", re.IGNORECASE)
    patternA = re.compile("area_this_epoch=[a-z]+", re.IGNORECASE)

    ls = len(rules)

    for i in range(ls):
        rule = "%s" % rules[i]
        mtch = patternYMA.match(rule)
        found = False
        if mtch!=None:
            aa = re.search("year=\d+",rule)
            a1 = aa.start() + 5
            a2 = aa.end()
            year = rule[a1:a2]
            aa = re.search("month=\d+",rule)
            a1 = aa.start() + 6
            a2 = aa.end()
            month = rule[a1:a2]
            found = True
        mtch = patternYA.match(rule)
        if mtch!=None:
            aa = re.search("year=\d+",rule)
            a1 = aa.start() + 5
            a2 = aa.end()
            year = rule[a1:a2]
            month = "99"
            found = True
        mtch = patternMA.match(rule)
        if mtch!=None:
            year = '9999'
            aa = re.search("month=\d+",rule)
            a1 = aa.start() + 6
            a2 = aa.end()
            month = rule[a1:a2]
            found = True
        mtch = patternA.match(rule)
        if mtch!=None:
            year = "9999"
            month = "99"
            found = True
        if found:
            aa = re.search("area_this_epoch=\w+",rule)
            a1 = aa.start() + 16
            a2 = aa.end()
            area_this_epoch = rule[a1:a2]
            for j in range(narea):
                if area_name[j]==area_this_epoch:
                    lonm_this_epoch = area_lon_wmc[j]
                    latm_this_epoch = area_lat_wmc[j]
                    break
            aa = re.search("area_next_epoch=\w+",rule)
            a1 = aa.start() + 16
            a2 = aa.end()
            area_next_epoch = rule[a1:a2]
            for j in range(narea):
                if area_name[j]==area_next_epoch:
                    lonm_next_epoch = area_lon_wmc[j]
                    latm_next_epoch = area_lat_wmc[j]
                    break

            nr = nr + 1
            fid.write("%f\t" % rules[i].support)
            fid.write("%f\t" % rules[i].confidence)
            fid.write("%f\t" % rules[i].coverage)
            fid.write("%s\t" % year)
            fid.write("%s\t" % month)
            fid.write("%s\t" % area_this_epoch)
            fid.write("%f\t" % lonm_this_epoch)
            fid.write("%f\t" % latm_this_epoch)
            fid.write("%s\t" % area_next_epoch)
            fid.write("%f\t" % lonm_next_epoch)
            fid.write("%f\n" % latm_next_epoch)
        else:
            pass
            #print "%d\t%s" % (i,rules[i])

    fid.close()        
                    
    print ("Number of matched rules: %d" % nr)

    return (len(rules),nr)


def show_digraph(fvert,fname,minsup,minconf,scale,labels,fontsize,year,month):

    (verts,area_name,area_lon_wmc,area_lat_wmc) = load_vertices(fvert)

    narea = len(verts)
    
    #(narea,area_name,area_lon_wmc,area_lat_wmc,area_poly)=read_areas(fvert)

    (supp,conf,cove,year_this_epoch,month_this_epoch,area_this_epoch,lonm_this_epoch,latm_this_epoch,area_next_epoch,lonm_next_epoch,latm_next_epoch,nrules)=read_assoc_rules_areas(fname)

    ## MFS frame ##
    F_lon_min = -5.5
    F_lon_max = 36.0
    F_lat_min = 30.24
    F_lat_max = 45.9587212

    map = Basemap(llcrnrlon=F_lon_min,llcrnrlat=F_lat_min,urcrnrlon=F_lon_max,urcrnrlat=F_lat_max,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc')
    
    for i in range(narea):
        x = []
        y = []
        polygon = list(verts[i].exterior.coords)
        lp = len(polygon)
        for j in range(lp):
            x.append(polygon[j][0])
            y.append(polygon[j][1])
        (xp,yp)=map(x,y)
        plt.fill(xp,yp,facecolor='none', edgecolor='k')

    cdict1 = {'red':   ((0.0, 1.0, 1.0),
                        (1.0, 1.0, 1.0)),
              'green': ((0.0, 1.0, 1.0),
                        (1.0, 0.0, 0.0)),
              'blue':  ((0.0, 1.0, 1.0),
                        (1.0, 0.0, 0.0))
              }

    c_self = LinearSegmentedColormap('c_self', cdict1)
    plt.register_cmap(cmap=c_self)

    patches = []
    colors = zeros(1000)
    k = 0

    for i in range(nrules):
        if (year_this_epoch[i] == year and month_this_epoch[i] == month and area_this_epoch[i] == area_next_epoch[i]):
            for j in range(narea):
                if area_name[j] == area_this_epoch[i]:
                    if (verts[j].geom_type == 'Polygon'):
                        x = []
                        y = []
                        polygon = list(verts[j].exterior.coords)
                        lp = len(polygon)
                        for kk in range(lp):
                            x.append(polygon[kk][0])
                            y.append(polygon[kk][1])
                        (xp, yp) = map(x, y)
                        xpyp = zeros((lp, 2), dtype=float)                   
                        for kk in range(lp):
                            xpyp[kk, 0] = xp[kk]
                            xpyp[kk, 1] = yp[kk]
                        poly = ptch.Polygon(xpyp, closed=True)
                        patches.append(poly)
                        colors[k] = 100 * conf[i]
                        k = k + 1
                        if (labels == 2):
                            (xt, yt) = map(lonm_this_epoch[i], latm_this_epoch[i])
                            plt.text(xt, yt, str(round(conf[i] * 100, 1)), fontsize=fontsize, fontweight='normal', ha='center', va='center', color='b')
                        break
    p = PatchCollection(patches, cmap=c_self, alpha=1) 
    p.set_array(colors)
    p.set_clim([0,100])    
    ax = plt.gca()
    ax.add_collection(p)
    cbar = map.colorbar(p, location='bottom')
    cbar.set_label('Probability (%)')

    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='white', lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every five degrees.
    map.drawmeridians(arange(0,360,5),labels=[0,0,1,0],fontsize=10)
    map.drawparallels(arange(-90,90,5),labels=[1,0,0,0],fontsize=10)

    for i in range(nrules):
        if (year_this_epoch[i] == year and month_this_epoch[i] == month and area_this_epoch[i] != area_next_epoch[i] and
            supp[i] >= minsup and conf[i] >= minconf):
            (xbm, ybm) = map(lonm_this_epoch[i], latm_this_epoch[i])
            (xm, ym) = map(lonm_next_epoch[i], latm_next_epoch[i])
            arr = ptch.FancyArrow(xbm, ybm, xm - xbm, ym - ybm, width=scale * conf[i], length_includes_head=True, color='g')
            ax = plt.gca()
            ax.add_patch(arr)
            if (labels == 2):
                att = plt.text(xbm + (xm - xbm) / 3, ybm + (ym - ybm) / 3, str(round(conf[i] * 100, 1)), fontsize=fontsize,
                               fontweight='normal', ha='center', va='center', color='b')
                ax.add_artist(att)

    # legend of edges    
    
    LConfStart = 10                    # legend confidence minimum
    LConfStep = 10                     # legend confidence step
    LXPosStart = 29.5                  # legend edge start x position
    LXPosEnd = 34.5                    # legend edge end x position
    LYPosStart = 46.5                  # legend edge start y position
    LYPosEnd = LYPosStart              # legend edge end y position
    LYStepFactor = -1. / 5.             # legend edges step factor
    LTextDis = 1.0                     # legend text displacement
    
    confm = int(array(conf).max(0) * 100)
    confm = 50
    for confl in range(LConfStart, confm, LConfStep):
        xpos1 = LXPosStart
        ypos1 = LYPosStart + confl * LYStepFactor
        xpos2 = LXPosEnd
        ypos2 = ypos1
        (xbm, ybm) = map(xpos1, ypos1)
        (xm, ym) = map(xpos2, ypos2)
        arr = ptch.FancyArrow(xbm, ybm, xm - xbm, ym - ybm, width=scale * confl / 100, length_includes_head=True, color='g')
        ax = plt.gca()
        ax.add_patch(arr)
        
        xpost = xpos2 + LTextDis
        ypost = ypos1
        conf1 = confl - (LConfStep / 2)
        conf2 = confl + (LConfStep / 2)
        if conf2 > 100:
            conf = 100
        (xp, yp) = map(xpost, ypost)
        att = plt.text(xp, yp, '%d %s' % (confl, '%'), fontsize=fontsize,
                       fontweight='normal', ha='center', va='center', color='k')
        ax.add_artist(att)
    
    ypost = ypost + 10 * LYStepFactor       
    (xp, yp) = map(xpost, ypost)
    att = plt.text(xp, yp, 'Confidence', fontsize=fontsize, fontweight='normal', ha='center', va='center', color='k')            
    
    plt.show()        

def read_assoc_rules(arfile):

    lines = open(arfile, "r").readlines()

    header = lines[0].split()
    data = lines[1:]
    lh = len(header)
    ld = len(data)

    for j in range(lh):
        str = "%s=[]" % (header[j])
        exec str

    for i in range(ld):
        list1 = data[i].split()
        list = []
        for j in range(lh):
            try:
                list.append(float(list1[j]))
                str = "%s.append(%f)" % (header[j],float(list1[j]))
            except:
                list.append(list1[j])
                str = "%s.append(\'%s\')" % (header[j],list1[j])
            else:
                pass
            exec str

    nrules = ld
                    
    return (supp,conf,cove,year,month,area_this_epoch,lonm_this_epoch,latm_this_epoch,area_next_epoch,lonm_next_epoch,latm_next_epoch,nrules)

def generate_graphs(model, vertfile, yearmin, yearmax, monthmin, monthmax, timeinter, depthmin, \
                     depthmax, graphpath,  rulefile, statfile):
    
    str1 = re.search(model+'_rules_', rulefile)
    str2 = re.search('\d\d\d\ddays', rulefile[str1.end():])
    granulation = rulefile[str1.end():str1.end()+str2.start()-1]

    (verts,name,long_name,lon_wmc,lat_wmc) = load_vertices(vertfile)    
    
    (supp,conf,cove,year_this_epoch,month_this_epoch,area_this_epoch,lonm_this_epoch,latm_this_epoch, \
     area_next_epoch,lonm_next_epoch,latm_next_epoch,nrules)=read_assoc_rules(rulefile)

    years = unique(year_this_epoch)
    months = unique(month_this_epoch)
    g = []
    dict1 = dict()
    ngraph = 0
    for year in years:
        for month in months:
            g1 = igraph.Graph(len(verts),directed=True)
            g1["name"] = "%04d%02d" % (year,month)
            g1.vs["name"] = name
            g1.vs["long_name"] = long_name
            g1.vs["lon_wmc"] = lon_wmc
            g1.vs["lat_wmc"] = lat_wmc
            dict1[g1["name"]] = ngraph
            ngraph = ngraph + 1
            g.append(g1)

    for irule in range(len(supp)):
        key1 = "%04d%02d" % (year_this_epoch[irule],month_this_epoch[irule])
        ind1 = dict1[key1]
        vs = g[ind1].vs["name"]
        fr = vs.index(area_this_epoch[irule])
        to = vs.index(area_next_epoch[irule])
        g[ind1].add_edges((fr,to))
        eid = g[ind1].get_eid(fr,to)
        g[ind1].es[eid]["conf"] = conf[irule]
        g[ind1].es[eid]["supp"] = supp[irule]

    fileout=open(statfile,'w')
    fileout.write("year\tmonth\tvertices\tedges\tminsupp\tmaxsupp\tminconf\tmaxconf\n")    

    for i in range(ngraph):
        graphfile = graphpath + '/' + model + '_' + granulation + '_' + 'rules_'
        graphfile = graphfile + "%04d" % timeinter + "days_"
        year1 = g[i]["name"][0:4]
        month1 = g[i]["name"][4:6]
        if (year1 == '9999'):
            graphfile = graphfile + "%04d_%04d_" % (yearmin,yearmax)
            if (month1 == '99'):
                pass
            else:
                graphfile = graphfile + "%s_" % (month1)
        else:
            graphfile = graphfile + "%s_" % (year1)
            if (month1 == '99'):
                pass
            else:
                graphfile = graphfile + "%s_" % (month1)
        graphfile = graphfile + 'depth' + '_' + "%04d_%04d" % (abs(depthmin),abs(depthmax)) + '.gml'
        g[i].write(graphfile,format='graphml')
        if (g[i].ecount() > 0):
            minsupp = min(g[i].es.get_attribute_values('supp'))
            maxsupp = max(g[i].es.get_attribute_values('supp'))
            minconf = min(g[i].es.get_attribute_values('conf'))
            maxconf = max(g[i].es.get_attribute_values('conf'))
            fileout.write("%s\t%s\t%d\t%d\t%f\t%f\t%f\t%f\n" %
                          (year1,month1,g[i].vcount(),g[i].ecount(),minsupp, maxsupp, minconf, maxconf))

    fileout.close()        

    return ngraph 

def calc_elapsed_time(start_time):

    elapsed_time = time.time() - start_time
    m, s = divmod(elapsed_time, 60)
    h, m = divmod(m, 60)

    return (h, m, s)

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
