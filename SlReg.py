import os
import numpy as np
import sys
import cf
from netCDF4 import Dataset
from amrfile import io as amrio
import mapping_class as mp
import math
#from arg_to_file_exist import arg_to_file_exist
from CfReg import CfReg


class SlReg(CfReg):

  def __init__(self,argin,argout,argin_hdf):
    #err = 0
    #self.hdfplotfile, err  = arg_to_file_exist(argin_hdf, mandatory=True, err=err)
    #if err > 0:
    #  print 'Error opening ', argin_hdf
    #  sys.exit(2)
    if not os.path.exists(argin_hdf):
      print 'Error opening ', argin_hdf
      sys.exit(2)
    self.hdfplotfile = argin_hdf 
    super(SlReg,self).__init__(argin,argout)
  
  
  # diambil dari TEMP/ExtrParam.py
  def Distance(self,lat1,lon1,lat2,lon2):
    rad = 6371229.0
    lat1=lat1*math.pi/180 ;  lon1=lon1*math.pi/180  ;  lat2=lat2*math.pi/180 ;   lon2=lon2*math.pi/180
    dlat = lat2-lat1
    midlat =  (lat2+lat1)/2
  
    if lon1*lon2 >= -4:  # apart from Ross & Fimbul, must be pos. Fimbul will be neg small   ????
       dlon = lon2-lon1
       #dlon = locmax[0,1]-locmin[0,1]
    else: # Ross shelf
       dlon = (math.pi - np.sign(lon1)*lon1) + (math.pi - np.sign(lon2)*lon2)
  
    dist = rad*(dlat**2 + (math.cos(midlat)*dlon)**2)**0.5
    return dist
  
  
  def regrid(self):
    bisicles_to_nemo_mapping_file = '/projects/jomp/asiaha/TEMP/DataSementara/u-bh845/bisicles-AIS_lev0_to_eORCA.map2d'
    
    level=0 # level of grid refinement
    order=1 # interpolation order, 0 for piecewise constant, 1 for linear
    l_verb = True
    
    if (l_verb): print 'reading hdf file'
    amrID = amrio.load(self.hdfplotfile)
    lo,hi = amrio.queryDomainCorners(amrID, level)
    _,_,isf = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
    _,_,bathy = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)
    bathy=np.transpose(bathy)
    isf=np.transpose(isf)
    
    bn_map=mp.load(bisicles_to_nemo_mapping_file)
    
    latb = Dataset('/projects/ukesm/rsmith/ancils/cf_gridfile_BISICLES_lev0-AIS.nc').variables['latitude'][:]
    lonb = Dataset('/projects/ukesm/rsmith/ancils/cf_gridfile_BISICLES_lev0-AIS.nc').variables['longitude'][:]
    
    
    ny_n=100
    #ny_n=np.shape(bn_map.x)[0]
    nx_n=np.shape(bn_map.x)[1]
    
    bikemelt = np.zeros_like(isf)
    bikeheat = np.zeros_like(isf)
    
    
    if (l_verb): print "remap geometry to NEMO"
    for j in range(ny_n):
      if (l_verb): print '\r',j,'  ',
      sys.stdout.flush()
      print 'j, sum(bikemelt) = ', j, np.sum(bikemelt)
      for i in range(nx_n):
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           slope = np.zeros(ncontrib)
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               if (np.abs(isf[jj,ji] - bathy[jj,ji])  > 1E-3) \
               & (np.abs(isf[jj,ji]) > 1E-3):
                   if (i == 187) & (j == 53):
                      print 'jcont, ji, jj = ', jcont, ji, jj
                   slplon = [];slplat=[]
                   if (np.abs(isf[jj,ji+1] - bathy[jj,ji+1])  > 1E-3) & (np.abs(isf[jj,ji+1]) > 1E-3) :
                      slplonpl = np.abs(isf[jj,ji] - isf[jj,ji+1])/self.Distance(latb[jj,ji],lonb[jj,ji],latb[jj,ji+1],lonb[jj,ji+1])
                      slplon.append(slplonpl)
                   #else:
                   #   slplonpl = 0.0
    
                   if (np.abs(isf[jj,ji-1] - bathy[jj,ji-1])  > 1E-3) & (np.abs(isf[jj,ji-1]) > 1E-3) :
                      slplonne = np.abs(isf[jj,ji] - isf[jj,ji-1])/self.Distance(latb[jj,ji],lonb[jj,ji],latb[jj,ji-1],lonb[jj,ji-1])
                      slplon.append(slplonne)
                   #else:
                   #   slplonne = 0.0
    
                   if (np.abs(isf[jj+1,ji] - bathy[jj+1,ji])  > 1E-3) & (np.abs(isf[jj+1,ji]) > 1E-3) :
                      slplatpl = np.abs(isf[jj,ji] - isf[jj+1,ji])/self.Distance(latb[jj,ji],lonb[jj,ji],latb[jj+1,ji],lonb[jj+1,ji])
                      slplat.append(slplatpl)
                   #else:
                   #  slplatpl = 0.0
    
                   if (np.abs(isf[jj-1,ji] - bathy[jj-1,ji])  > 1E-3) & (np.abs(isf[jj-1,ji]) > 1E-3) :
                     slplatne = np.abs(isf[jj,ji] - isf[jj-1,ji])/self.Distance(latb[jj,ji],lonb[jj,ji],latb[jj-1,ji],lonb[jj-1,ji])
                     slplat.append(slplatne)
                   #else:
                   #  slplatne = 0.0
    
                   if len(slplon) > 0:
                      slplon = sum(slplon)/len(slplon)
                   else:
                      slplon = 0.0
    
                   if len(slplat) > 0:
                      slplat = sum(slplat)/len(slplat)
                   else:
                      slplat = 0.0
                   #slope[jcont] = np.max(np.array([slplonpl,slplonne,slplatpl,slplatne]))
                   slope[jcont] = math.sqrt(slplon**2 + slplat**2)
                   if math.isnan(slope[jcont]):
                      print 'i, j, jcont, ji, jj Nan = ', i, j, jcont, ji, jj
               else:
                 slope[jcont] = 0.0
    
    
           sigslope = np.sum(slope)
           totcav = np.sum(slope>0.0)
           scaling = totcav/sigslope*self.melt_water[j,i]
           scalingheat = totcav/sigslope*self.melt_heat[j,i]
           if math.isnan(scaling):
              print 'i, j, ncontrib = ', i, j, ncontrib, totcav, sigslope, self.melt_water[j,i]
    
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               bikemelt[jj,ji] = slope[jcont]*scaling   # this may be duplicated
               bikeheat[jj,ji] = slope[jcont]*scalingheat   # this may be duplicated
               if totcav == 0:
                  bikemelt[jj,ji] = 0.0
                  bikeheat[jj,ji] = 0.0
    
    print sigslope, totcav, scaling
    
    print bikemelt.shape
    
    print 'sum(bikemelt) = ', np.sum(bikemelt)
    
    bikemelt = bikemelt * self.water_unit_factor
    maskisf=np.where(bikemelt==0,True,False)
    
    #melt=np.ma.array(bikemelt,mask=maskisf)
    #melt = np.ma.transpose(melt)
    self.melt = np.ma.transpose(bikemelt)
    
    bikeheat = bikeheat * self.heat_unit_factor
    #heat=np.ma.array(bikeheat,mask=maskisf)
    #heat = np.ma.transpose(heat)
    self.heat = np.ma.transpose(bikeheat)
    print 'sum(melt) = ', np.sum(self.melt)
    #heat.fill_value=0.
    #heat=heat.filled()*heat_unit_factor
    #melt.fill_value=0.
    #melt=melt.filled()*water_unit_factor
    
    #print "heat[80,330], melt[80,330] = " ,  heat[80,330].fill_value , melt[80,330].fill_value
    
    
    #x_bike=np.load("bike_xcoords.dump")
    #y_bike=np.load("bike_ycoords.dump")
    self.x_bike=np.load("/projects/ukesm/rsmith/ancils/bike_xcoords-AIS.dump")
    self.y_bike=np.load("/projects/ukesm/rsmith/ancils/bike_ycoords-AIS.dump")
  



if __name__ == '__main__':
   filein = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple.nc'
   #filein = '/projects/jmmp/asiaha/TEMP/DataSementara/u-be235/nemo_u-be235o_1y_20050101-20060101_grid-T.nc'
   #fileout = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple-AIS.hdf5'
   fileout='/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/hehe.hdf5'
   filehdfin = '/home/d03/asiaha/cylc-run/u-bk240/bisicles_bh845c_18541201_plot-AIS.hdf5'
   RegObj = SlReg(filein,fileout,filehdfin)
   RegObj.regrid()
   RegObj.makeout()

