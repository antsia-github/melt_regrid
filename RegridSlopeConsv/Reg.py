import os
import numpy as np
import sys
import cf
from netCDF4 import Dataset
from amrfile import io as amrio
import mapping_class as mp
import math

Clust = {}
#Clust = {'PIG':[(52,186), (52,187), (53,186), (53,187), (54,186), (54,187), (55,186), (55,187), (53,188)], 'Thwaites':[(53,180),(53,181), (53,182), (52,183), (54,180), (54,181), (54,182)]}

class Reg(object):
  def __init__(self,argin, argout):
    good_files=0
    print 'argin = ', argin
    for file in [argin]:
      print 'file = ', file
      if file:
        if  os.path.isfile(file):
          good_files+=1
        else:
            print ""
            print "ERROR: specified file does not exist:",file

    if good_files==1:
      self.nctgridfile=argin
    else:
      print ""
      print "ERROR: problem with one or more input files. I want ALL of them"
      print ""
      parser.print_help()
      sys.exit(2)

    if argout != None:
      self.regrid_hdf5file=argout
    else:
      print ""
      print "ERROR: explicitly specifiy an output file"
      print ""
      parser.print_help()
      sys.exit(2)

    self.isf = None; self.bathy = None; self.latb = None; self.lonb = None # defined later if necessary  
    self.hdfplotfile = None
    self.regrid_ncfile  =os.path.splitext(os.path.basename(self.regrid_hdf5file))[0]+".nc"

    # nanti utk slope extrapolation harus beda 
    #self.bike_ncgridfile="cf_bikegridfile.nc"
    #self.nemo_ncgridfile="cf_gridfile.nc"
    self.bike_ncgridfile="/projects/ukesm/rsmith/ancils/cf_gridfile_BISICLES_lev0-AIS.nc" # 768 x 768  ===> input
    self.nemo_ncgridfile="/projects/ukesm/rsmith/ancils/cf_gridfile_eORCA1_v2.2x.nc" # 362 x 332  ==> input 

    self.water_unit_factor=60*60*24*360/1e3 #==> input
    self.heat_unit_factor=60*60*24*360  # ==> input 

    self.x_bike=np.load("/projects/ukesm/rsmith/ancils/bike_xcoords-AIS.dump") #==> input 
    self.y_bike=np.load("/projects/ukesm/rsmith/ancils/bike_ycoords-AIS.dump") #==> input

    print 'self.nctgridfile = ', self.nctgridfile
    h=cf.read(self.nctgridfile)

    self.melt_water=h.select('long_name:Ice shelf melting')[0].array.squeeze()
    self.melt_heat =h.select('long_name:Ice shelf heat content flux')[0].array.squeeze()

    #self.melt_water=np.ma.masked_equal(self.melt_water,0)
    #self.melt_heat=np.ma.masked_equal(self.melt_heat,0)
    #need to do better than this - diagnostics on ORCA grids may have junk in the halos/overlap points, it seems

    print 'dim shape = ', np.shape(self.melt_water), self.melt_water.ndim
    #time average
    if self.melt_water.ndim==3:
      print "Doing time average",np.shape(self.melt_water)
      self.melt_water=np.mean(self.melt_water,axis=0)
      self.melt_heat =np.mean(self.melt_heat,axis=0)
      print np.shape(self.melt_water)


  def gen_bathy_isf(self):
     if self.hdfplotfile is None:
        print 'Error, choose the right file '

     level=0 # level of grid refinement
     order=1 # interpolation order, 0 for piecewise constant, 1 for linear
     print 'reading hdf file'
     amrID = amrio.load(self.hdfplotfile) #--this needs to be defined beforehands
     lo,hi = amrio.queryDomainCorners(amrID, level)
     _,_,self.isf = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
     _,_,self.bathy = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)
     self.bathy=np.transpose(self.bathy)
     self.isf=np.transpose(self.isf)

  def gen_latlon_bike(self):
    self.latb = Dataset('/projects/ukesm/rsmith/ancils/cf_gridfile_BISICLES_lev0-AIS.nc').variables['latitude'][:]
    self.lonb = Dataset('/projects/ukesm/rsmith/ancils/cf_gridfile_BISICLES_lev0-AIS.nc').variables['longitude'][:]


  def adjust_slope(self, mapfile=None, meltpat=None, heatpat=None, l_verb=True):
    if mapfile is None:
      mapfile = '/projects/jomp/asiaha/TEMP/DataSementara/u-bh845/bisicles-AIS_lev0_to_eORCA.map2d'

    if self.isf is None:
      self.gen_bathy_isf()
    if self.latb is None:
      self.gen_latlon_bike() #generate self.latb and self.lonb

    if meltpat is None:
      bikemelt = np.ma.transpose(self.melt) # transpose because of the 2D mapping
      bikeheat = np.ma.transpose(self.heat)
    else:
      bikemelt = np.ma.transpose(meltpat) # does it need to be multiplied with unit factor ?
      bikeheat = np.ma.transpose(heatpat)

    bn_map=mp.load(mapfile)
    ny_n=100
    #ny_n=np.shape(bn_map.x)[0]
    nx_n=np.shape(bn_map.x)[1]

    AllClustCell = [] # container for all NEMO cells we group in CLust dictionary 
    for idr in Clust.keys(): #for every region
      AllClustCell.append(Clust[idr])

    if (l_verb): print "remap geometry to NEMO for slope adjustment"
    for j in range(ny_n):
      if (l_verb): print '\r',j,'  ',
      sys.stdout.flush()
      print 'j, sum(bikemelt) = ', j, np.sum(bikemelt)
      for i in range(nx_n):
       if (j,i) not in AllClustCell:
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           slope = np.zeros(ncontrib)
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               if (np.abs(self.isf[jj,ji] - self.bathy[jj,ji])  > 1E-3) \
               & (np.abs(self.isf[jj,ji]) > 1E-3):
                    #ingat self.meltmelt[jj,ji] dan self.heat[jj,ji] sudah di-transpose
                   #-----here is the place to edit for different slope technique ----  
                   slplon = self.gradslplon(ji,jj)
                   slplat = self.gradslplat(ji,jj)
                   #---------------------------------- 
                   slope[jcont] = math.sqrt(slplon**2 + slplat**2)
           sigslope = np.sum(slope); totcav = np.sum(slope>0.0)
           avgslope = sigslope/totcav

           if math.isnan(avgslope):
              print 'i, j, ncontrib = ', i, j, ncontrib, totcav, sigslope, self.melt_water[j,i]

           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               bikemelt[jj,ji] = slope[jcont]/avgslope * bikemelt[jj,ji] #will make grounded meltrate into zero, as its slope is zero
               bikeheat[jj,ji] = slope[jcont]/avgslope * bikeheat[jj,ji]
               if totcav == 0: #----sepertinya ini dibutuhkan utk mencegah NaN di dalam NetCDF file
                  bikemelt[jj,ji] = 0.0  ;  bikeheat[jj,ji] = 0.0
    #Transpose back into the original form
    for idr in Clust.keys(): #for every region    
      ListNemoCell = Clust[idr]
      ListSlope = []
      for (j,i) in ListNemoCell:
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               if (np.abs(self.isf[jj,ji] - self.bathy[jj,ji])  > 1E-3) \
               & (np.abs(self.isf[jj,ji]) > 1E-3):
                    #ingat self.meltmelt[jj,ji] dan self.heat[jj,ji] sudah di-transpose
                   #-----here is the place to edit for different slope technique ----  
                   slplon = self.gradslplon(ji,jj)
                   slplat = self.gradslplat(ji,jj)
                   #---------------------------------- 
                   ListSlope.append(math.sqrt(slplon**2 + slplat**2))
               else:
                   ListSlope.append(0.0)
      print 'Collecting slope idr, len = ', idr, len(ListSlope)
      slope = np.array(ListSlope)
      sigslope = np.sum(slope); totcav = np.sum(slope>0.0)
      avgslope = sigslope/totcav

      if math.isnan(avgslope):
          print 'i, j, ncontrib = ', i, j, ncontrib, totcav, sigslope, self.melt_water[j,i]

      totcont = 0 
      for (j,i) in ListNemoCell:
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               bikemelt[jj,ji] = slope[totcont]/avgslope * bikemelt[jj,ji]
               bikeheat[jj,ji] = slope[totcont]/avgslope * bikeheat[jj,ji]
               if totcav == 0: #---may be necessary to avoid NaN
                  bikemelt[jj,ji] = 0.0  ;  bikeheat[jj,ji] = 0.0
               totcont = totcont + 1  
      print 'Adjusting the slope idr, totcont = ', idr, totcont
    self.melt = np.ma.transpose(bikemelt)
    self.heat = np.ma.transpose(bikeheat)


  def consv(self,mapfile=None,cluster=None,hdfile=None): 
    #Clust = {'PIG':[(52,186), (52,187), (53,186), (53,187), (54,186), (54,187), (55,186), (55,187), (53,188)], 'Thwaites':[(53,180),(53,181), (53,182), (52,183), (54,180), (54,181), (54,182)]}
    if mapfile is None:
      mapfile = '/projects/jomp/asiaha/TEMP/DataSementara/u-bh845/bisicles-AIS_lev0_to_eORCA.map2d'
    if hdfile is not None: #--needed by CF regridding. Slope adjustment already provided one. 
      self.hdfplotfile = hdfile 
    if self.isf is None:
      self.gen_bathy_isf()

    e1t=Dataset('/projects/ukesm/rsmith/NEMO_RSS/mesh_mask-be235.nc').variables['e1t'][0]
    e2t=Dataset('/projects/ukesm/rsmith/NEMO_RSS/mesh_mask-be235.nc').variables['e2t'][0]
    area = e1t*e2t
    AllClustCell = [] # container for all NEMO cells we group in CLust dictionary
    for idr in Clust.keys(): #for every region
      AllClustCell.append(Clust[idr])

    #bikemelt = np.ma.transpose(self.melt) # transpose because of the 2D mapping
    #bikeheat = np.ma.transpose(self.heat)

    bn_map=mp.load(mapfile)
    ny_n=100
    #ny_n=np.shape(bn_map.x)[0]
    nx_n=np.shape(bn_map.x)[1]
#----------> SlReg.regrid(self,np.ma.transp se(self.melt),np.ma.transpose(self.heat))
    for j in range(ny_n):
      for i in range(nx_n):
       if (j,i) not in AllClustCell:
        if (np.abs(self.melt_water[j,i]) > 0): 
           ncontrib = np.int(bn_map.nmap[j,i])
           meltsub = np.zeros(ncontrib)
           heatsub = np.zeros(ncontrib)
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               if (np.abs(self.isf[jj,ji] - self.bathy[jj,ji])  > 1E-3) \
               & (np.abs(self.isf[jj,ji]) > 1E-3):   # rememer, isf is the transpose of melt
                 meltsub[jcont] = self.melt[ji,jj]   #---using [ji,jj] as it is not transposed 
                 heatsub[jcont] = self.heat[ji,jj]
               else:  #make it zero if set wrongly in the regridding stage
                 self.melt[ji,jj]= 0.0;  self.heat[ji,jj] = 0.0 
           totcav = np.sum(np.abs(meltsub)>1e-10)
           if  totcav< 1: #----to avoid NaN, sometimes the mapping 2D is not accurate ?? A Nemo cell with all grounded BIKE cell ??
               continue 
           #   print 'total meltsub = 0.0, i, j, ji,jj = ', i, j, ji, jj 
           avgmeltsub = np.sum(meltsub)/np.sum(np.abs(meltsub)>1e-10)
           avgheatsub = np.sum(heatsub)/np.sum(np.abs(heatsub)>1e-10)
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
 #-----> perhatikan kalo belon di-transpose 
               self.melt[ji,jj] = self.melt[ji,jj] * self.melt_water[j,i] * self.water_unit_factor/avgmeltsub
               self.heat[ji,jj] = self.heat[ji,jj] * self.melt_heat[j,i] * self.heat_unit_factor/avgheatsub

    for idr in Clust.keys(): #for every region
      ListNemoCell = Clust[idr]
      ListMelt = []; ListHeat=[]
      totarea = 0.0; meltxarea = 0.0; heatxarea = 0.0
      for (j,i) in ListNemoCell:
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               if (np.abs(self.isf[jj,ji] - self.bathy[jj,ji])  > 1E-3) \
               & (np.abs(self.isf[jj,ji]) > 1E-3):   # remember, isf is the transpose of melt
                 ListMelt.append(self.melt[ji,jj])
                 ListHeat.append(self.heat[ji,jj])
               else:
                 self.melt[ji,jj]= 0.0;  self.heat[ji,jj] = 0.0
           meltxarea = meltxarea + area[j,i]*self.melt_water[j,i]*self.water_unit_factor
           heatxarea = heatxarea + area[j,i]*self.melt_heat[j,i]*self.heat_unit_factor
           totarea = totarea + area[j,i]
      meltsub = np.array(ListMelt)
      heatsub = np.array(ListHeat)
      totcav = np.sum(np.abs(meltsub)>1e-10)
      if  totcav< 1: #----to avoid NaN, sometimes the mapping 2D is not accurate ?? A Nemo cell with all grounded BIKE cell ??
          continue #the next group if no cavity BIKEin that group. Very unlikely for a group, though ! 
      avgmeltsub = np.sum(meltsub)/np.sum(np.abs(meltsub)>1e-10)
      avgheatsub = np.sum(heatsub)/np.sum(np.abs(heatsub)>1e-10)
      avgnemomelt = meltxarea/totarea
      avgnemoheat = heatxarea/totarea
      print 'Averaging the melt idr, totcont = ', idr, len(ListMelt)

#           coefconsv = sigslope*self.melt_water[j,i]*self.water_unit_factor/sigslopexmelt 
      totcont = 0
      for (j,i) in ListNemoCell:
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               self.melt[ji,jj] = self.melt[ji,jj] * avgnemomelt/avgmeltsub
               self.heat[ji,jj] = self.heat[ji,jj] * avgnemoheat/avgheatsub
               totcont = totcont + 1
      print 'Conserving the melt idr, totcont = ', idr, totcont

  def gradslplon(self,ji,jj):
    slplon = []
    if (np.abs(self.isf[jj,ji+1] - self.bathy[jj,ji+1])  > 1E-3) & (np.abs(self.isf[jj,ji+1]) > 1E-3) :
       slplonpl = np.abs(self.isf[jj,ji] - self.isf[jj,ji+1])/self.Distance(self.latb[jj,ji],self.lonb[jj,ji],self.latb[jj,ji+1],self.lonb[jj,ji+1])
       slplon.append(slplonpl)
    if (np.abs(self.isf[jj,ji-1] - self.bathy[jj,ji-1])  > 1E-3) & (np.abs(self.isf[jj,ji-1]) > 1E-3) :
       slplonne = np.abs(self.isf[jj,ji] - self.isf[jj,ji-1])/self.Distance(self.latb[jj,ji],self.lonb[jj,ji],self.latb[jj,ji-1],self.lonb[jj,ji-1])
       slplon.append(slplonne)
    if len(slplon) > 0:
       slplon = sum(slplon)/len(slplon)
    else:
       slplon = 0.0
    return slplon

  def gradslplat(self,ji,jj):
    slplat = []
    if (np.abs(self.isf[jj+1,ji] - self.bathy[jj+1,ji])  > 1E-3) & (np.abs(self.isf[jj+1,ji]) > 1E-3) :
       slplatpl = np.abs(self.isf[jj,ji] - self.isf[jj+1,ji])/self.Distance(self.latb[jj,ji],self.lonb[jj,ji],self.latb[jj+1,ji],self.lonb[jj+1,ji])
       slplat.append(slplatpl)
    if (np.abs(self.isf[jj-1,ji] - self.bathy[jj-1,ji])  > 1E-3) & (np.abs(self.isf[jj-1,ji]) > 1E-3) :
       slplatne = np.abs(self.isf[jj,ji] - self.isf[jj-1,ji])/self.Distance(self.latb[jj,ji],self.lonb[jj,ji],self.latb[jj-1,ji],self.lonb[jj-1,ji])
       slplat.append(slplatne)
    if len(slplat) > 0:
       slplat = sum(slplat)/len(slplat)
    else:
       slplat = 0.0
    return slplat

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



  def makeout(self, ncfile=None):
    from netCDF4 import Dataset
  #add Steph's calving-edge-enforced-by-melt criterion - how now?
  #melt[np.where()]=-1.0e+4

    if ncfile is None:
      ncfile_out = Dataset(self.regrid_ncfile,'w',format='NETCDF3_CLASSIC')
    else:
      ncfile_out = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')
    ncfile_out.createDimension('x',self.x_bike.shape[0])
    ncfile_out.createDimension('y',self.x_bike.shape[1])

    x_nc=ncfile_out.createVariable('x',np.dtype('float64').char,('x'))
    y_nc=ncfile_out.createVariable('y',np.dtype('float64').char,('y'))

    x_nc[:]=self.x_bike[0,:]
    y_nc[:]=self.y_bike[:,0]

    heat_nc=ncfile_out.createVariable('melt_heat',np.dtype('float64').char,('y','x'))
    melt_nc=ncfile_out.createVariable('melt_water',np.dtype('float64').char,('y','x'))
    heat_nc[:]=self.heat
    melt_nc[:]=self.melt

    ncfile_out.close()

    print " calling nctoamr"
    #os.system('nctoamr2d.ex %(self.regrid_ncfile)s %(self.regrid_hdf5file)s self.melt_water self.melt_heat' % locals())
    regrid_ncfile = self.regrid_ncfile ; regrid_hdf5file = self.regrid_hdf5file

    #os.system('nctoamr2d.ex %(regrid_ncfile)s %(regrid_hdf5file)s melt_water melt_heat' % locals())
    os.system('/projects/ukesm/rsmith/unicicles_rss_ukesm/BISICLES/code/filetools/nctoamr2d.Linux.64.CC.ftn.OPT.INTEL.ex %(regrid_ncfile)s %(regrid_hdf5file)s melt_water melt_heat' % locals())




