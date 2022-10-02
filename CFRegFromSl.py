import os
import numpy as np
import sys
import cf
from netCDF4 import Dataset

class CfReg(object):

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
  
  
    self.regrid_ncfile  =os.path.splitext(os.path.basename(self.regrid_hdf5file))[0]+".nc"
   
    # nanti utk slope extrapolation harus beda 
    #self.bike_ncgridfile="cf_bikegridfile.nc"
    #self.nemo_ncgridfile="cf_gridfile.nc"
    self.bike_ncgridfile="/projects/ukesm/rsmith/ancils/cf_gridfile_BISICLES_lev0-AIS.nc" # 768 x 768
    self.nemo_ncgridfile="/projects/ukesm/rsmith/ancils/cf_gridfile_eORCA1_v2.2x.nc" # 362 x 332
  
    #self.water_unit_factor=60*60*24*360/1e3
    #self.heat_unit_factor=60*60*24*360
    self.water_unit_factor=1.0
    self.heat_unit_factor=1.0
 
    print 'self.nctgridfile = ', self.nctgridfile 
    h=cf.read(self.nctgridfile)
  
    #self.melt_water=h.select('melt_water')[0].array.squeeze()
    #self.melt_heat =h.select('melt_heat')[0].array.squeeze()
  
    self.melt_water=Dataset(self.nctgridfile).variables['melt_water'][:]
    self.melt_heat=Dataset(self.nctgridfile).variables['melt_heat'][:]

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
  
    self.regrid()
  
    #create output
    self.makeout()
  
  
  def regrid(self):
    g=cf.read(self.bike_ncgridfile)[0]
    print 'self.bike_ncgridfile =', self.bike_ncgridfile 
    print " regrid 1"
    #f=cf.read(self.nemo_ncgridfile)[0]
    f=cf.read(self.bike_ncgridfile)[0]
    print 'self.nemo_ncgridfile = ', self.nemo_ncgridfile
    f.insert_data(cf.Data(self.melt_water,units='kg/m2/s'), axes=('Y', 'X'))
    w_regrid=f.regrids(g,src_cyclic=True,method='bilinear')
    #w_regrid=g
  
    print " regrid 2"
    #f=cf.read(self.nemo_ncgridfile)[0]
    f=cf.read(self.bike_ncgridfile)[0]
    f.insert_data(cf.Data(self.melt_heat,units='kg/m2/s'), axes=('Y', 'X'))
    h_regrid=f.regrids(g,src_cyclic=True,method='bilinear')
    #h_regrid=g
  
  
    #i=cf.FieldList()
    #i.append(w_regrid)
    #i.append(h_regrid)
    #cf.write(i,regrid_ncfile)
  
    #like BIKEtoNEMO needs a rollaxis (and I don't get why)
    #we need to switch x and y here to get it the way BISICLES wants it
    self.heat=np.swapaxes(h_regrid.array,0,1)
    self.melt=np.swapaxes(w_regrid.array,0,1)
  
    #self.heat.fill_value=0.
    #self.heat=self.heat.filled()*self.heat_unit_factor
    #self.melt.fill_value=0.
    #self.melt=self.melt.filled()*self.water_unit_factor
  
    #self.x_bike=np.load("bike_xcoords.dump")
    #self.y_bike=np.load("bike_ycoords.dump")
  
    self.x_bike=np.load("/projects/ukesm/rsmith/ancils/bike_xcoords-AIS.dump")
    self.y_bike=np.load("/projects/ukesm/rsmith/ancils/bike_ycoords-AIS.dump")
  
  
  def makeout(self):
    from netCDF4 import Dataset
  #add Steph's calving-edge-enforced-by-melt criterion - how now?
  #melt[np.where()]=-1.0e+4
  
    ncfile_out = Dataset(self.regrid_ncfile,'w',format='NETCDF3_CLASSIC')
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

if __name__ == '__main__':
   filein = 'slopeReg.nc'
   #filein = '/projects/jmmp/asiaha/TEMP/DataSementara/u-be235/nemo_u-be235o_1y_20050101-20060101_grid-T.nc'
   #fileout = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple-AIS.hdf5'
   fileout='slopeRegtoCF.hdf5'
   RegObj = CfReg(filein,fileout)
   RegObj.regrid()
   RegObj.makeout()


