import os
import numpy as np
import sys
import cf
from Reg import Reg
import mapping_class as mp

class BinReg(Reg):
  def __init__(self,argin,argout,argin_hdf):
    #err = 0
    #self.hdfplotfile, err  = arg_to_file_exist(argin_hdf, mandatory=True, err=err)
    #if err > 0:
    #  print 'Error opening ', argin_hdf
    #  sys.exit(2)
    if not os.path.exists(argin_hdf):
      print 'Error opening ', argin_hdf
      sys.exit(2)
    super(BinReg,self).__init__(argin,argout)
    self.hdfplotfile = argin_hdf
    
 
  def regrid(self,mapfile=None):
    if mapfile is None:
      mapfile = '/projects/jomp/asiaha/TEMP/DataSementara/u-bh845/bisicles-AIS_lev0_to_eORCA.map2d'

      #level=0 # level of grid refinement
      #order=1 # interpolation order, 0 for piecewise constant, 1 for linear
    l_verb = True

    #if (l_verb): print 'reading hdf file'
    #amrID = amrio.load(self.hdfplotfile)
    #lo,hi = amrio.queryDomainCorners(amrID, level)
    #_,_,isf = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
    #_,_,bathy = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)
    #bathy=np.transpose(bathy)
    #isf=np.transpose(isf)
    if self.isf is None:
      self.gen_bathy_isf()  

    bn_map=mp.load(mapfile) 
    ny_n=100
    #ny_n=np.shape(bn_map.x)[0]
    nx_n=np.shape(bn_map.x)[1]  
    bikemelt = np.zeros_like(self.isf)
    bikeheat = np.zeros_like(self.isf)


    if (l_verb): print "remap geometry to NEMO"
    for j in range(ny_n):
      if (l_verb): print '\r',j,'  ',
      sys.stdout.flush()
      print 'j, sum(bikemelt) = ', j, np.sum(bikemelt)
      for i in range(nx_n):
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]        
               if (np.abs(self.isf[jj,ji] - self.bathy[jj,ji])  > 1E-3) \
               & (np.abs(self.isf[jj,ji]) > 1E-3):
                  bikemelt[jj,ji] = self.melt_water[j,i]*self.water_unit_factor #not yet transposed
                  bikeheat[jj,ji] = self.melt_heat[j,i]*self.heat_unit_factor # not yet transposed
    #we need to switch x and y here to get it the way BISICLES wants it, such as self.heat=np.swapaxes(h_regrid.array,0,1) in CF regridding
    self.melt = np.ma.transpose(bikemelt) # sejak disini sudah transpose 
    self.heat = np.ma.transpose(bikeheat)



  #def consv(self,mapfile=None):
  #    print 'This binning type always satisfy conservation.'


if __name__ == '__main__':
   filein = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple.nc'
   #filein = '/projects/jmmp/asiaha/TEMP/DataSementara/u-be235/nemo_u-be235o_1y_20050101-20060101_grid-T.nc'
   #fileout = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple-AIS.hdf5'
   fileout='hehe.hdf5'
   filehdfin = '/home/d03/asiaha/cylc-run/u-bk240/bisicles_bh845c_18541201_plot-AIS.hdf5'
   RegObj = BinReg(filein,fileout,filehdfin)
   RegObj.regrid()
   #------making it constant, cancelling the above regridding---------
   RegObj.melt = np.where(np.abs(RegObj.melt)> 1e-10,-10.0,RegObj.melt)
   RegObj.adjust_slope()
   #RegObj.consv()
   RegObj.makeout()

