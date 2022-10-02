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
from SlReg import SlReg

class CFconsvReg(SlReg):
  def __init__(self,argin,argout,argin_hdf):
    SlReg.__init__(self,argin,argout,argin_hdf)

  def regrid(self):
    CfReg.regrid(self)
    #self.melt_water = self.melt 
    #self.melt_heat =  self.heat
    #self.water_unit_factor = 1.0
    #self.heat_unit_factor = 1.0
    #SlReg.regrid(self,self.melt,self.heat)
    bisicles_to_nemo_mapping_file = '/projects/jomp/asiaha/TEMP/DataSementara/u-bh845/bisicles-AIS_lev0_to_eORCA.map2d'
    bn_map=mp.load(bisicles_to_nemo_mapping_file)
    ny_n=100
    nx_n=np.shape(bn_map.x)[1]

#----------> SlReg.regrid(self,np.ma.transp se(self.melt),np.ma.transpose(self.heat))
    for j in range(ny_n):
      for i in range(nx_n):
        if (np.abs(self.melt_water[j,i]) > 0):
           ncontrib = np.int(bn_map.nmap[j,i])
           meltsub = np.zeros(ncontrib)
           heatsub = np.zeros(ncontrib)
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               meltsub[jcont] = self.melt[ji,jj]
               heatsub[jcont] = self.heat[ji,jj]
           avgmeltsub = np.sum(meltsub)/np.sum(np.abs(meltsub)>1e-10) 
           avgheatsub = np.sum(heatsub)/np.sum(np.abs(heatsub)>1e-10) 
           for jcont in range(ncontrib):
               ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
               self.melt[ji,jj] = self.melt[ji,jj] * self.melt_water[j,i] * self.water_unit_factor/avgmeltsub
               self.heat[ji,jj] = self.heat[ji,jj] * self.melt_heat[j,i] * self.heat_unit_factor/avgheatsub


if __name__ == '__main__':
   filein = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple.nc'
   #filein = '/projects/jmmp/asiaha/TEMP/DataSementara/u-be235/nemo_u-be235o_1y_20050101-20060101_grid-T.nc'
   #fileout = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple-AIS.hdf5'
   fileout='/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/CFconsv.hdf5'
   filehdfin = '/home/d03/asiaha/cylc-run/u-bk240/bisicles_bh845c_18541201_plot-AIS.hdf5'
   RegObj = CFconsvReg(filein,fileout,filehdfin)
   RegObj.regrid()
   RegObj.makeout()
