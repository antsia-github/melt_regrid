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

class CfRegConst(CfReg):
  #def __init__(self,argin,argout,argin_hdf):
  #  SlReg.__init__(self,argin,argout,argin_hdf)

  def regrid(self):
    self.melt_water = np.where(np.abs(self.melt_water)> 1e-10,-10.0/self.water_unit_factor,self.melt_water)
    CfReg.regrid(self)
    #self.melt_water = self.melt 
    #self.melt_heat =  self.heat
    #self.water_unit_factor = 1.0
    #self.heat_unit_factor = 1.0
    #SlReg.regrid(self,self.melt,self.heat)
    #self.melt_water = np.where(np.abs(self.melt_water)> 1e-10,-10.0/self.water_unit_factor,self.melt_water)
    #SlReg.regrid(self)
    #SlReg.regrid(self,np.ma.transpose(self.melt),np.ma.transpose(self.heat))



if __name__ == '__main__':
   filein = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple.nc'
   #filein = '/projects/jmmp/asiaha/TEMP/DataSementara/u-be235/nemo_u-be235o_1y_20050101-20060101_grid-T.nc'
   #fileout = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple-AIS.hdf5'
   fileout='/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/CfRegConst.hdf5'
   RegObj = CfRegConst(filein,fileout)
   RegObj.regrid()
   RegObj.makeout()
