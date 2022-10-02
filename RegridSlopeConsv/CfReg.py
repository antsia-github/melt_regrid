import os
import numpy as np
import sys
import cf
from Reg import Reg
import mapping_class as mp

class CfReg(Reg):

  def regrid(self):
    g=cf.read(self.bike_ncgridfile)[0]
    print 'self.bike_ncgridfile =', self.bike_ncgridfile
    print " regrid 1"
    f=cf.read(self.nemo_ncgridfile)[0]
    print 'self.nemo_ncgridfile = ', self.nemo_ncgridfile
    f.insert_data(cf.Data(self.melt_water,units='kg/m2/s'), axes=('Y', 'X'))
    w_regrid=f.regrids(g,src_cyclic=True,method='bilinear')
    #w_regrid=g

    print " regrid 2"
    f=cf.read(self.nemo_ncgridfile)[0]
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

    self.heat.fill_value=0.
    self.heat=self.heat.filled()*self.heat_unit_factor
    self.melt.fill_value=0.
    self.melt=self.melt.filled()*self.water_unit_factor

    #self.x_bike=np.load("bike_xcoords.dump")
    #self.y_bike=np.load("bike_ycoords.dump")

    #self.x_bike=np.load("/projects/ukesm/rsmith/ancils/bike_xcoords-AIS.dump")
    #self.y_bike=np.load("/projects/ukesm/rsmith/ancils/bike_ycoords-AIS.dump")


if __name__ == '__main__':
   filein = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple.nc'
   #filein = '/projects/jmmp/asiaha/TEMP/DataSementara/u-be235/nemo_u-be235o_1y_20050101-20060101_grid-T.nc'
   #fileout = '/projects/jmmp/asiaha/TEMP/FileCobaNemo/u-bl374/nemo_bl374c_P1Y_20081201-20091201_icecouple-AIS.hdf5'
   fileout='hehe.hdf5'
   filehdfin = '/home/d03/asiaha/cylc-run/u-bk240/bisicles_bh845c_18541201_plot-AIS.hdf5'
   RegObj = CfReg(filein,fileout)
   RegObj.regrid()
   #------making it constant, cancelling the above regridding---------
   #RegObj.melt = np.where(np.abs(RegObj.melt)> 1e-10,-10.0,RegObj.melt)

   RegObj.hdfplotfile = filehdfin;   RegObj.adjust_slope()
   RegObj.consv(hdfile=filehdfin)
   RegObj.makeout()

