__author__ = 'Pushkar Kumar Jain'

import pandas as pd
import numpy as np
import read_amr


class ReadAmrForLevel(read_amr.ReadAmr):
    def __init__(self, filename, amrl):
        super(ReadAmrForLevel,self).__init__(filename)
        self.amrl = amrl
        self.mxv, self.myv = self.get_meshgrid()
        
        #Read raw columns
        self.total_height = self.get_mycolumn("height", amrlevel=amrl)
        self.momx = self.get_mycolumn("xvel", amrlevel=amrl)
        self.eta = self.get_mycolumn("eta", amrlevel=amrl)

        #Get eta with value in land = 0
        self.eta_with_land = self.get_eta_with_land0()
        #Get land with water masked
        self.land = self.get_land()
        #Get water with land masked
        self.water = self.get_water()

    def get_meshgrid(self):
        xlower = self.x_low[0]
        ylower = self.y_low[0]
        dx = self.dx[0]
        dy = self.dy[0]
        mx = self.mx[0]
        my = self.my[0]

        xupper = xlower + (mx+1)*dx
        yupper = ylower + (my+1)*dy

        x_cell = np.linspace(xlower + dx/2.0, xupper - dx/2.0, mx)
        y_cell = np.linspace(yupper - dy/2.0, ylower + dy/2.0, my)
        mxv, myv = np.meshgrid(x_cell, y_cell)
        return(mxv, myv)


    def get_land(self):
        masked_eta_land = np.ma.array(self.total_height, mask=self.momx==0.0E0)
        reshaped_masked_eta_land = np.reshape(masked_eta_land,(self.mx[0], self.my[0]))
        return reshaped_masked_eta_land

    def get_water(self):
        masked_eta_water = np.ma.array(self.eta, mask=self.total_height==0.0E0)
        reshaped_masked_eta_water = np.reshape(masked_eta_water,(self.mx[0], self.my[0]))
        return reshaped_masked_eta_water

    def get_mycolumn(self, column, amrlevel=1.0):
        mycolumn_data = self.pandas_dataframe[column][self.pandas_dataframe.amrlevel==self.amrl]
        return mycolumn_data
    
    def get_eta_with_land0(self):
        eta2 = self.eta
        eta2[self.total_height == 0.0] = 0.0
        #eta_with_land_matrix = eta2.as_matrix()
        eta_with_land_matrix = eta2.values
        reshaped_eta_with_land0 = np.reshape(eta_with_land_matrix, (self.mx[0], self.my[0]))
        return reshaped_eta_with_land0


if __name__ == "__main__":
    read_verification_output = "./ens_1_0/fort.q0012"
    test_case = ReadAmrForLevel(read_verification_output,1)
    print test_case.mx[0]


