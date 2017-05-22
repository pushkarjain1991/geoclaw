__author__ = 'Pushkar Kumar Jain'

import pandas as pd
import numpy as np


class ReadAmr2(object):

    def __init__(self, filename, sort_in_frame=True):
        #self.filename = filename
        #self.sort_in_frame = sort_in_frame
        self.all_data = pd.read_table(filename, header=None, names = ["height","xvel","yvel","eta"], index_col=False, sep=r"\s+")
        self.Grid_level = self.capture_data("grid_number")
        self.AMR_level = self.capture_data("AMR_level")
        self.mx = self.capture_data("mx")
        self.my = self.capture_data("my")
        self.x_low = self.capture_data("xlow")
        self.y_low = self.capture_data("ylow")
        self.dx = self.capture_data("dx")
        self.dy = self.capture_data("dy")
        self.pandas_dataframe = self.amrdataframe(sort_in_frame)


    def get_mycolumn(self, column,amrl):
        mycolumn_data = self.pandas_dataframe[column][self.pandas_dataframe.amrlevel==amrl]
        return mycolumn_data


    def capture_data(self, data_string):
        out = self.all_data[self.all_data.xvel==data_string].height.values
        #print data_string, out
        if(data_string in ["grid_number", "AMR_level", "mx", "my"]):
            out = out.astype('int32')
           
        return out

    def amrdataframe(self, sort_in_frame):
        # Read all levels grid
        data = self.all_data.dropna()
        data = data.reset_index(drop=True)

        for num,levelnum in enumerate(self.AMR_level):
            if num == 0:
                firstpoint = 0
            else:
                firstpoint = firstpoint + self.mx[num-1]*self.my[num-1]
            secondpoint = firstpoint + self.mx[num]*self.my[num]-1
            data.loc[firstpoint:secondpoint, 'amrlevel'] = levelnum 
            x_left = self.x_low[num] + self.dx[num]/2.0
            x_right = self.x_low[num] + self.dx[num]*self.mx[num] - self.dx[num]/2.0
            y_down = self.y_low[num] + self.dy[num]/2.0
            y_up = self.y_low[num] + self.dy[num]*self.my[num] - self.dy[num]/2.0
            xrow = np.linspace(x_left , x_right, num = self.mx[num], dtype='float64')
            yrow = np.linspace(y_down , y_up, num = self.my[num], dtype='float64')
            xmesh,ymesh = np.meshgrid(xrow, yrow)
            data.loc[firstpoint:secondpoint,'xcoord']=np.ravel(xmesh)
            data.loc[firstpoint:secondpoint,'ycoord']=np.ravel(ymesh)

        #Rearranging xcoord and ycoord as per domain mesh
        if sort_in_frame:
            data.sort_values(by=['amrlevel','ycoord', 'xcoord'], ascending=[True,True,True],inplace=True)
            data = data.reset_index(drop=True)
        return data

     
if __name__=="__main__":
    #hello = ReadAmr("./ens_1_1/fort.q0012")
    hello = ReadAmr2("../_output/fort.q0001", sort_in_frame=False)
    yoyo = hello.pandas_dataframe
    #print yoyo[(yoyo.xcoord==47.0) & (yoyo.ycoord==5.0)]
    #print yoyo[(yoyo.xcoord < 47.01) & (yoyo.xcoord > 46.9)]
    #print yoyo[(yoyo.xcoord == 47.00)]
    #print yoyo[np.isclose(yoyo.xcoord, 47.0)]
    #print yoyo[(yoyo.ycoord==5.0)]
    #print yoyo.dtypes
    print yoyo
    #yoyo.to_csv("hello", sep='\t')
    #print yoyo.xcoord
