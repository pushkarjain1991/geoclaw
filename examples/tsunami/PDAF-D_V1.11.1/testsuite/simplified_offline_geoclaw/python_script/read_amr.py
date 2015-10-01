__author__ = 'Pushkar Kumar Jain'

import pandas as pd
import numpy as np


class ReadAmr(object):
    def __init__(self, filename):
        self.filename = filename
        self.Grid_level, self.Grid_lines = self.capture_data("grid_number")
        self.AMR_level, self.AMR_lines = self.capture_data("AMR_level")
        self.mx, self.mx_lines = self.capture_data("mx")
        self.my, self.my_lines = self.capture_data("my")
        self.x_low, self.x_low_lines = self.capture_data("xlow")
        self.y_low, self.y_low_lines = self.capture_data("ylow")
        self.dx, self.dx_lines = self.capture_data("dx")
        self.dy, self.dy_lines = self.capture_data("dy")
        self.pandas_dataframe = self.amrdataframe()

    def get_mycolumn(self, column,amrl):
        mycolumn_data = self.pandas_dataframe[column][self.pandas_dataframe.amrlevel==amrl]
        return mycolumn_data

    #def get_land(self):
    #    masked_eta_land = np.ma.array(

    def capture_data(self, data_string):
        value_list = []
        line_number_list = []
        with open(self.filename, "r") as f1:
            for unwanted_line_number, j in enumerate(f1.readlines()):
                if data_string in j:
                    value = j.split()
                    if value[0].isdigit():
                        value_list.append(int(value[0]))
                    else:
                        value_list.append(float(value[0]))
                    line_number_list.append(unwanted_line_number)
        return value_list, line_number_list

    def amrdataframe(self):
        # Read all levels grid
        data = pd.read_table(self.filename, header=None, names = ["height","xvel","yvel","eta"], index_col=False, sep=r"\s+")
        data = data.dropna()
        data=data.reset_index(drop=True)
        xmain = np.array([])
        ymain = np.array([])

        for num,levelnum in enumerate(self.AMR_level):
            if num == 0:
                firstpoint = 0
            else:
                firstpoint = firstpoint + self.mx[num-1]*self.my[num-1]
            secondpoint = firstpoint + self.mx[num]*self.my[num]
            #data['amrlevel'] = levelnum 
            data.loc[firstpoint:secondpoint, 'amrlevel'] = levelnum 
            x_left = self.x_low[num] + self.dx[num]/2.0
            x_right = self.x_low[num] + self.dx[num]*self.mx[num] - self.dx[num]/2.0
            y_down = self.y_low[num] + self.dy[num]/2.0
            y_up = self.y_low[num] + self.dy[num]*self.my[num] - self.dy[num]/2.0
            xrow = np.linspace(x_left , x_right, num = self.mx[num])
            yrow = np.linspace(y_down , y_up, num = self.my[num])
            #yrow = np.linspace( self.y_low[num] + self.dy[num]*self.my[num], self.y_low[num], num = self.my[num]+1)
            xmesh,ymesh = np.meshgrid(xrow,yrow)
            #flatx = np.ravel(xmesh)
            #flaty = np.ravel(ymesh)
            #print np.shape(flatx)
            #data.insert(firstpoint, 'xcoord',np.ravel(xmesh))
            xmain = np.append(xmain, np.ravel(xmesh))
            ymain = np.append(ymain, np.ravel(ymesh))
            #xseries = pd.Series(np.ravel(xmesh), name='xcoord')
            #yseries = pd.Series(np.ravel(ymesh), name='ycoord')
            #print xrow
            #print yrow
            #print xseries
            #print yseries
            #data = data.append(xseries, ignore_index=True) 
            #data.ix[firstpoint:secondpoint, 'xcoord'] = xseries
            #data['xcoord'][firstpoint:lastpoint] = xseries
        xseries = pd.Series(xmain, name='xcoord')
        yseries = pd.Series(ymain, name='ycoord')
        data = pd.concat([data,xseries,yseries], axis=1) 
        return data
     
def print_full(x, filename):
    pd.set_option('display.max_rows', len(x))
    #print(x)
    pd.reset_option('display.max_rows')
    x.to_csv(filename,sep='\t')

if __name__=="__main__":
    hello = ReadAmr("./ens_1/fort.q0012")
    #print hello.AMR_level
    #print hello.x_low
    yoyo = hello.amrdataframe()
    #print yoyo.keys()
    #print yoyo["eta"]
    #print yoyo["amrlevel"]
    #print print_full(yoyo)
    #print print_full(yoyo['xcoord'])

    # In the most refined mesh, for a category, check 
    #print yoyo["eta"][(yoyo.xcoord == 72.5) & (yoyo.ycoord == -17.5)] 
    mama = yoyo["eta"][yoyo.amrlevel == 1.0]
    mama2 = mama.as_matrix()
    print mama2



