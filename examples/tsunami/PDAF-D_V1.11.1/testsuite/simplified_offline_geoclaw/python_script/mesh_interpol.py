import numpy as np

def interpol(in_2d_array):
    """ Interpolates from cell center to nodal points. The returned array has extra row and column. np.pad adds rows of zeros and columns before and after the (num1, num2). First tuple represents num1 rows added before the first row and num2 row added after the last row. Second tuple represents for the column.

    ARGUMENTS:
         in_2d_array : Numpy array that represents values at cell center

    RETURN:
         inter_value : Numpy array that represents values at nodes
"""

    inter_value = ((np.pad(in_2d_array,((0,1),(1,0)), mode='constant')) + 
(np.pad(in_2d_array,((0,1),(0,1)), mode='constant')) +  
(np.pad(in_2d_array,((1,0),(0,1)), mode='constant')) + 
(np.pad(in_2d_array,((1,0),(1,0)), mode='constant')))/4.0
    return inter_value

if __name__=='__main__':
    input_array = np.array([[1,2,3],[3,4,5],[2,4,4],[8,6,5]])
    print "input array is"
    print input_array
    output_array = interpol(input_array)
    print "output array is"
    print output_array
    
