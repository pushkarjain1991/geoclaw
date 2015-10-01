import numpy as np
import geoclaw_input_format as gcif

def pdaf_to_geoclaw(xv, yv, infile, outfile):
    """ Read pdaf format file and convert to geoclaw input format file"""
    z = np.loadtxt(infile)
    gcif.geoclaw_input_format(xv, yv, z, outfile)
    #print z
