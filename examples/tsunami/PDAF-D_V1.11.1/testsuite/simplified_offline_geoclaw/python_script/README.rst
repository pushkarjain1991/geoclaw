README
======

Overview
++++++++

The python script **clean_connector.py** does the following::
    #. Constructs initial water surface elevation (analogous to **hump.xyz** of bowl_radial)
    #. Constructs an ensemble of initial water surface elevations
    #. Constructs observation data at **dx** * **dy** inetervals on mesh for dtobs list of times
    #. Constructs individual directories with name ``ens_``, for running individual ensembles
    #. Performs forecast step
    #. Performs assimilation step


PDAF input format
+++++++++++++++++

--- --- --- --- (nx times) --- --- --- ---
--- --- --- --- (nx times) --- --- --- ---|
--- --- --- --- (nx times) --- --- --- ---| (ny times)
--- --- --- --- (nx times) --- --- --- ---|
--- --- --- --- (nx times) --- --- --- ---


GeoClaw input format
++++++++++++++++++++

# x   y   z
 --- --- --- 
 --- --- ---|
 --- --- ---| (nx * ny times)
 --- --- ---|
 --- --- ---|
 --- --- ---


.. note:: 
   #. Define the model parameters and the DA parameters and run the script 
   #. Make sure that the num_ens in Python script matches the num_ens in python script
   #. Compile PDAF and copy the executable in your current folder
   #. Currently, it is assumed that we have observation data at the nodal points 

Plotting
++++++++
python $CLAW/visclaw/src/python/visclaw/plotclaw.py ./ ./ setplot.py

