BP_IF.m is the function to predict the influence function (the process parameters, material properties, etc. are all defined in the function, which can be also adjusted with support of experiments)

To predict the entire removal map, convolution based calculation can be employed, but one can also run TopoGen.m, which allows predict the removal map on workpiece even if the influence function is varying from position to position. 
For parallel calculation, run TopoGen_ParallelCal.m
