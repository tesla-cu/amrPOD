### Efficient Algorithm to Compute POD on AMR Grids
Code is used to evaluate computational efficiency of a new algorithm to compute POD on AMR grids.


# How to check if we are computing POD correctly
Questioning whether we did POD correctly with all this strange reshaping and convoluted operation skipping? Don't worry, we did too, so that is why we created some functions to check whether we actually did it correctly! To see for yourself, run `driver_check.py` with updated information (i.e., your directories where the AMR data is, the parameters of the simulation, etc.). Then run the corresponding Matlab script, again with updated parameters, where we do the standard operations. At the end of execution, Matlab will output the maximum error between the two computations, and the difference should be small.
