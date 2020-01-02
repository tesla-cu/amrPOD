# Efficient Algorithm to Compute POD on AMR Grids
Code is used to evaluate computational efficiency of a new algorithm to compute POD on AMR grids.

Background, discussion, goal, briefly how to use

writes data out in data folder


## Organization
This repository is organized in the following manner:
  * tests/
  
    Contains codes to test the algorithms used in this repository. Specifically, we test:
    * grid_generation/
    
      This tests our synthetically grid generation technique.
      
      _To run_ : 
      - `cd tests/grid_generation/`
      - set parameters in `grid_generation.py`
      - `python grid_generation.py`
      - check if `ls` is the same as `l_comp`, and if `lc` is the same as `lc_comp` (Note, that with small `nt`, `lc` may be larger than set)
    * python_vs_matlab/
    
      This compares our python code to Matlab.
      
      _To run_ : 
      - `cd tests/python_vs_matlab/`
      - set parameters and AMR directory in `POD_pvm.py` and `POD_pvm.m`
      - `python POD_pvm.py`
      - in Matlab, run `POD_pvm.m`
      - check that the error is sufficiently small
      
    * reshaping/
    
      This tests that the data that is a result from the POD computation is the same between the unshaped data and the reshaped.
      
      _To run_ : 
      - `cd tests/reshaping/`
      - set parameters in `POD_rshp.py`
      - `python POD_rshp.py`
      
      _To test column-major order reshaping_ :
      - in Matlab, run `column_order_shape.m`
      - check that the repeated cells in the column vector are contiguous and the final grid is the same as the original grid

## How to use

## Citation
CITATION. COPY amrex.

## License
LICENSE INFORMATION

### How to check if we are computing POD correctly
Questioning whether we did POD correctly with all this strange reshaping and convoluted operation skipping? Don't worry, we did too, so that is why we created some functions to check whether we actually did it correctly! To see for yourself, run `driver_check.py` with updated information (i.e., your directories where the AMR data is, the parameters of the simulation, etc.). Then run the corresponding Matlab script, again with updated parameters, where we do the standard operations. At the end of execution, Matlab will output the maximum error between the two computations, and the difference should be small.
