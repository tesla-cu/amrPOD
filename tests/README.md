_To run_ : 
- `cd tests/grid_generation/`
- set parameters in `grid_generation.py`
- `python grid_generation.py`
- check if `ls` is the same as `l_comp`, and if `lc` is the same as `lc_comp` (Note, that with small `nt`, `lc` may be larger than set)



_To run_ : 
- `cd tests/python_vs_matlab/`
- set parameters and AMR directory in `POD_pvm.py` and `POD_pvm.m`
- `python POD_pvm.py`
- in Matlab, run `POD_pvm.m`
- check that the error is sufficiently small


_To run_ : 
- `cd tests/reshaping/`
- set parameters in `POD_rshp.py`
- `python POD_rshp.py`

_To test column-major order reshaping_ :
- in Matlab, run `column_order_shape.m`
- check that the repeated cells in the column vector are contiguous and the final grid is the same as the original grid 