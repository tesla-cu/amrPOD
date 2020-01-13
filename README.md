# Efficient Algorithm to Compute POD on AMR Grids
Code is used to evaluate computational efficiency of a new algorithm to compute POD on AMR grids.

Background, discussion, goal, briefly how to use

writes data out in data folder



## How to use


## Organization
This repository is organized in the following manner:
  * `tests/` contains codes to test the algorithms used in this repository. Specifically, we test:
    * our synthetically grid generation technique in `grid_generation/`
    * our python code to Matlab  in `python_vs_matlab/`
    * POD using our reshaping technique and computation of POD against standard matrix operations in `reshaping/`

## Citation
CITATION. COPY amrex.

## Recreate data
To recreate the synthetic data in PAPER, simply do:
- `python write_synthetic_data.py` in `miscellaneous/`
- `./MasterRun.sh` in `drivers/`
- `python MasterPlot.py` in `plotting/` (note, comment out lines using genuine data)

The genuine data used for this paper can be found at URL. Download this data, then STUFF
- need to get `POD.batch`
- compile
- etc.

## License
LICENSE INFORMATION

### How to check if we are computing POD correctly
Questioning whether we did POD correctly with all this strange reshaping and convoluted operation skipping? Don't worry, we did too, so that is why we created some functions to check whether we actually did it correctly! To see for yourself, run `driver_check.py` with updated information (i.e., your directories where the AMR data is, the parameters of the simulation, etc.). Then run the corresponding Matlab script, again with updated parameters, where we do the standard operations. At the end of execution, Matlab will output the maximum error between the two computations, and the difference should be small.
