# Efficient Algorithm to Compute POD on AMR Grids
This repo hosts code to compute proper orthogonal decomposition (POD) on a set of grid that utilized adaptive mesh refinement and evaluate the efficiency of the algorithm compared to the standard algorithm. 

To evaluate the efficiecy of the algorithm, we provide a variety of parameters sweeps in FOLDER which use the python source code to count the number of operations and outputs the data into a new folder `data/` that sits one level up from these code directories. This is particular useful to see how the AMR algorithm works (all documented code is in files with `_CPU.py`) and to develop new algorithms in a user friendly environment.

To actually compute POD, we recommend using the Fortran version found in `fortran/parallel/` which computes POD in an efficient, scalable manner which also can bypass memory issues since storing all data in RAM can be quite costly. The code here can easily compute POD using the standard or AMR algorithm, so it is not restricted to only AMR data.

## Organization
This repository is organized in the following manner:
  * `drivers/` STUFF
  * `fortran/` contains codes to both:
    * CPU times comparing the standard and AMR algorithms in `CPU/`
    * fully compute POD in an efficient parallel manner in `parallel/`
  * `miscellaneous/` contains codes for miscellaneous tasks
  * `plotting/` contains all codes to generate figures
  * `source/` contains all python code that actually computes each operation of POD
  * `tests/` contains codes to test the algorithms used in this repository. Specifically, we test:
    * our synthetically grid generation technique in `grid_generation/`
    * our python code to Matlab  in `python_vs_matlab/`
    * POD using our reshaping technique and computation of POD against standard matrix operations in `reshaping/`

## How to use
If you are interested in testing a new algorithm, make those changes in the source code then use the code in `reshaping/` to test if the result is faster.

If you are interested in spanning new parameter spaces, simply copy one of the DRIVERS similar to your needs and edit to your liking, then run with `python code.py` in DRIVERS folder. Weights of the operations can be changed in `source/Compute_POD.py`.

If you would like to compute CPU time, use code in `fortran/CPU/`. Compile lines are at the top of `POD.f90`.

If you would like to compute POD as fast and efficiently as possible, use code in `fortran/parallel/`. Compile lines are at the top of `POD.f90`. This code is a hydridization MPI/OMP code for HPC use and alleviates memory issues by partial loading.

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
