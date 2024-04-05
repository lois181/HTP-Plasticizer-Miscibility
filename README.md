# HTP-Plasticizer-Miscibility
A high throughput procedure (HTP) to relate the miscibility in a polyisoprene melt of simplified plasticizers (PLs) of varying topologies to their geometrical and thermodynamic properties. The key idea to this procedure is to test various PL 'descriptors' that are cheap to calculate compared to their phase behaviour in the polymer melt and, as such, prepare the data neccessary for learning methods such as the Decision Tree or Logistic Regression which can create desgin paths to miscible PLs that require minimal computational resources. 

# Software and Installation

We use the Molecular Dynamics (MD) software GROMACS 2016.4 for all simulations. This version is required for the use of the 'group' cut-off, which is removed in later versions. For instructions on how to install it, see here: `https://manual.gromacs.org/2016.4/download.html`

To manage the HTP, we use the workflow manager `signac` and `signac-flow` which are avaliable as python packages. The following instructions assume a working knowledge of the software. The required Signac packages can be installed using `pip`

    pip install Signac
    pip install Signac-flow

and to clone the repository:

    git clone https://github.com/lois181/HTP-Plasticizer-Miscibility.git


# Running Procedure 

## Preparing Plasticizer Input 

The `signac` statepoints allow the user to change geometric features of PLs of interest, along with their concentration and NVT production run times. A `.txt` file (`parameters.txt`) is provided to the script which initialises the `signac` procedure (`init.py`). This utilises `plasticizer_class.py` to provide each job with the correct GROMACS input files to run each type of simulation. `parameters.txt` should be formatted according to the following example (column titles provided for user):

    backbone length    side chain length   side chain frequency    flexibility    simulation time (ps)  simulation type    PL beads in PL/PL box    PL concentration in PI/PL simulations (phr)
           10                  5                    5                   f               10000                C5_C5                  6500                                 5

The `side chain frequency` refers to the number of beads separating side chains, which are always evenly spaced. Note that the side chain placement always begins on the second bead of the backbone and a side chain is never placed on the last bead of the backbone. Setting the side chain length to zero will make this state variable redundant but it should still be included in the file. The simulation type can be `C_C5` (box of PLs in polyisoprene), `C5_C5` (box of PLs), `C_C` (box of polyisoprene) or `C5` (box of one PL); depending on the system of interest. All types of simulations can be performed in the same `signac` project but note that certain state variables are redundant depending on simulation type (for example the PL concentration when selecting the `C5` simulation type). 

The project can then be initialised: 

    python init.py run 

which creates a `my_project` directory containing the `workspace` directory. Each job directory contains the GROMACS input files neccessary for the `simulation type` specified in `parameters.txt`. A list of the job directory names and associated statepoints is also printed. 

## Preparing cluster submission

The procedure is designed to be compatible with the 2016.4 MPI version of GROMACS (signle precision) with a SLURM scheduler. We have modified the default cluster submission template provided by `signac`, please refer to the `signac` documentation for guidance on how to submit on your own scheduler: (`https://docs.signac.io/en/latest/cluster_submission.html`). The number of cores provided to each `signac` operation can be modified using the `directives` argument in `project.py`. 

To perform the simulations, first move into the initialised project root directory

        cd my_project

## Signac-flow procedure for PL/PI simulations (C_C5) and miscibility analysis 

The flow project performs the PI/PL box creation, minimisation, NPT equilibration, NVT production run and miscibility analysis based on the PL/PL cumulative radial distribution number in one operation. In the case of a failed simulation, resubmission of the operation will commence from the last successful step. This behaviour can be modified by changing the operation post conditions. 

To submit the operation:

    python project_c_c5.py submit -o c_sim

To check its progress:

    python project_c_c5.py status 

Upon success, the result is provided in binary format in `signac_job_document.json` of each relevant job directory (0- miscible, 1- immiscible).

## PL/PL Simulations (C5_C5) and configurational entropy analysis 

To submit the operation:

    python project_c5_c5.py submit -o c5_c5_sim 

Upon success, the result for the average PL configurational entropy in 600 ns blocks is printed as a dictionary to `signac_job_document.json`. 

## PL/vacuum simulations and square radius of gyration, acylindricity and configurational entropy of a single PL analysis 

To submit the operation:

    python project_c5.py submit -o c5_sim

Upon success, the square radii of gyration, acylindricities and configurational entropies are printed as a dictionary to `signac_job_document.json`

## Collecting results 

Results from successful simulations can currently be collected by integrating the project with `pandas DataFrames`. For example:

    import pandas as pd

    df = my_project.to_dataframe()
    

# Known Issues

1. The allocation of resources is constant through each operation which can lead to sub-optimal resource management for long simulations. Users can currently perform the GROMACS simulations and analyses in separate `signac-flow` projects to circumvent this. 
2. The ease of collation and presentation of results is not yet optimised, users can see the `signac dashboard` package for ease of searching by state variable

