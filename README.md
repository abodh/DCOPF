# DC Optimal Power Flow
This repository contains a generalized DCOPF formulation for standard IEEE test cases using concrete modeling approach in Pyomo. The current version is tested for IEEE-14 bus, IEEE-39 bus, IEEE-118 bus, and IEEE-300 bus system with the MATPOWER solution. I will test bigger test cases considering my time availability this Summer :) 

You can use this DCOPF formulation in any of your pyomo-based models for restoration, planning, unit committment, and so forth. Please cite this repository if used for academic or commericial research. 

# DCOPF Objective
The current objective minimizes the cost of generation for all the existing generators in the system (economic dispatch). 

# DCOPF Constraints
1. Power balance in each of the nodes
2. Power flow in each of the lines (including transformer turn ratio and phase shifts)
3. Generator's power output limit
4. Line thermal limit (power carrying capability)
5. Limit in the angle differences between the buses
6. Reference bus has 0 voltage angle

# Notes
1. The standard test cases are mat files of the standard MATPOWER cases
2. This code assumes the same data columns as there are in the MATPOWER test cases
3. Although Gurobi is used for the current version, you can use any open source solvers like glpk and cbc as long as you do not introduce integers, and other complexities. The model has also been tested successfully on CPLEX.    

# Dependencies
Please make sure that you install and import the following packages
1. `numpy`
2. `scipy`
3. `pyomo`
4. `gurobi`/`cplex`/`glpk` or any other optimization solver
