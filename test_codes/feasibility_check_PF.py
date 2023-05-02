# import libraries
import networkx
from pyomo.environ import *
import numpy as np
import scipy.io
import pdb
import math
import os


# get the current working directory
dir_file = os.getcwd()

# extract the scenarios and their probabilities from external file
scenario_file = dir_file + r"/scenario_files/failure_scenarios.mat"
probability_file = dir_file + r"/scenario_files/scenario_probabilities_49.csv"

# testcase to use for the optimization
test_case = '/14bus/'

# degree to radian conversion factor
deg_to_rad = math.pi / 180

# create a concrete pyomo model
model = ConcreteModel()

# load the bus, generation, and line data
model.bus = np.loadtxt(dir_file + test_case + r'bus_data.txt')
model.line = np.loadtxt(dir_file + test_case + r'line_data.txt')
model.gen = np.loadtxt(dir_file + test_case + r'gen_data.txt')
model.gen_cost = np.loadtxt(dir_file + test_case + r'gen_cost.txt')


# initialize the parameters
model.nEdges = len(model.line)  # total number of edges
model.nNodes = len(model.bus)  # total number of nodes
model.nGen = len(model.gen)  # total number of generators

# gen dataset also has specific MVA base; use that instead when updating this version
model.Pbase = 100  # MVA base 100 MVA to convert the system to pu system

# Pmax for bounds set
if max(model.line[:, 5]) == 0:
    Pmax_line = 100000/model.Pbase
else:
    Pmax_line = max(model.line[:, 5]) / model.Pbase

model.x_ij = range(0, model.nEdges)  # edges variable range
model.b_i = range(0, model.nNodes)  # buses variable range

###############################################################################################################
####################################### Variables #############################################################
###############################################################################################################

# declaring pyomo variables

# first we declare steady state variables i.e. power flow needs to be maintained while in the steady state
# since these should be the same for all scenarios, they are first stage variables
# they have ss at the end for representation

# although bounds are mentioned here to maintain the standard, they will be redefined as per gen bus
model.bus_gen_ss = Var(model.b_i, bounds=(0, 1000), within=Reals, initialize=0)  # bus generation variable
model.Pij_ss = Var(model.x_ij, bounds=(-Pmax_line, Pmax_line), within=Reals, initialize=0)  # active power flowing through each lines
model.theta_ss = Var(model.b_i, bounds=(-2*math.pi, 2*math.pi), within=Reals, initialize=0)  # angle of each bus
# model.load_shed = Var(model.b_i, bounds=(-math.inf, math.inf), within=Reals)  # real active power shed at each bus

###############################################################################################################
####################################### Constraints ###########################################################
###############################################################################################################

# pyomo constraints
model.c = ConstraintList()  # creates a list of constraints as placeholders

#################### bus power balance constraints ############################
# bus data col 3: active power demand, col 5: shunt conductance
for bus_num in range(model.nNodes):
    # identify the list of generators connected to each bus
    gens = np.where(model.gen[:, 0] == bus_num + 1)[0].tolist()
    to_bus_list = np.where(model.line[:, 1] == bus_num + 1)[0].tolist()
    from_bus_list = np.where(model.line[:, 0] == bus_num + 1)[0].tolist()

    model.c.add(sum(model.bus_gen_ss[gen_num] for gen_num in gens) +
                sum(model.Pij_ss[to_bus] for to_bus in to_bus_list) -
                sum(model.Pij_ss[from_bus] for from_bus in from_bus_list) ==
                model.bus[bus_num, 2] / model.Pbase + model.bus[bus_num, 4] / model.Pbase)

################## generator power limit constraint ###########################
# generator should generate power between its min and max active power limit
# col 9: PMAX and col 10: PMIN (Note: in Python number starts from 0)
for gen_num in range(model.nGen):
    model.c.add(model.bus_gen_ss[model.gen[gen_num, 0] - 1] <= model.gen[gen_num, 8] / model.Pbase)
    model.c.add(model.bus_gen_ss[model.gen[gen_num, 0] - 1] >= model.gen[gen_num, 9] / model.Pbase)

# make sure non-generating bus do not generate anything
for bus_num in range(model.nNodes):
    if not np.any(np.equal(model.gen[:, 0], bus_num+1)):
        model.c.add(model.bus_gen_ss[bus_num] == 0)

####################### active power flow constraint on each line ################
'''
Note: in Python number starts from 0

linedata:
col 4: reactance (X)
col 9: transformer tap ratio
col 10: transformer phase shift (in degrees)

busdata:
col 9: voltage angle (in degrees) -> this is a variable here so no need to use as parameter
'''

for line_num in range(model.nEdges):
    reciprocal_term = 1 / (model.line[line_num, 3] * model.line[line_num, 8])
    model.c.add(model.Pij_ss[line_num] == reciprocal_term * (model.theta_ss[model.line[line_num, 0] - 1] -
                                                             model.theta_ss[model.line[line_num, 1] - 1] -
                                                             (model.line[line_num, 9] * deg_to_rad)))

################### thermal limit (MVA_limits) ############################
# since the flow can be bi-directional, limits range from neg to positive value
# col 6: max MVA limit of the line (0 means unlimited capacity)
for line_num in range(model.nEdges):
    if model.line[line_num, 5] == 0:
        model.line[line_num, 5] = 10000
    model.c.add(model.Pij_ss[line_num] <= model.line[line_num, 5] / model.Pbase)
    model.c.add(model.Pij_ss[line_num] >= - model.line[line_num, 5] / model.Pbase)


################### angle difference between two buses on each line ################
# from bus and to bus reference is obtained via line
# col 12: min angle difference (degree), col 13: max angle difference (degree)
for angle_num in range(model.nEdges):
    model.c.add((model.theta_ss[model.line[angle_num, 0] - 1] - model.theta_ss[model.line[angle_num, 1] - 1]) <=
                model.line[angle_num, 12] * deg_to_rad)
    model.c.add((model.theta_ss[model.line[angle_num, 0] - 1] - model.theta_ss[model.line[angle_num, 1] - 1]) >=
                model.line[angle_num, 11] * deg_to_rad)

# the angle can be anywhere from -2pi to 2pi hence we need to maintain 0 angle at reference (slack) bus

# identifying slack bus
slack_bus = np.where(model.bus[:, 1] == 3)[0][0]

# ensure the angle at reference bus is 0
model.c.add(model.theta_ss[slack_bus] == 0)

# # ensure generation from G1 and G2 to match MATPOWER solution
# model.c.add(model.bus_gen_ss[0] == 2.21)
# model.c.add(model.bus_gen_ss[1] == 0.38)

############################# Overall Objective ###########################

def overall_objective(model):
    expr = 0
    expr = sum(model.gen_cost[i, 4] * (model.bus_gen_ss[model.gen[i, 0] - 1] * model.Pbase) ** 2 +
               model.gen_cost[i, 5] * (model.bus_gen_ss[model.gen[i, 0] - 1] * model.Pbase) +
               model.gen_cost[i, 6] for i in range(model.nGen))
    return expr

model.min_gen_cost = Objective(rule=overall_objective, sense=minimize)


####################### Solve the economic dispatch problem ################

# create lp file for debugging the model
model.write('check.lp', io_options={'symbolic_solver_labels': True})
solver = SolverFactory('gurobi')
results = solver.solve(model, tee=True)

# Observing the results

# for k in range(0, model.nGen):
    # print("G[%d] = %f" % (model.gen[k, 0], value(model.bus_gen_ss[model.gen[k, 0] - 1])))

for k in range(0, model.nNodes):
    print("G[%d] = %f" % (k+1, value(model.bus_gen_ss[k])))

sum_P = 0;
for k in range(0, model.nEdges):
    print("Pij[%d] = %f" % (k + 1, value(model.Pij_ss[k])))
    sum_P = sum_P + value(model.Pij_ss[k])

for k in range(0, model.nNodes):
    print("theta[%d] = %f" % (k + 1, value(model.theta_ss[k]) * 1 / deg_to_rad))
