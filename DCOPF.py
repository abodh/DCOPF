# import libraries
from pyomo.environ import *
import numpy as np
import scipy.io as readmat
import pdb
import math
import os

# get the current working directory
dir_file = os.getcwd()

# testcase to use for the optimization
test_case = '300'

# degree to radian conversion factor
deg_to_rad = math.pi / 180

# create a concrete pyomo model
model = ConcreteModel()

# load matpower test case in mat format
matpower_mat_file = readmat.loadmat(dir_file + '/power_system_test_cases/case' + test_case + '.mat',
                                    struct_as_record=False,
                                    squeeze_me=True)
# ensure that all the saved mat file are saved under workspace var name 'matpower_testcase'
test_case = matpower_mat_file['matpower_testcase']

# load the bus, generation, and line data
model.bus = test_case.bus
model.line = test_case.branch
model.gen = test_case.gen
model.gen_cost = test_case.gencost

# initialize the parameters
model.nEdges = len(model.line)  # total number of edges
model.nNodes = len(model.bus)  # total number of nodes
model.nGen = len(model.gen)  # total number of generators

# gen dataset also has specific MVA base; use that instead when updating this version
model.Pbase = 100  # MVA base 100 MVA to convert the system to pu system

# let us assume that infinite capacity is equal to 100 GW
Inf_transfer_Pmax = 10e6
Pmax_line = 10e6/model.Pbase

Gmax = max(model.gen[:, 8]) / model.Pbase
max_load = max(model.bus[:, 2]) / model.Pbase

# variable ranges
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
model.bus_gen_ss = Var(model.b_i, bounds=(0, Gmax), within=Reals, initialize=0)  # bus generation variable
model.Pij_ss = Var(model.x_ij, bounds=(-Pmax_line, Pmax_line), within=Reals,
                   initialize=0)  # active power flowing through each lines
model.theta_ss = Var(model.b_i, bounds=(-2 * math.pi, 2 * math.pi), within=Reals, initialize=0)  # angle of each bus

###############################################################################################################
####################################### Constraints ###########################################################
###############################################################################################################

# pyomo constraints
# creates a list of constraints as placeholders
#################### bus power balance constraints ############################
model.power_balance = ConstraintList()

# bus data col 3: active power demand, col 5: shunt conductance
for bus_num_idx in range(model.nNodes):
    bus_num = model.bus[bus_num_idx, 0]

    # identify the list of generators connected to each bus
    gens = np.where(model.gen[:, 0] == bus_num)[0].tolist()
    to_bus_list = np.where(model.line[:, 1] == bus_num)[0].tolist()
    from_bus_list = np.where(model.line[:, 0] == bus_num)[0].tolist()

    model.power_balance.add(sum(model.bus_gen_ss[np.where(model.bus[:, 0] == model.gen[gen_num, 0])[0][0]]
                                for gen_num in gens) +
                            sum(model.Pij_ss[to_bus] for to_bus in to_bus_list) -
                            sum(model.Pij_ss[from_bus] for from_bus in from_bus_list) ==
                            model.bus[bus_num_idx, 2] / model.Pbase + model.bus[bus_num_idx, 4] / model.Pbase)

################## generator power limit constraint ###########################
model.gen_limit = ConstraintList()
# generator should generate power between its min and max active power limit
# col 9: PMAX and col 10: PMIN (Note: in Python number starts from 0)
for gen_num in range(model.nGen):
    model.gen_limit.add(model.bus_gen_ss[np.where(model.bus[:, 0] == model.gen[gen_num, 0])[0][0]] <=
                        model.gen[gen_num, 8] / model.Pbase)
    model.gen_limit.add(model.bus_gen_ss[np.where(model.bus[:, 0] == model.gen[gen_num, 0])[0][0]] >=
                        model.gen[gen_num, 9] / model.Pbase)

# make sure non-generating bus do not generate anything
for bus_num_idx in range(model.nNodes):
    bus_num = model.bus[bus_num_idx, 0]

    if not np.any(np.equal(model.gen[:, 0], bus_num)):
        model.gen_limit.add(model.bus_gen_ss[bus_num_idx] == 0)

####################### active power flow constraint on each line ################
model.power_flow = ConstraintList()
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

    # MATPOWER keeps 0 for transmission lines without transformer
    # here we need to ensure tap ratio for transmission line is 1
    if model.line[line_num, 8] == 0:
        model.line[line_num, 8] = 1

    reciprocal_term = 1 / (model.line[line_num, 3] * model.line[line_num, 8])
    model.power_flow.add(model.Pij_ss[line_num] ==
                         reciprocal_term * (model.theta_ss[np.where(model.bus[:, 0] == model.line[line_num, 0])[0][0]] -
                                            model.theta_ss[np.where(model.bus[:, 0] == model.line[line_num, 1])[0][0]] -
                                            (model.line[line_num, 9] * deg_to_rad)))

################### thermal limit (MVA_limits) ############################
# since the flow can be bi-directional, limits range from neg to positive value
# col 6: max MVA limit of the line (0 means unlimited capacity)
# this constraint tightens the -inf, inf bound set during variable initialization
for line_num in range(model.nEdges):
    if model.line[line_num, 5] == 0:
        model.line[line_num, 5] = Pmax_line
    model.power_flow.add(model.Pij_ss[line_num] <= model.line[line_num, 5] / model.Pbase)
    model.power_flow.add(model.Pij_ss[line_num] >= - model.line[line_num, 5] / model.Pbase)

################### angle difference between two buses on each line ################
model.angle_limit = ConstraintList()
# from bus and to bus reference is obtained via line
# col 12: min angle difference (degree), col 13: max angle difference (degree)
for angle_num in range(model.nEdges):
    model.angle_limit.add((model.theta_ss[np.where(model.bus[:, 0] == model.line[angle_num, 0])[0][0]] -
                           model.theta_ss[np.where(model.bus[:, 0] == model.line[angle_num, 1])[0][0]])
                          <= model.line[angle_num, 12] * deg_to_rad)
    model.angle_limit.add((model.theta_ss[np.where(model.bus[:, 0] == model.line[angle_num, 0])[0][0]] -
                           model.theta_ss[np.where(model.bus[:, 0] == model.line[angle_num, 1])[0][0]])
                          >= model.line[angle_num, 11] * deg_to_rad)

# the angle can be anywhere from -2pi to 2pi hence we need to maintain 0 angle at reference (slack) bus

# identifying slack bus
slack_bus = np.where(model.bus[:, 1] == 3)[0][0]

# ensure the angle at reference bus is 0
model.angle_limit.add(model.theta_ss[slack_bus] == 0)


############################# Overall Objective ###########################

def overall_objective(model):
    expr = sum(model.gen_cost[gen_num, 4] * (model.bus_gen_ss[np.where(model.bus[:, 0] == model.gen[gen_num, 0])[0][0]]
                                             * model.Pbase) ** 2 +
               model.gen_cost[gen_num, 5] * (model.bus_gen_ss[np.where(model.bus[:, 0] == model.gen[gen_num, 0])[0][0]]
                                             * model.Pbase) +
               model.gen_cost[gen_num, 6] for gen_num in range(model.nGen))
    return expr


model.min_gen_cost = Objective(rule=overall_objective, sense=minimize)

####################### Solve the economic dispatch problem ################

# create lp file for debugging the model
model.write('check.lp', io_options={'symbolic_solver_labels': True})
solver = SolverFactory('gurobi')
results = solver.solve(model, tee=True)

# Observing the results
for k in range(0, model.nGen):
    print("G[%d] = %f" % (model.gen[k, 0], value(model.bus_gen_ss[np.where(model.bus[:, 0] ==
                                                                           model.gen[k, 0])[0][0]])))

sum_P = 0;
for k in range(0, model.nEdges):
    print("Pij[%d] = %f" % (k + 1, value(model.Pij_ss[k])))
    sum_P = sum_P + value(model.Pij_ss[k])

for k in range(0, model.nNodes):
    print("theta[%d] = %f" % (k + 1, value(model.theta_ss[k]) * 1 / deg_to_rad))
