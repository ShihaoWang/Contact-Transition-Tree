import sys, os
from random import randint
sys.path.insert(0, '/home/shihao/trajOptLib')
from trajOptLib.io import getOnOffArgs
from trajOptLib import trajOptCollocProblem
from trajOptLib.snoptWrapper import directSolve
from trajOptLib.libsnopt import snoptConfig, probFun, solver
import functools
import numpy as np
import math
from OwnLib import *
from Terrain_Fun import *

# This function is used to setup the initial optimization problem
# A modification will be made to have two optimizations: first for configuration, second for velocity.

"""
        Optimization for Robot Configuration
"""
def Initial_Config_Constraint_Eqn(world, config_ref, contact_link_dictionary, contact_status_dictionary, terr_model, x):
    # This function is used to generate the constraint equations for optimized configuration.
    sim_robot = world.robot(0)
    # Here x is only the configuration variable list
    sim_robot.setConfig(x)
    y_val = [];                     y_type = []
    """
    The objective is the deviation from the initial given configuration to the optimized configuration
    """
    obj_val = 0
    for i in range(0,len(x)):
        obj_val = obj_val + (x[i] - config_ref[i]) * (x[i] - config_ref[i])
    y_val.append(obj_val);          y_type.append(1)
    """
    # Constraints
    Constraints on the position and velocity of contact extremities
    Constraints on distance to the edge of the object
    """
    contact_link_list = contact_status_dictionary.keys()
    for i in contact_link_list:
        # For each link, there are some contact points
        link_i_contact_points = contact_link_dictionary[i]
        link_i_contact_status = contact_status_dictionary[i]
        Contact_Link_i_PosNVel_List = Contact_Link_i_PosNVel(sim_robot, i, link_i_contact_points)         # a dictionary with keys (Pos, Vel)
        for j in range(0, len(link_i_contact_points)):
            Contact_Link_i_PosNVel_Point_j_status = link_i_contact_status[j]
            Contact_Link_i_PosNVel_Point_j_Pos = Contact_Link_i_PosNVel_List["Pos"][j]
            Contact_Link_i_PosNVel_Point_j_Face_Dist, Contact_Link_i_PosNVel_Point_j_Edge_Dist, _ = Robot_Link_2_All_Terr_Dist(Contact_Link_i_PosNVel_Point_j_Pos, terr_model)

            if Contact_Link_i_PosNVel_Point_j_status == 1:
                # In this case, the constraint has to be active
                # Vertical distance
                y_val.append(Contact_Link_i_PosNVel_Point_j_Face_Dist);                 y_type.append(0)
                # Distance to the edge of contact surface
                y_val.append(Contact_Link_i_PosNVel_Point_j_Edge_Dist);                 y_type.append(1)
            else:
                # In this case, the constraint has to be inactive
                # However, the robot has to remain not having any penetration
                y_val.append(Contact_Link_i_PosNVel_Point_j_Face_Dist);                 y_type.append(1)
    return y_val, y_type

class Robot_Init_Config_Opt_Prob(probFun):
    # This is the class type used for the initial condition validation optimization problem
    def __init__(self, world, config_ref, contact_link_dictionary, contact_status_dictionary_init, terr_model):
        self.world = world;                                             self.config_ref = config_ref
        self.contact_link_dictionary = contact_link_dictionary;         self.contact_status_dictionary = contact_status_dictionary_init
        self.terr_model = terr_model

        # Here the derivative is set to be finite difference method
        self.grad = False

        nx = len(config_ref)
        y_val, y_type = Initial_Config_Constraint_Eqn(world, config_ref, contact_link_dictionary, contact_status_dictionary_init, terr_model, config_ref)
        nc = len(y_type)
        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Initial_Config_Constraint_Eqn(self.world, self.config_ref, self.contact_link_dictionary, self.contact_status_dictionary, self.terr_model, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

def Robot_Init_Config_Opt_fn(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model):
    # The inputs to this function is the WorldModel, initial contact mode, and the initial state
    # The output is the feasible robot configuration
    qmin, qmax = world.robot(0).getJointLimits();
    # Bounds on the optimized variables
    xlb = qmin;                                 xub = qmax

    # Optimization problem setup
    Robot_Init_Opt = Robot_Init_Config_Opt_Prob(world, state_init[0:len(qmin)], contact_link_dictionary, contact_status_dictionary, terr_model)
    Robot_Init_Opt.xlb = xlb;                   Robot_Init_Opt.xub = xub

    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Initial_Config_Constraint_Eqn(world, state_init[0:len(qmin)], contact_link_dictionary, contact_status_dictionary, terr_model, state_init[0:len(qmin)])
    lb, ub = Constraint_Bounds(y_type)
    Robot_Init_Opt.lb = lb;                     Robot_Init_Opt.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 0;
    # cfg.printFile = "result.txt"
    cfg.majorIterLimit = 100
    slv = solver(Robot_Init_Opt, cfg)
    # rst = slv.solveRand()
    rst = slv.solveGuess(np.array(state_init[0:len(qmin)]))
    #  Here a basic solution separation will be conducted to get the configuraiton and velocity
    Initial_Config_List = []
    for i in range(0, rst.sol.size):
        Initial_Config_List.append(rst.sol[i])
    return Initial_Config_List

"""
        Optimization for Robot Velocity
"""

class Robot_Init_Velocity_Opt_Prob(probFun):
    # This is the class type used for the initial condition validation optimization problem
    def __init__(self, world, config_ref, contact_link_dictionary, contact_status_dictionary_init, terr_model):
        self.world = world;                                             self.config_ref = config_ref
        self.contact_link_dictionary = contact_link_dictionary;         self.contact_status_dictionary = contact_status_dictionary_init
        self.terr_model = terr_model

        # Here the derivative is set to be finite difference method
        self.grad = False

        nx = len(config_ref)
        y_val, y_type = Initial_Velocity_Constraint_Eqn(world, config_ref, contact_link_dictionary, contact_status_dictionary_init, terr_model, config_ref)
        nc = len(y_type)
        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Initial_Velocity_Constraint_Eqn(self.world, self.config_ref, self.contact_link_dictionary, self.contact_status_dictionary, self.terr_model, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

def Initial_Velocity_Constraint_Eqn(world, config_ref, contact_link_dictionary, contact_status_dictionary, terr_model, x):
    # This function is used to generate the constraint equations for the optimization
    # Due to the specific usage of this function, it can only be used inside the structure where self is defined
    sim_robot = world.robot(0)
    # Now the x should be only the velocity list
    sim_robot.setConfig(config_ref)
    sim_robot.setVelocity(x)

    y_val = [];                     y_type = []

    """
    The objective for the velocity optimization is not very important.
    """
    obj_val = 0.0
    y_val.append(obj_val);          y_type.append(1)

    """
    # Constraints
    Constraints on the velocity of contact extremities
    Constraints on kinetic energy
    """
    contact_link_list = contact_status_dictionary.keys()

    for i in contact_link_list:
        # For each link, there are some contact points
        link_i_contact_points = contact_link_dictionary[i]
        link_i_contact_status = contact_status_dictionary[i]
        Contact_Link_i_PosNVel_List = Contact_Link_i_PosNVel(sim_robot, i, link_i_contact_points)         # a dictionary with keys (Pos, Vel)
        for j in range(0, len(link_i_contact_status)):
            Contact_Link_i_PosNVel_Point_j_status = link_i_contact_status[j]
            Contact_Link_i_PosNVel_Point_j_Vel = Contact_Link_i_PosNVel_List["Vel"][j]
            if Contact_Link_i_PosNVel_Point_j_status == 1:
                # In this case, the constraint has to be active
                # Velocity
                y_val.append(Contact_Link_i_PosNVel_Point_j_Vel[0]);                   y_type.append(0)
                y_val.append(Contact_Link_i_PosNVel_Point_j_Vel[1]);                   y_type.append(0)
                y_val.append(Contact_Link_i_PosNVel_Point_j_Vel[2]);                   y_type.append(0)

    # Constraint to set the initial kinetic energy
    KE_ref = 50
    KE_val = sim_robot.getKineticEnergy()
    y_val.append(KE_val - KE_ref);                                y_type.append(0)
    return y_val, y_type

def Robot_Init_Velocity_Opt_fn(world, config_ref, contact_link_dictionary, contact_status_dictionary, terr_model):
    # The inputs to this function is the WorldModel, initial contact mode, and the initial state
    # The output is the feasible robot configuration
    dqmax_val = world.robot(0).getVelocityLimits()
    # Bounds on the optimized variables
    xlb = []
    xub = []
    for i in range(0,len(dqmax_val)):
        xlb.append(-dqmax_val[i]); xub.append(dqmax_val[i])
    # Optimization problem setup
    Robot_Init_Opt = Robot_Init_Velocity_Opt_Prob(world, config_ref, contact_link_dictionary, contact_status_dictionary, terr_model)
    Robot_Init_Opt.xlb = xlb;                   Robot_Init_Opt.xub = xub

    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Initial_Velocity_Constraint_Eqn(world, config_ref, contact_link_dictionary, contact_status_dictionary, terr_model, config_ref)
    lb, ub = Constraint_Bounds(y_type)
    Robot_Init_Opt.lb = lb;                     Robot_Init_Opt.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 0;
    # cfg.printFile = "result.txt"
    cfg.majorIterLimit = 100
    slv = solver(Robot_Init_Opt, cfg)
    rst = slv.solveRand()
    # rst = slv.solveGuess(np.array(state_init[0:len(config_ref)]))
    #  Here a basic solution separation will be conducted to get the configuraiton and velocity
    Initial_Velocity_List = []
    for i in range(0, rst.sol.size):
        Initial_Velocity_List.append(rst.sol[i])
    return Initial_Velocity_List

def Robot_Init_Opt_fn(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model):

    # The inputs to this function is the WorldModel, initial contact mode, and the initial state
    # The output is the feasible robotstate
    Config_Init = Robot_Init_Config_Opt_fn(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model)
    Velocity_Init = Robot_Init_Velocity_Opt_fn(world, Config_Init, contact_link_dictionary, contact_status_dictionary, terr_model)
    return Config_Init, Velocity_Init
