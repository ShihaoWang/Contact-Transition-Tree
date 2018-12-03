import sys, os
from random import randint
sys.path.insert(0, '/home/shihao/trajOptLib')
from trajOptLib.io import getOnOffArgs
from trajOptLib import trajOptCollocProblem
from trajOptLib.snoptWrapper import directSolve
from trajOptLib.libsnopt import snoptConfig, probFun, solver
import functools
import numpy as np
import math, ipdb
from OwnLib import *
from Terrain_Fun import *
from General_Opt import *

# This function is used to setup the initial optimization problem

class Robot_Init_Opt_Prob(probFun):
    # This is the class type used for the initial condition validation optimization problem
    def __init__(self, world, state_init, contact_link_dictionary, contact_status_dictionary_init, terr_model):

        self.world = world;                                             self.state_init = state_init
        self.contact_link_dictionary = contact_link_dictionary;         self.contact_status_dictionary = contact_status_dictionary_init
        self.terr_model = terr_model

        # Here the derivative is set to be finite difference method
        self.grad = False

        nx = len(state_init)
        y_val, y_type = Initial_Constraint_Eqn(world, state_init, contact_link_dictionary, contact_status_dictionary_init, terr_model, state_init)
        nc = len(y_type)

        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Initial_Constraint_Eqn(self.world, self.state_init, self.contact_link_dictionary, self.contact_status_dictionary, self.terr_model, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass
        # y[0] = x[0] ** 2 + x[1] ** 2
        # y[1] = x[0] + x[1]
        # G[:2] = 2 * x
        # G[2:] = 1.0
        # if rec:
        #     row[:] = [0, 0, 1, 1]
        #     col[:] = [0, 1, 0, 1]

def Initial_Constraint_Eqn(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model, x):
    # This function is used to generate the constraint equations for the optimization
    # Due to the specific usage of this function, it can only be used inside the structure where self is defined
    sim_robot = world.robot(0)
    Robot_State_Update(sim_robot, x)

    y_val = [];                     y_type = []

    """
    The objective is the deviation from the initial given configuration to the optimized configuration
    """
    obj_val = 0
    for i in range(0,len(x)/2):
        obj_val = obj_val + (x[i] - state_init[i]) * (x[i] - state_init[i])

    y_val.append(obj_val);          y_type.append(1)

    # ipdb.set_trace()

    """
    # Constraints
    Constraints on the position and velocity of contact extremities
    Constraints on distance to the edge of the object
    """
    robot_contact_link_list = contact_status_dictionary.keys()
    # ipdb.set_trace()

    for i in robot_contact_link_list:
        # For each link, there are some contact points
        robot_link_i_contact_points = contact_link_dictionary[i]
        robot_link_i_contact_status = contact_status_dictionary[i]
        Contact_Link_i_PosNVel_List = Contact_Link_i_PosNVel(sim_robot, i, robot_link_i_contact_points)         # a dictionary with keys (Pos, Vel)
        for j in range(0, len(robot_link_i_contact_points)):
            Contact_Link_i_PosNVel_Point_j_status = robot_link_i_contact_status[j]
            Contact_Link_i_PosNVel_Point_j_Pos = Contact_Link_i_PosNVel_List["Pos"][j]
            Contact_Link_i_PosNVel_Point_j_Vel = Contact_Link_i_PosNVel_List["Vel"][j]
            Contact_Link_i_PosNVel_Point_j_Vel_SoS = Dot_Product(Contact_Link_i_PosNVel_Point_j_Vel, Contact_Link_i_PosNVel_Point_j_Vel)        # Sum of Squares

            Contact_Link_i_PosNVel_Point_j_Face_Dist, Contact_Link_i_PosNVel_Point_j_Edge_Dist = Robot_Link_2_All_Terr_Dist(Contact_Link_i_PosNVel_Point_j_Pos, terr_model)

            if Contact_Link_i_PosNVel_Point_j_status == 1:
                # In this case, the constraint has to be active
                # Position
                y_val.append(Contact_Link_i_PosNVel_Point_j_Face_Dist);                 y_type.append(0)
                # Velocity
                y_val.append(Contact_Link_i_PosNVel_Point_j_Vel_SoS);                   y_type.append(0)
                # Distance to the edge of contact surface
                y_val.append(Contact_Link_i_PosNVel_Point_j_Edge_Dist);                 y_type.append(1)
            else:
                # In this case, the constraint has to be inactive
                # However, the robot has to remain not having any penetration
                y_val.append(Contact_Link_i_PosNVel_Point_j_Face_Dist);                 y_type.append(1)

    # Constraint to set the initial kinetic energy
    KE_ref = 30
    KE_val = Kinetic_Energy_fn(sim_robot, x)

    y_val.append((KE_val - KE_ref) * (KE_val - KE_ref));                                y_type.append(0)
    return y_val, y_type

def Robot_Init_Opt_fn(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model):

    # The inputs to this function is the WorldModel, initial contact mode, and the initial state
    # The output is the feasible robotstate

    qmin, qmax = world.robot(0).getJointLimits();    dqmax_val = world.robot(0).getVelocityLimits()
    # Bounds on the optimized variables
    xlb = qmin;                                 xub = qmax
    for i in range(0,len(dqmax_val)):
        xlb.append(-dqmax_val[i]);              xub.append(dqmax_val[i])

    # Optimization problem setup
    Robot_Init_Opt = Robot_Init_Opt_Prob(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model)
    Robot_Init_Opt.xlb = xlb;                   Robot_Init_Opt.xub = xub

    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Initial_Constraint_Eqn(world, state_init, contact_link_dictionary, contact_status_dictionary, terr_model, state_init)
    lb, ub = Constraint_Bounds(y_type)
    Robot_Init_Opt.lb = lb;                     Robot_Init_Opt.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 0;                         cfg.printFile = "result.txt"
    cfg.majorIterLimit = 300
    # ipdb.set_trace()

    slv = solver(Robot_Init_Opt, cfg)
    # rst = slv.solveRand()
    rst = slv.solveGuess(np.array(state_init))

    #  Here a basic solution separation will be conducted to get the configuraiton and velocity
    Initial_State_List = []
    for i in range(0, rst.sol.size):
        Initial_State_List.append(rst.sol[i])

    Init_Config = Initial_State_List[0:rst.sol.size/2]
    Init_Velocity = Initial_State_List[rst.sol.size/2:]
    return Init_Config, Init_Velocity
