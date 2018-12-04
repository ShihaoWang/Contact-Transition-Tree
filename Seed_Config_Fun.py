import sys, os
from random import randint
import numpy as np
import math, ipdb

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *

class Seed_Conf_Optimization_Prob(probFun):
    # This is the class type used for the initial condition validation optimization problem
    def __init__(self, world, state_init, contact_link_dictionary, contact_status_dictionary_init, terr_model):
        self.world = world;                                             self.state_init = state_init
        self.contact_link_dictionary = contact_link_dictionary;         self.contact_status_dictionary = contact_status_dictionary_init
        self.terr_model = terr_model
        # Here the derivative is set to be finite difference method
        self.grad = False
        nx = len(state_init)
        y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(world, state_init, contact_link_dictionary, contact_status_dictionary_init, terr_model, state_init)
        nc = len(y_type)
        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(self.world, self.state_init, self.contact_link_dictionary, self.contact_status_dictionary, self.terr_model, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

def Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, seed_state):
    # This function is useed to optimize a robot state such that it is on the constraint manifold

    # Update the robotstate according to the seed_state
    sim_robot = world.robot(0);                             Robot_State_Update(sim_robot, seed_state)
    # Constraint value and type
    y_val = [];                                             y_type = []
    # The objective function is the kinetic energy of the robot
    seed_conf_opt_obj_val = Kinetic_Energy_fn(sim_robot, seed_state)
    y_val.append(seed_conf_opt_obj_val);                    y_type.append(1)

    """
    Constraints
    1. Active constraints have to be maintained: Position and Velocity
    2. Inactive constraints have to be strictly away from the environment obstacles
    3*.If a contact will be added, its certain relative distance to terrain will be set to zero
    """

    Parant_Contact_PosNVel = treenode_parent["Contact_PosNVel"]
    Child_Contact_PosNVel = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)

    # Type of the contact at the child node
    contact_type, contact_link_cmp_dict = TreeNode_Status_CMP(treenode_parent, treenode_child)

    if contact_type == -1:
        # Contact retraction
        # The targeted contact status is the parent contact status
        

    else:
        # Contact addition or maintenance
        # The targeted contact status is the child contact status



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
    cfg.printLevel = 1;
    # cfg.printFile = "result.txt"
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
