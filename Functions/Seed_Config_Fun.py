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
    def __init__(self, world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_state):
        self.world = world;
        self.treenode_parent = treenode_parent;                         self.treenode_child = treenode_child;
        self.contact_link_dictionary = contact_link_dictionary;         self.terr_model = terr_model
        self.robot_ref_state = robot_ref_state
        # Here the derivative is set to be finite difference method
        self.grad = False

        nx = len(treenode_parent["State"])
        y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_state, robot_ref_state)
        nc = len(y_type)
        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(self.world, self.treenode_parent, self.treenode_child, self.contact_link_dictionary, self.terr_model, self.robot_ref_state, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

def Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, ref_state, seed_state):
    # This function is useed to optimize a robot state such that it is on the constraint manifold
    # Update the robotstate according to the seed_state
    sim_robot = world.robot(0);                             Robot_State_Update(sim_robot, seed_state)
    # Constraint value and type
    y_val = [];                                             y_type = []
    # The objective function is the kinetic energy of the robot or the difference between the parent state and the seed_state
    # Seed_Conf_Opt_Obj_Val = Kinetic_Energy_fn(sim_robot, seed_state)

    State_Diff_List = List_Minus_fn(seed_state, ref_state)           # This step is actually not right
    Seed_Conf_Opt_Obj_Val = Dot_Product(State_Diff_List, State_Diff_List)

    y_val.append(Seed_Conf_Opt_Obj_Val);                    y_type.append(1)

    contact_link_list = contact_link_dictionary.keys()

    """
    Constraints
    1. Active contacts have to be maintained: Position and Velocity
    2. All contacts have to not have any penetration with the environment obstacles
    3*.If a contact will be added, its certain relative distance to terrain will be set to zero
    """

    Parent_Contact_PosNVel = treenode_parent["Contact_PosNVel"]                                     # list of dictionaries according to contact_link_list
    Child_Contact_PosNVel = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)

    # Type of the contact at the child node
    contact_type, contact_link_cmp_dict, contact_link_itc_dict = TreeNode_Status_CMP(treenode_parent, treenode_child)
    # ipdb.set_trace()

    # 1. Active Maintenance Constraint
    for i in range(0, len(contact_link_itc_dict)):
        contact_link_i_number = contact_link_list[i]
        Parent_Link_i_Dict = Parent_Contact_PosNVel[i]
        Child_Link_i_Dict = Child_Contact_PosNVel[i]
        for j in contact_link_itc_dict[contact_link_i_number]:
            # This is the order of Parent_Contact_PosNVel/Child_Contact_PosNVel list
            # If not empty, there exists active pos/vel
            Parent_Link_i_Pos_j = Parent_Link_i_Dict["Pos"][j];             Parent_Link_i_Vel_j = Parent_Link_i_Dict["Vel"][j]
            Child_Link_i_Pos_j = Child_Link_i_Dict["Pos"][j];               Child_Link_i_Vel_j = Child_Link_i_Dict["Vel"][j]

            Parent_Child_Link_i_Pos_j_Constraint = List_Minus_fn(Parent_Link_i_Pos_j, Child_Link_i_Pos_j)
            Parent_Child_Link_i_Vel_j_Constraint = List_Minus_fn(Parent_Link_i_Vel_j, Child_Link_i_Vel_j)

            List_Obj_Update(Parent_Child_Link_i_Pos_j_Constraint, 0, y_val, y_type)
            List_Obj_Update(Parent_Child_Link_i_Vel_j_Constraint, 0, y_val, y_type)

    # 2. Environment penetration constraint
    for i in range(0, len(contact_link_list)):
        Child_Contact_PosNVel_Link_i = Child_Contact_PosNVel[i]
        Child_Contact_PosNVel_Link_i_Pos = Child_Contact_PosNVel_Link_i["Pos"]
        for j in range(0, len(Child_Contact_PosNVel_Link_i_Pos)):
            Child_Contact_PosNVel_Link_i_Pos_j = Child_Contact_PosNVel_Link_i_Pos[j]
            Child_Contact_PosNVel_Link_i_Pos_j_Dist, _, _ = Robot_Link_2_All_Terr_Dist(Child_Contact_PosNVel_Link_i_Pos_j, terr_model)
            y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Dist)
            y_type.append(1)

    # 3. Contact Addition constraint
    if contact_type == 1:
        for i in range(0, len(contact_link_list)):
            Child_Contact_PosNVel_Link_i = Child_Contact_PosNVel[i]
            Child_Contact_PosNVel_Link_i_Pos = Child_Contact_PosNVel_Link_i["Pos"]
            for j in contact_link_cmp_dict[contact_link_list[i]]:
                Child_Contact_PosNVel_Link_i_Pos_j = Child_Contact_PosNVel_Link_i_Pos[j]
                Child_Contact_PosNVel_Link_i_Pos_j_Dist, Child_Contact_PosNVel_Link_i_Pos_j_Edge, _ = Robot_Link_2_All_Terr_Dist(Child_Contact_PosNVel_Link_i_Pos_j, terr_model)
                # ipdb.set_trace()
                y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Dist)
                y_type.append(0)
                y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Edge)
                y_type.append(1)
    return y_val, y_type

def Seed_Guess_Gene_Robotstate(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_state):
    """
    This function is used to project the reference state into a goal state that lies on the constrained manifold to initialize the upper level optimization
    """
    # The inputs to this function is the WorldModel, initial contact mode, and the initial state
    # The output is the feasible robotstate that lies on the constrained manifold

    qmin, qmax = world.robot(0).getJointLimits();    dqmax_val = world.robot(0).getVelocityLimits()
    # Bounds on the optimized variables
    xlb = qmin;                                 xub = qmax
    for i in range(0,len(dqmax_val)):
        xlb.append(-dqmax_val[i]);              xub.append(dqmax_val[i])

    # Optimization problem setup
    Seed_Conf_Opt = Seed_Conf_Optimization_Prob(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_state)
    Seed_Conf_Opt.xlb = xlb;                   Seed_Conf_Opt.xub = xub
    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_state, robot_ref_state)
    lb, ub = Constraint_Bounds(y_type)
    Seed_Conf_Opt.lb = lb;                     Seed_Conf_Opt.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 1;
    # cfg.printFile = "result.txt"
    cfg.majorIterLimit = 300

    slv = solver(Seed_Conf_Opt, cfg)
    # rst = slv.solveRand()
    rst = slv.solveGuess(np.array(robot_ref_state))

    #  Here a basic solution separation will be conducted to get the configuraiton and velocity
    Seed_State_List = []
    for i in range(0, rst.sol.size):
        Seed_State_List.append(rst.sol[i])

    Seed_Config = Seed_State_List[0:rst.sol.size/2]
    Seed_Velocity = Seed_State_List[rst.sol.size/2:]
    return Seed_Config, Seed_Velocity
