import sys, os
from random import randint
import numpy as np
import math

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *

class Seed_Conf_Optimization_Prob(probFun):
    # This is the class type used for the initial condition validation optimization problem
    def __init__(self, world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, config_ref):
        self.world = world;
        self.treenode_parent = treenode_parent;                         self.treenode_child = treenode_child;
        self.contact_link_dictionary = contact_link_dictionary;         self.terr_model = terr_model
        self.config_ref = config_ref
        # Here the derivative is set to be finite difference method
        self.grad = False

        nx = len(config_ref)
        y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, config_ref, config_ref)
        nc = len(y_type)
        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(self.world, self.treenode_parent, self.treenode_child, self.contact_link_dictionary, self.terr_model, self.config_ref, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

def Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, config_ref, config_opt):
    # This function is useed to optimize a robot state such that it is on the constraint manifold
    # Update the robotstate according to the seed_state
    sim_robot = world.robot(0);
    sim_robot.setConfig(config_opt)
    # Constraint value and type
    y_val = [];                                             y_type = []
    # The objective function is the kinetic energy of the robot or the difference between the parent state and the seed_state
    obj_val = 0
    for i in range(0,len(config_ref)):
        obj_val = obj_val + (config_opt[i] - config_ref[i]) * (config_opt[i] - config_ref[i])
    y_val.append(obj_val);          y_type.append(1)

    contact_link_list = contact_link_dictionary.keys()
    """
    Constraints
    1. Previous Active contacts have to be maintained.
    2. All contacts have to not have any penetration with the environment obstacles
    3*.If a contact will be added, its certain relative distance to terrain will be set to zero
    """

    Parent_Contact_PosNVel = treenode_parent["Contact_PosNVel"]                                     # list of dictionaries according to contact_link_list
    Child_Contact_PosNVel = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)

    # Type of the contact at the child node
    contact_type, contact_link_cmp_dict, contact_link_itc_dict = TreeNode_Status_CMP(treenode_parent, treenode_child)

    # 1. Active Maintenance Constraint
    for i in range(0, len(contact_link_itc_dict)):
        contact_link_i_number = contact_link_list[i]
        Parent_Link_i_Dict = Parent_Contact_PosNVel[i]
        Child_Link_i_Dict = Child_Contact_PosNVel[i]
        for j in contact_link_itc_dict[contact_link_i_number]:
            # This is the order of Parent_Contact_PosNVel/Child_Contact_PosNVel list
            # If not empty, there exists active pos/vel
            Parent_Link_i_Pos_j = Parent_Link_i_Dict["Pos"][j]
            Child_Link_i_Pos_j = Child_Link_i_Dict["Pos"][j]
            Parent_Child_Link_i_Pos_j_Constraint = List_Minus_fn(Parent_Link_i_Pos_j, Child_Link_i_Pos_j)
            List_Obj_Update(Parent_Child_Link_i_Pos_j_Constraint, 0, y_val, y_type)

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
                y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Dist)
                y_type.append(0)
                y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Edge)
                y_type.append(1)
    return y_val, y_type

def Seed_Guess_Gene_Robotstate(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_config):
    """
        This function is used to project the reference state into a goal state that lies on the constrained manifold to initialize the upper level optimization
        Since that zero velocity will always satisfy the constraint, then the default velocity of the goal state is set to be zero.
    """
    # The inputs to this function is the WorldModel, initial contact mode, and the initial state
    # The output is the feasible robotstate that lies on the constrained manifold

    qmin, qmax = world.robot(0).getJointLimits();
    # Bounds on the optimized variables: configuration
    xlb = qmin;                                 xub = qmax
    # Optimization problem setup
    Seed_Conf_Opt = Seed_Conf_Optimization_Prob(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_config)
    Seed_Conf_Opt.xlb = xlb;                   Seed_Conf_Opt.xub = xub
    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Seed_Conf_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_config, robot_ref_config)
    lb, ub = Constraint_Bounds(y_type)
    Seed_Conf_Opt.lb = lb;                     Seed_Conf_Opt.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 0;
    # cfg.printFile = "result.txt"
    cfg.majorIterLimit = 150
    slv = solver(Seed_Conf_Opt, cfg)
    # rst = slv.solveRand()
    rst = slv.solveGuess(np.array(robot_ref_config))
    #  Here a basic solution separation will be conducted to get the configuraiton and velocity
    Goal_Config_List = []
    for i in range(0, rst.sol.size):
        Goal_Config_List.append(rst.sol[i])
    return Goal_Config_List
