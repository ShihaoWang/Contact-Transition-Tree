#!/usr/bin/python

import osqp
import sys, os
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin
import numpy as np
import ipdb, string
import scipy as sp
import scipy.sparse as sparse
from klampt.model.trajectory import Trajectory
import time
import math

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *
from Seed_Config_Fun import *
from Seed_Guess_Fun import *
from Dynamics_Fun import *
from CTT_Fun import *

# Some system parameters
Terr_Model = []                                      # The generalized terrain model structure
Contact_Link_Dictionary = dict()                     # Information of the robot contact link extremities
All_TreeNodes = []                                   # All of the tree nodes
Frontier_Nodes = []                                  # Frontier_Nodes
Robot_Option =""
def main():

    # This funciton is used for the multi-contact humanoid push recovery
    # The default robot to be loaded is the HRP2 robot in this same folder
    Robot_Option = "./User_File/HRP2_Robot/"
    # Robot_Option = "./User_File/JQ_Robot/"

    print "This funciton is used for the 3-D Humanoid Multi-Contact Fall Mitigation"
    if len(sys.argv)<=1:
        print "USAGE: The default robot to be loaded is the JQ robot"
        exit()
    world = WorldModel()                    # WorldModel is a pre-defined class
    input_files = sys.argv[1:];             # sys.argv will automatically capture the input files' names
    for fn in input_files:
        result = world.readFile(fn)         # Here result is a boolean variable indicating the result of this loading operation
        if not result:
            raise RuntimeError("Unable to load model " + fn)


    # The definition of the environment into an efficient way
    Terr_Model = Terr_Model_Cal(world)

    # The following function can be used in two ways: the first way is to load the Config_Init.config file while the second way is to load two
    # DOF, Config_Init, Velocity_Init = State_Loader_fn("Init_Config.config", Robot_Option)
    # DOF, Config_Init, Velocity_Init = State_Loader_fn("Init_Config.txt", "Init_Velocity.txt", Robot_Option)
    DOF, Config_Init, Velocity_Init = State_Loader_fn("Opt_Init_Config.txt", "Opt_Init_Velocity.txt", Robot_Option)

    # State_Writer_fn(Config_Init, "Init_Config_from_txt.config", Robot_Option)
    # State_Writer_fn(Config_Init, Velocity_Init, "Inn_Config.txt", "Inn_Velo.txt",Robot_Option)

    Contact_Link_Dictionary = Contact_Link_Reader("Contact_Link.txt", Robot_Option)
    Contact_Status_Dictionary_Init = Contact_Status_Reader("Init_Contact.txt", Robot_Option)
    # According to the initial condition of the robot contact status, a basic optimization may have to be contacted to enforce the initial constraints.

    # Now it is the validation of the feasibility of the given initial condition
    State_Init = List_Append_fn(Config_Init, Velocity_Init)

    Config_Init, Velocity_Init = Robot_Init_Opt_fn(world, State_Init, Contact_Link_Dictionary, Contact_Status_Dictionary_Init, Terr_Model)
    # State_Writer_fn(Config_Init, Velocity_Init, "Opt_Init_Config.txt", "Opt_Init_Velocity.txt",Robot_Option)

    # Now it is the contact transition tree optimization
    # Root node initialization
    Root_Node = TreeNode_Dict_Init(world, Config_Init, Velocity_Init, Contact_Link_Dictionary, Contact_Status_Dictionary_Init, All_TreeNodes)
    Frontier_Add(Frontier_Nodes, Root_Node)

    Impact_Mapping_fn(world, State_Init, Root_Node, Contact_Link_Dictionary)


    while len(Frontier_Nodes)>0:
        treenode_parent = Frontier_Pop(Frontier_Nodes)
        """"
        * For the current node, first is the Node_Self_Opt to optimize a motion while maintain the current mode
        * if this does not work, then expand the current node into the adjacent nodes then do the Nodes_Connectivity_Opt
        """
        # Opt_Soln, Opt_Flag, Final_State = Nodes_Optimization_fn(world, treenode_parent, Current_Node, Contact_Link_Dictionary, Terr_Model, Robot_Option)
        Opt_Flag = False
        if Opt_Flag == True:
            # Here this means that the robot has already achieved a stabilized state with inertia shaping
            treenode_parent["IS"] = Opt_Soln
        else:
            # This means the node expansion has to be conducted to enable the stabilization with the modification of robot contacts
            Children_Nodes_Contact_Status_List = Node_Expansion_fn(treenode_parent, Contact_Link_Dictionary.keys())
            """
                Then the job is to test the connectivity of the children nodes from the current node
            """
            ipdb.set_trace()
            for i in range(0, len(Children_Nodes_Contact_Status_List)):
                child_node_i_contact_status = Children_Nodes_Contact_Status_List[i]
                # However, here the child_node_i has not been initialized
                treenode_child = TreeNode_Dict_Init_CS(child_node_i_contact_status)
                Opt_Soln, Opt_Flag, Final_State = Nodes_Optimization_fn(world, treenode_parent, treenode_child, Contact_Link_Dictionary, Terr_Model, Robot_Option)
                if Opt_Flag == True:
                    # This means that this child node is connectable
                    treenode_child["P2C"] = Opt_Soln
                    contact_type, contact_link_cmp_dict, contact_link_itc_dict = TreeNode_Status_CMP(treenode_parent, treenode_child)
                    if contact_type == 1:
                        # This means that there is impact mapping
                        pass

                else:
                    pass

if __name__ == "__main__":
    main()
