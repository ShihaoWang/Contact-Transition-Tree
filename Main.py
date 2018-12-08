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
    global Robot_Option;
    Robot_Option = "./User_File/HRP2_Robot/"
    # Robot_Option = "./User_File/JQ_Robot/"
    # ipdb.set_trace()

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
    global Terr_Model;                  Terr_Model = Terr_Model_Cal(world)

    # The following function can be used in two ways: the first way is to load the Config_Init.config file while the second way is to load two
    # DOF, Config_Init, Velocity_Init = State_Loader_fn("Init_Config.config", Robot_Option)
    # DOF, Config_Init, Velocity_Init = State_Loader_fn("Init_Config.txt", "Init_Velocity.txt", Robot_Option)
    DOF, Config_Init, Velocity_Init = State_Loader_fn("Opt_Init_Config.txt", "Opt_Init_Velocity.txt", Robot_Option)

    # State_Writer_fn(Config_Init, "Init_Config_from_txt.config", Robot_Option)
    # State_Writer_fn(Config_Init, Velocity_Init, "Inn_Config.txt", "Inn_Velo.txt",Robot_Option)

    global Contact_Link_Dictionary;
    Contact_Link_Dictionary = Contact_Link_Reader("Contact_Link.txt", Robot_Option)
    Contact_Status_Dictionary_Init = Contact_Status_Reader("Init_Contact.txt", Robot_Option)

    # According to the initial condition of the robot contact status, a basic optimization may have to be contacted to enforce the initial constraints.
    # ipdb.set_trace()
    # # # Now it is the validation of the feasibility of the given initial condition
    State_Init = List_Append_fn(Config_Init, Velocity_Init)

    Config_Init, Velocity_Init = Robot_Init_Opt_fn(world, State_Init, Contact_Link_Dictionary, Contact_Status_Dictionary_Init, Terr_Model)
    # State_Writer_fn(Config_Init, Velocity_Init, "Opt_Init_Config.txt", "Opt_Init_Velocity.txt",Robot_Option)

    # ipdb.set_trace()
    ## Now it is the contact transition tree optimization
    # Root node initialization
    Root_Node = TreeNode_Dict_Init(world, Config_Init, Velocity_Init, Contact_Link_Dictionary, Contact_Status_Dictionary_Init, All_TreeNodes)
    Frontier_Add(Frontier_Nodes, Root_Node)
    # Seed_Conf_Optimization_ObjNConstraint(world, Root_Node, Root_Node, Contact_Link_Dictionary, Terr_Model, State_Init)
    # Seed_Guess_Gene(world, Root_Node, Root_Node, Contact_Link_Dictionary, Terr_Model, Robot_Option, 1.0, 8)
    while len(Frontier_Nodes)>0:
        pass

        """"
        * For the current node, first is the Node_Self_Opt to optimize a motion while maintain the current mode
        * if this does not work, then expand the current node into the adjacent nodes then do the Nodes_Connectivity_Opt
        """
        Nodes_Optimization_fn(world, Root_Node, Root_Node, Contact_Link_Dictionary, Terr_Model, Robot_Option)


    # Given the pre-optimized robot state, we could directly load em in
    robot_viewer = MyGLViewer(world, Config_Init, Velocity_Init)
    ipdb.set_trace()

    robot_viewer.drawContacts = True
    robot_viewer.drawSensors = True
    vis.setPlugin(robot_viewer)
    vis.show()

    robot_viewer_time = time.time()

    robot_sim = robot_viewer.sim
    sim_robot_controller = robot_sim.controller(0)

    # Up to here the green shadow is no longer there anymore
    Robot_State_Update(world.robot(0), State_Init)

    while vis.shown():
        time.sleep(0.1)
        robot_sim.updateWorld()
        sim_robot_controller.setMilestone(Config_Init, Velocity_Init)

if __name__ == "__main__":
    main()
