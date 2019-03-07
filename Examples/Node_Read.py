# This file provides a tutorial on how to read a node file and show the visualization.
import sys, os
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin
import numpy as np
import string
import scipy as sp
import scipy.sparse as sparse
from klampt.model.trajectory import Trajectory
import time
import math

sys.path.insert(0, '../Functions')

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *
from Seed_Config_Fun import *
from Seed_Guess_Fun import *
from Dynamics_Fun import *
from CTT_Fun import *


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

    All_TreeNodes = []                                   # All of the tree nodes
    Frontier_Nodes = []                                  # Frontier_Nodes

    """
        Some system variables: Grids_Number
    """
    Grids_Number = 8

    # The definition of the environment into an efficient way
    Terr_Model = Terr_Model_Cal(world)

    # Here the contact transition tree has already been finished!
    # Write soln to file!
    # CTT_Nodes_Write(Time_Duration_List, Solution_List, Grids_Number, Robot_Option)

    # Read soln from file!
    Txt_File_Name =
    Time_Duration_List, Solution_List = CTT_Nodes_Read(Txt_File_Name, Robot_Option)

    # Total robot animation
    Total_Robot_Motion_Plot(world, DOF, Contact_Link_Dictionary, Terr_Model, Robot_Option, Grids_Number, Time_Duration_List, Solution_List)

if __name__ == "__main__":
    main()
