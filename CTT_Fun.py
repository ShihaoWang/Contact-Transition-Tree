import sys, os
from random import randint
import numpy as np
import math, ipdb

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *

# This files contains all the files needed for the Contact Transition Tree (CTT) optimization

def Nodes_Optimization_fn(treenode_parent, treenode_child):
    """
    This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
    """
    # The outer operation runs over the enumeration of the duration guess
    Duration_min = 0.1
    Duration_max = 2.5
    Duration_Grids_number = 41
    Duration_array = np.linspace(Duration_min, Duration_max, Duration_Grids_number)
    Duration_list = Duration_array.tolist()

    for Duration_i in Duration_list:
        # Inner operation
        Nodes_Optimization_Inner_Opt(treenode_parent, treenode_child, Duration_i)

def Nodes_Optimization_Inner_Opt(treenode_parent, treenode_child, duration):
