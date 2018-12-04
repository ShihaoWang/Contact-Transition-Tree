import sys, os
from random import randint
import numpy as np
import math, ipdb

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *

# This file contains all the functions needed for the initialization of the seed for nonlinear optimization


def Seed_Guess_Gene(treenode_parent, treenode_child, duration):
    """
    This function will generate the optimization variables needed for the further optimization
    """

    Seed_Guess_Gene_Robotstate(treenode_parent, treenode_child);




def Seed_Guess_Gene_Robotstate(treenode_parent, treenode_child):
    """
    This function is used to generate a configuration to initialize the upper level optimization
    """
