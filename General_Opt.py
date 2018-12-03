
# This is used to serve as a library such that some functions do not need to be rewritten over and ove again
import sys, os
from random import randint
import math, ipdb

# This function saves some general functions to be used for the optimization

def Constraint_Bounds(y_type):
    # This function is used generate the bounds for the constraint equations
    lb = []
    ub = []
    High_Bd_Val = float('Inf')
    for i in range(0, len(y_type)):
        if(y_type[i]>0):
            lb.append(0)
            ub.append(High_Bd_Val)
        else:
            lb.append(0)
            ub.append(0)
    return lb, ub
