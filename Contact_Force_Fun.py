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

# This function is used to calculate the contact force related function

mu = 0.35

def Contact_Force_Constraint(contact_force, contact_link_list, contact_link_status, terr_model, y_val, y_type):
    # This function is used to take care of the contact force constraint given the contact force, contact link list and contact link status
    # The order of contact force is arranged according to the contact link list
    contact_force_index = 0
    for i in range(0, len(contact_link_list)):
        # This is the main force categories
        contact_link_i = contact_link_list[i]
        contact_link_i_status = contact_link_status[contact_link_i]
        for j in range(0, len(contact_link_i_status)):
            contact_link_i_point_j_status = contact_link_i_status[j]
            if contact_link_i_point_j_status == 1:
                # This means that the constraint force is active so should be feasible
                contact_force_of_link_i_at_point_j = Contact_Force_Element_from_Index(contact_force, contact_force_index)
                # Here this force is a 3 by 1 list
                _, _, Contact_Link_i_Point_j_Normal = Robot_Link_2_All_Terr_Dist(contact_force_of_link_i_at_point_j, terr_model)

                # Positivity constraint
                contact_force_of_link_i_at_point_j_normal = Dot_Product(contact_force_of_link_i_at_point_j, Contact_Link_i_Point_j_Normal)
                y_val.append(contact_force_of_link_i_at_point_j_normal)
                y_type.append(1)

                # Friction cone constraint
                contact_force_of_link_i_at_point_j_SoS = Dot_Product(contact_force_of_link_i_at_point_j, contact_force_of_link_i_at_point_j)
                contact_force_of_link_i_at_point_j_tang_sq = contact_force_of_link_i_at_point_j_SoS - contact_force_of_link_i_at_point_j_normal * contact_force_of_link_i_at_point_j_normal
                y_val.append(mu * mu * contact_force_of_link_i_at_point_j_normal * contact_force_of_link_i_at_point_j_normal - contact_force_of_link_i_at_point_j_tang_sq)
                y_type.append(1)
            else:
                # This means that the constraint force is notactive so should be zero
                contact_force_of_link_i_at_point_j = Contact_Force_Element_from_Index(contact_force, contact_force_index)
                List_Obj_Update(contact_force_of_link_i_at_point_j, 0, y_val, y_type)
            contact_force_index = contact_force_index + 1

def Contact_Force_Element_from_Index(contact_force, index_i):
    # This function gets the contact force element from certain given index_i
    start_index = 3 * index_i
    end_index = start_index + 3
    return contact_force[start_index:end_index]
