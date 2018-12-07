import sys, os
from random import randint
import numpy as np
import math, ipdb

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *
from Seed_Guess_Fun import *

# This files contains all the files needed for the Contact Transition Tree (CTT) optimization

def Nodes_Optimization_fn(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option):
    """
    This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
    """
    # The outer operation runs over the enumeration of the duration guess
    Duration_min = 0.1
    Duration_max = 2.5
    Duration_Grids_number = 41
    Duration_array = np.linspace(Duration_min, Duration_max, Duration_Grids_number)
    Duration_list = Duration_array.tolist()

    Grids_Number = 8

    for Duration_i in Duration_list:
        # Inner operation
        Nodes_Optimization_Inner_Opt(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option,  Duration_i, Grids_Number)

def Nodes_Optimization_Inner_Opt(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids):
    Seed_Guess_List, DOF, Control_Force_Len = Seed_Guess_Gene(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids)

    Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, Seed_Guess_List)



def Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, control_force_len, seed_guess_list):
    # This is the core optimization function of this contact transition tree project
    sim_robot = world.robot(0);
    # Constraint value and type
    y_val = [];                                                 y_type = []

    Duration, State_List_Array, Control_List_Array, Contact_Force_List_Array = Seed_Guess_Unzip(seed_guess_list, DOF, control_force_len, grids)


    # Objective function initialization
    y_val.append(0);                                            y_type.append(1)

    """
        Constraint 1: Initial matching constraint condition
    """
    Parent_State = treenode_parent["State"];                    State_0 = State_List_Array[0]
    Initial_Match_Constraint = State_0 - Parent_State
    List_Obj_Update(Initial_Match_Constraint, 0, y_val, y_type)

    """
        Constraint 2: Dynamics constraints
                                            at collocation point (middle point)
    """

    for i in range(0, grids - 1):
        # Robot state dynamics at the front side
        State_Front = State_List_Array[i];
        Robot_State_Update(sim_robot, State_Front)
        D_q_Front, CG_q_qdot_Front, Jac_Full_Front, Jac_Full_Trans_Front = Dynamics_Matrices(sim_robot, contact_link_dictionary)
        Control_Front = Control_List_Array[i]
        Contact_Force_Front = Control_List_Array[i]
        Acc_Front = Dynamics_To_Acc(sim_robot, DOF, State_Front, Contact_Force_Front, Control_Front, contact_link_dictionary)

        # Robot state dynamics at the back side
        State_Back = State_List_Array[i + 1];
        Robot_State_Update(sim_robot, State_Back)
        D_q_Back, CG_q_qdot_Back, Jac_Full_Back, Jac_Full_Trans_Back = Dynamics_Matrices(sim_robot, contact_link_dictionary)
        Control_Back = Control_List_Array[i + 1]
        Contact_Force_Back = Control_List_Array[i + 1]
        Acc_Back = Dynamics_To_Acc(sim_robot, DOF, State_Back, Contact_Force_Back, Control_Back, contact_link_dictionary)








    # Update the robotstate according to the seed_state

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
            Child_Contact_PosNVel_Link_i_Pos_j_Dist = Robot_Link_2_All_Terr_Dist(Child_Contact_PosNVel_Link_i_Pos_j, terr_model)
            if Child_Contact_PosNVel_Link_i_Pos_j_Dist<0:
                y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Dist)
                y_type.append(1)

    # 3. Contact Addition constraint
    if contact_type == 1:
        for i in range(0, len(contact_link_list)):
            Child_Contact_PosNVel_Link_i = Child_Contact_PosNVel[i]
            Child_Contact_PosNVel_Link_i_Pos = Child_Contact_PosNVel_Link_i["Pos"]
            for j in contact_link_cmp_dict[contact_link_list[i]]:
                Child_Contact_PosNVel_Link_i_Pos_j = Child_Contact_PosNVel_Link_i_Pos[j]
                Child_Contact_PosNVel_Link_i_Pos_j_Dist = Robot_Link_2_All_Terr_Dist(Child_Contact_PosNVel_Link_i_Pos_j, terr_model)
                y_val.append(Child_Contact_PosNVel_Link_i_Pos_j_Dist)
                y_type.append(0)
    return y_val, y_type
