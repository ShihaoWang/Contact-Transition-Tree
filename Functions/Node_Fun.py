import sys, os, copy
import numpy as np
import math
from compiler.ast import flatten
from collections import Counter
import pickle
from OwnLib import *
import datetime

# This file containts the functions of the node operation
def TreeNode_Dict_Init(world, config_i, velocity_i, contact_link_dictionary, contact_Status_Dictionary_i, all_treenode):
    # This function is used to initialize a Tree_Node dictionary

    """
    dictionary keys             type
    -------------------------------------------
    0.Node Index                int
    1.Config:                   list of double
    2.Velocity:                 list of double
    3.State:                    list of double
    4.Kinetic Energy(KE):       double
    5.Parent TreeNode           int
    6.Children TreeNodes        list of int
    7.Contact PosNVel info      list of dictionaries
    8.Contact Status info       list of dictionaries
    9.Contact Normal info       list of dictionaries
    10. Inertial Shaping(IS)    a list of optimized solution
    11. Parent to Child traj    a list of optimized solution
    -------------------------------------------
    """
    state_i = List_Append_fn(config_i, velocity_i)

    Robot_State_Update(world.robot(0), state_i)

    TreeNode = dict()
    TreeNode["Index"] = len(all_treenode)
    TreeNode["Config"] = config_i
    TreeNode["Velocity"] = velocity_i
    TreeNode["State"] = List_Append_fn(config_i, velocity_i)
    TreeNode["KE"] = world.robot(0).getKineticEnergy()
    TreeNode["Parent"] = -1
    TreeNode["Children"] = []
    TreeNode["Contact_PosNVel"] = Contact_Link_PosNVel(world.robot(0), contact_link_dictionary, -1)
    TreeNode["Contact_Status"] = contact_Status_Dictionary_i
    TreeNode["IS"] = []
    TreeNode["P2C"] = []

    all_treenode.append(TreeNode)

    return TreeNode

def TreeNode_Dict_Init_CS(contact_Status_Dictionary_i):
    # This function is used to initialze a treenode especially child node with only the information of the contact status
    TreeNode = dict()
    TreeNode["Index"] = -1
    TreeNode["Config"] = []
    TreeNode["Velocity"] = []
    TreeNode["State"] = []
    TreeNode["KE"] = 0
    TreeNode["Parent"] = -1
    TreeNode["Children"] = []
    TreeNode["Contact_PosNVel"] = {}
    TreeNode["Contact_Status"] = contact_Status_Dictionary_i
    TreeNode["IS"] = []
    TreeNode["P2C"] = []

    return TreeNode

def TreeNode_Parent_Update(treenode_child, treenode_parent):
    # This function is used to update the first treenode with the second treenode as its parent
    if treenode_child["Parent"] == -1:
        treenode_child["Parent"] = treenode_parent["Index"]
    else:
        raise RuntimeError("The given child node has already an unempty parent!\n")

def TreeNode_Child_Update(treenode_parent, treenode_child):
    # This function is used to update the first treenode with the second treenode as its child
    treenode_parent["Children"].append(treenode_child["Index"])

def Frontier_Add(frontier_nodes, treenode_i):
    # This function is used to add treenode to the frontier
    frontier_nodes.append(treenode_i)

def Frontier_Pop(frontier_nodes):
    # This function is used to pop out the treenode whose KE is the lowest
    if len(frontier_nodes)<1:
        raise RuntimeError("There is not node in Frontier!\n")
    else:
        frontier_nodes_len = len(frontier_nodes)
        KE_list = []
        for i in range(0, frontier_nodes_len):
            KE_list.append(frontier_nodes[i]["KE"])
        KE_min_index = KE_list.index(min(KE_list))
        treenode_i = frontier_nodes[KE_min_index]
        del frontier_nodes[KE_min_index]
        return treenode_i

def TreeNode_Status_CMP(treenode_parent, treenode_child):
    """
    This function is used to compare the two tree nodes
    There are three possible results:
                                        Contact Addition  (1)
                                        Contact Reduction (-1)
                                        Contact Unchanged (0)
    """
    treenode_parent_status = treenode_parent["Contact_Status"];
    treenode_child_status = treenode_child["Contact_Status"]

    treenode_parent_status_values = treenode_parent_status.values();
    treenode_child_status_values = treenode_child_status.values()

    treenode_parent_status_value_list = List_Flatten_fn(treenode_parent_status_values)
    treenode_child_status_value_list = List_Flatten_fn(treenode_child_status_values)

    # contact_status_cmp_res_list = List_Minus_fn(treenode_child_status_value_list, treenode_parent_status_value_list)
    # contact_status_cmp_res_max = max(contact_status_cmp_res_list)
    # contact_status_cmp_res_min = min(contact_status_cmp_res_list)

    contact_status = 0
    contact_status_cmp_dict = dict()
    contact_status_itc_dict = dict()

    contact_link_list = treenode_parent_status.keys()

    for i in range(0, len(treenode_parent_status_values)):
        treenode_parent_status_i = treenode_parent_status[contact_link_list[i]]
        treenode_child_status_i = treenode_child_status[contact_link_list[i]]
        contact_status_cmp_dict[contact_link_list[i]] = []
        contact_status_itc_dict[contact_link_list[i]] = []
        for j in range(0, len(treenode_parent_status_i)):
            treenode_parent_status_i_j = treenode_parent_status_i[j]
            treenode_child_status_i_j = treenode_child_status_i[j]
            treenode_parent2child_status_i_j = treenode_child_status_i_j - treenode_parent_status_i_j

            if treenode_parent_status_i_j == 1 and treenode_child_status_i_j == 1:
                contact_status_itc_dict[contact_link_list[i]].append(j)

            if treenode_parent2child_status_i_j == 0:
                continue
            else:
                if treenode_parent2child_status_i_j == 1:
                    contact_status = 1
                    contact_status_cmp_dict[contact_link_list[i]].append(j)
                else:
                    contact_status = -1
                    contact_status_cmp_dict[contact_link_list[i]].append(j)

    return contact_status, contact_status_cmp_dict, contact_status_itc_dict

def Node_Expansion_fn(treenode_parent, contact_link_list):
    """
        This function is used to conduct the node expansion for a given parent node
        The basic consideration is not to have the flying-in-air phase
    """
    parent_contact_status = treenode_parent["Contact_Status"].copy()       # Here parent_contact_status is a dictionary
    Children_nodes_contact_status_list = []
    for i in range(0, len(contact_link_list)):
        parent_contact_link_i_status = parent_contact_status[contact_link_list[i]]
        if len(parent_contact_link_i_status) == 1:
            # This means that this link has a point contact, the rule is to revert its value
            child_node_i_contact_status = copy.deepcopy(parent_contact_status)
            child_node_i_contact_status[contact_link_list[i]][0] = int(not(parent_contact_link_i_status[0]))
            Children_nodes_contact_status_list.append(child_node_i_contact_status)
        else:
            """
                Two possibilities:
                                    1. No contact:              ======>         In contact
                                    2. In contact:              ======>         Not in contact
            """
            vertices_number = len(parent_contact_link_i_status)
            child_node_i_contact_status_sum = sum(parent_contact_link_i_status)
            if child_node_i_contact_status_sum == 0:
                child_node_i_contact_status = copy.deepcopy(parent_contact_status)
                child_node_i_contact_status[contact_link_list[i]] = [1] * vertices_number
                Children_nodes_contact_status_list.append(child_node_i_contact_status)
            else:
                # ADD No contact
                child_node_i_contact_status = copy.deepcopy(parent_contact_status)
                child_node_i_contact_status[contact_link_list[i]] = [0] * vertices_number
                Children_nodes_contact_status_list.append(child_node_i_contact_status)
    return Children_nodes_contact_status_list

def Final_Node_2_Root_Node(Final_Node, All_Nodes):
    # This function is used to get the trajectories from the final node back to the root node
    Time_Duration_List = []
    Solution_List = []
    Time_Duration_List.append(Final_Node["IC"][0])
    Solution_List.append(Final_Node["IC"][1:])
    Current_Node = copy.deepcopy(Final_Node)
    while Current_Node["Parent"] != -1:
        # Attach to the front of the lists
        Time_Duration_List.insert(0, Current_Node["P2C"][0])
        Solution_List.insert(0,Current_Node["P2C"][1:])
        Current_Node = copy.deepcopy(All_Nodes[Current_Node["Parent"]])
    return Time_Duration_List, Solution_List

def CTT_Nodes_Write(Time_Duration_List, Solution_List, Grids_Number, Robot_Option):
    # This function is used to write the CTT nodes to file after get them
    """
        The optimized txt file should be of the following format:
        Grids8 (8 is a example of the grids number)
        Duration1
        Solution1
        Grids8
        Duration2
        Solution2
        .
        .
        .
    """
    # Get the proper name of the optimized solution
    currentDT = datetime.datetime.now()     # Current time

    file_name = str(currentDT.month) + "_" + str(currentDT.day) + "_" + str(currentDT.hour) + "_" + str(currentDT.minute) + "_Soln.txt"
    file_path_name = Robot_Option + file_name
    Txt_File = open(file_path_name, 'w')
    for i in range(0, len(Time_Duration_List)):
        # For each soln there is a default format
        soln_i_header = "Grids" + str(Grids_Number)
        print>>Txt_File, soln_i_header
        soln_i_time = Solution_List[i]
        print>>Txt_File, soln_i_time
        for j in Solution_List:
            print>>Txt_File, j
    Txt_File.close()

def CTT_Nodes_Read(Txt_File_Name, Robot_Option):
    """
        This function is used to read the CTT nodes solution
    """
    Time_Duration_List = []
    Solution_List = []
    file_name = Robot_Option + Txt_File_Name
    with open(file_name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        for Txt_File_Str_i in Txt_File_Str:
            if "Grids" in Txt_File_Str_i:
                # This indicates the header of a certain solution
                # Then the next few points would be the local coordinates of the contact extremities
                Grids_Number = Txt_File_Str_i.translate(None, 'Grids')     # This step is to get the link number out
                if len(Solution_List)>0:
                    Time_Duration_List_i = Raw_Soln_List_i[0]
                    Solution_List_i = Raw_Soln_List_i[1:]
                    Time_Duration_List.append(Time_Duration_List_i)
                    Solution_List.append(Solution_List_i)
                Raw_Soln_List_i = []
                continue
            Raw_Soln_List_i.append(float(Txt_File_Str_i))
    return Time_Duration_List, Solution_List, Grids_Number
