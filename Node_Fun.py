import sys, os
import numpy as np
import math, ipdb
import ipdb
from compiler.ast import flatten
from collections import Counter
import pickle
from OwnLib import *

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
    -------------------------------------------
    """
    state_i = List_Append_fn(config_i, velocity_i)

    Robot_State_Update(world.robot(0), state_i)

    TreeNode = dict()
    TreeNode["Index"] = len(all_treenode)
    TreeNode["Config"] = config_i
    TreeNode["Velocity"] = velocity_i
    TreeNode["State"] = List_Append_fn(config_i, velocity_i)
    TreeNode["KE"] = Kinetic_Energy_fn(world.robot(0), state_i)
    TreeNode["Parent"] = -1
    TreeNode["Children"] = []
    TreeNode["Contact_PosNVel"] = Contact_Link_PosNVel(world.robot(0), contact_link_dictionary, -1)
    TreeNode["Contact_Status"] = contact_Status_Dictionary_i

    all_treenode.append(TreeNode)

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
