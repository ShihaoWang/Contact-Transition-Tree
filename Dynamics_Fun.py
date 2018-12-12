import sys, os
from random import randint
import numpy as np
import math, ipdb

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *

# This function is used to get the dynamics function
def Dynamics_Matrices(robot, contact_link_dictionary):
    # Here the robot object should have already been updated according to the new state
    D_q = robot.getMassMatrix()                             # List of list
    C_q_qdot = robot.getCoriolisForces()                    # List
    G_q = robot.getGravityForces([0,0,-9.8])                # List
    # ipdb.set_trace()
    CG_q_qdot = List_Sum_fn(C_q_qdot, G_q)
    Jac_Full = []
    for contact_link_index_i in contact_link_dictionary.keys():
        robot_link_i = robot.link(contact_link_index_i)
        contact_link_index_i_extremities = contact_link_dictionary[contact_link_index_i]
        for j in range(0, len(contact_link_index_i_extremities)):
            robot_link_i_point_j = contact_link_index_i_extremities[j]
            jac_link_i_point_j = robot_link_i.getPositionJacobian(robot_link_i_point_j)
            for jac_link_i_point_j_coord_z in jac_link_i_point_j:
                Jac_Full.append(jac_link_i_point_j_coord_z)
    Jac_Full_Trans = map(list, zip(*Jac_Full))
    return D_q, CG_q_qdot, Jac_Full, Jac_Full_Trans

def Dynamics_RHS_Matrix_fn(DOF, Jac_Full_Trans):
    # This function will generate the right hand side matrix needed to conduct the matrix inversion
    Dynamics_RHS_Matrix = []
    for i in range(0, DOF):
        Jac_Full_Trans_ith_row = Jac_Full_Trans[i]
        B_ith_row = [0] * (DOF - 6)
        if i>=6:
            B_ith_row[i - 6] = 1
        Dynamics_RHS_Matrix_i_row = Jac_Full_Trans_ith_row + B_ith_row
        Dynamics_RHS_Matrix.append(Dynamics_RHS_Matrix_i_row)
    return Dynamics_RHS_Matrix

def Dynamics_To_Acc(sim_robot, DOF, robot_state, contact_force, control, contact_link_dictionary):
    # This function will calculate the right hand side of the dynamics function
    D_q, CG_q_qdot, Jac_Full, Jac_Full_Trans = Dynamics_Matrices(sim_robot, contact_link_dictionary)
    Jac_Trans_Mul_Contact_Force = np.matmul(Jac_Full_Trans, contact_force)
    Dynamics_RHS = Jac_Trans_Mul_Contact_Force.copy()
    for i in range(0, DOF):
        if i>=6:
            Dynamics_RHS[i] = Dynamics_RHS[i] + control[i - 6]
    D_q_Acc = Dynamics_RHS - CG_q_qdot
    Acc = np.matmul(np.linalg.pinv(D_q), D_q_Acc)
    return Acc

def Control_2_Full_List(DOF, control_i):
    control_list = [0] * DOF
    for i in range(0, DOF - 6):
        control_list[6 + i] = control_i[i]
    return control_list

def Robot_Total_Mass(world):
    # This function is used to calculate the given robot's total mass
    sim_robot = world.robot(0)
    Total_Mass = 0.0
    for i in range(0, sim_robot.numLinks()):
        sim_robot_link_i_Mass_structure = sim_robot.link(i).getMass()
        sim_robot_link_i_mass = sim_robot_link_i_Mass_structure.getMass()
        Total_Mass = Total_Mass +  sim_robot_link_i_mass
    return Total_Mass

def Impact_Mapping_fn(world, final_state, treenode_child, contact_link_dictionary):
    """
        This function is used to calculate the robot state after the impact mapping.
        The two assumptions:
                            1. The collision is infinitesimal.
                            2. The collsiion is inelastic.
    """
    # q_new, qdot_new = Impact_Mapping_fn(world, State_Init, Root_Node, Contact_Link_Dictionary)

    sim_robot = world.robot(0)
    Robot_State_Update(sim_robot, final_state)
    D_q, CG_q_qdot, Jac_Full, Jac_Full_Trans = Dynamics_Matrices(sim_robot, contact_link_dictionary)
    DOF = len(D_q)
    q_old = final_state[0:DOF];    qdot_old = final_state[DOF:]
    Jac_Act_Array, Jac_Act_Trans_Array = Jac_Act_Selection(Jac_Full, contact_link_dictionary, treenode_child["Contact_Status"])
    J_D_inv_J_Trans = np.matmul(np.matmul(Jac_Act_Array, np.linalg.inv(D_q)), Jac_Act_Trans_Array)
    J_qdot_old = np.matmul(Jac_Act_Array, qdot_old)
    Negative_Impulse = np.matmul(np.linalg.pinv(J_D_inv_J_Trans), J_qdot_old)
    Impulse = np.dot(Negative_Impulse,-1)
    qdot_offset = np.matmul(np.matmul(np.linalg.inv(D_q), Jac_Act_Trans_Array), Impulse)
    qdot_new = qdot_old[:]
    for i in range(0, len(qdot_new)):
        qdot_new[i] = qdot_new[i] + qdot_offset[i]
    return q_old, qdot_new


def Jac_Act_Selection(jac_full_matrix_raw, contact_link_dictionary, contact_status_dictionary):
    # This function is used to select the active Jacobian matrix given the contact status
    # The sequence order of the jac_full_matrix is arranged accoridng to contact_link_dictionary.keys()

    Jac_Act_List = []
    contact_link_list = contact_link_dictionary.keys()
    jacobian_matrix_index = 0
    Jac_Act_Index_List = []
    for i in range(0, len(contact_link_list)):
        link_i_contact_status = contact_status_dictionary[contact_link_list[i]]
        for j in range(0, len(link_i_contact_status)):
            link_i_point_j_contact_status = link_i_contact_status[j]
            if link_i_point_j_contact_status == 1:
                Jac_Act_Index_List.append(jacobian_matrix_index)
            jacobian_matrix_index = jacobian_matrix_index + 1

    jac_full_matrix = jac_full_matrix_raw[:]
    # Based on the Jac_Act_Index_List, certain rows of jac_full_matrix have to be added.
    for i in Jac_Act_Index_List:
        Jac_Act_Start_Index = 3 * i
        Jac_Act_Final_Index = Jac_Act_Start_Index + 3
        for j in range(Jac_Act_Start_Index, Jac_Act_Final_Index):
            Jac_Act_List.append(jac_full_matrix[j])
    Jac_Act_Array = np.array(Jac_Act_List)
    Jac_Act_Trans_Array = a = copy.deepcopy(Jac_Act_Array.T)
    return Jac_Act_Array, Jac_Act_Trans_Array
