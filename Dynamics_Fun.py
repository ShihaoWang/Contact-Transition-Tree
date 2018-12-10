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

def Contact_Force_2_Full_List(DOF, contact_force):
    contact_force_list = [0] * DOF
    for i in range(0, DOF - 6):
        contact_force_list[6 + i] = contact_force[i]
    return contact_force_list

def Robot_Total_Mass(world):
    # This function is used to calculate the given robot's total mass
    sim_robot = world.robot(0)
    Total_Mass = 0.0
    for i in range(0, sim_robot.numLinks()):
        sim_robot_link_i_Mass_structure = sim_robot.link(i).getMass()
        sim_robot_link_i_mass = sim_robot_link_i_Mass_structure.getMass()
        Total_Mass = Total_Mass +  sim_robot_link_i_mass
    return Total_Mass
