import sys, os
from random import randint
import numpy as np
import math, ipdb

from OwnLib import *
from Terrain_Fun import *
from Visualization_Fun import *
from Initial_Setup_Opt import *
from Node_Fun import *
from Seed_Config_Fun import *
from Dynamics_Fun import *

# This file contains all the functions needed for the initialization of the seed for nonlinear optimization

# This file contains all the functions needed for the initialization of the seed for nonlinear optimization
def Seed_Guess_Gene(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids):
    """
    This function will generate the optimization variables needed for the further optimization
    """

    """
    The seed guess list will be of the following format

    seed_guess = [  duration,
                    State_List_1, ............, State_List_grids,
                    Ctrl_List_1, ............., Ctrl_List_grids,
                    Contact_Force_List_1, ...., Contact_Force_List_grids]
    """
    sim_robot = world.robot(0)

    # Get the length of the optimization variables
    State_List_Length = len(treenode_parent["State"]);      Robot_DOF = State_List_Length/2
    Ctrl_List_Length = State_List_Length - 6                  # All the internal joints are assumed to be actuated
    Contact_Force_Length = 0
    for contact_link_index_i in contact_link_dictionary.keys():
        contact_link_index_i_extremities = contact_link_dictionary[contact_link_index_i]
        for contact_link_index_i_extremity in contact_link_index_i_extremities:
            Contact_Force_Length = Contact_Force_Length + 1
    Contact_Force_List_Length = Contact_Force_Length * 3        # 3-dimensional force

    # First step is to calculate the reference robot state from the treenode_parent["State"]
    robot_init_state = treenode_parent["State"]
    robot_ref_state = robot_init_state[:]
    for i in range(0, Robot_DOF):
        state_i_acc = robot_init_state[Robot_DOF + i]/duration
        state_i_vel = robot_init_state[Robot_DOF + i]
        state_i_pos = robot_init_state[i] + state_i_vel * duration + 0.5 * state_i_acc * duration * duration
        robot_ref_state[i] = state_i_pos
        robot_ref_state[i + Robot_DOF] = 0

    # State_Writer_fn(robot_ref_state[0:Robot_DOF], "Ref_Config.config", robot_option)

    # Second step is to project this reference robot state into the constrained manifold
    # Goal_Config, Goal_Velocity = Seed_Guess_Gene_Robotstate(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_ref_state)
    # State_Writer_fn(Goal_Config, Goal_Velocity, "Goal_Config.txt", "Goal_Velocity.txt",robot_option)
    DOF, Goal_Config, Goal_Velocity = State_Loader_fn("Goal_Config.txt", "Goal_Velocity.txt", robot_option)
    # State_Writer_fn(Goal_Config, "Goal_Config.config", robot_option)

    # Third step is to construct a parabolic spline to connect the initial state to the goal configuration
    # For each DOF, its path should be of the following form: y(s) = a*s^2 + b*s + c
    DOF_Parabolic_Parameter_List = []
    for i in range(0, Robot_DOF):
        c_i = robot_init_state[i]
        b_i = robot_init_state[i + Robot_DOF] * duration
        a_i = Goal_Config[i] - b_i - c_i
        DOF_Parabolic_Parameter_List.append([a_i, b_i, c_i])
    path_s_array = np.linspace(0.0, 1.0, num=grids);    path_s_list = path_s_array.tolist()
    Pos_List, Vel_List, Acc_List = Parabolic_Para_All_Evaluation(DOF_Parabolic_Parameter_List, path_s_array, duration)

    """
    Initialization of the control and contact force with a pseudo inverse method
    """
    Contact_Force_List = []
    Control_List = []
    for i in range(0, grids):
        # The corresponding robot config, velocity, acceleration at i is
        Pos_i = Pos_List[i];        Vel_i = Vel_List[i];        Acc_i = Acc_List[i]
        # To get the system dynamics equation, please update the robot state at each iteration
        State_i = Pos_i + Vel_i

        Robot_State_Update(sim_robot, State_i)
        D_q_i, CG_q_qdot_i, Jac_Full_i, Jac_Full_Trans_i = Dynamics_Matrices(sim_robot, contact_link_dictionary)

        Dynamics_LHS_i = List_Sum_fn(List_Mat_Multi_List_Vec_fn(D_q_i, Acc_i),CG_q_qdot_i)
        Dynamics_RHS_Matrix_i = Dynamics_RHS_Matrix_fn(Robot_DOF, Jac_Full_Trans_i)
        Dynamics_RHS_Inverse_Matrix_i = np.linalg.pinv(Dynamics_RHS_Matrix_i)
        Contact_Force_N_Control_i = np.matmul(Dynamics_RHS_Inverse_Matrix_i, Dynamics_LHS_i)
        Contact_Force_i = Contact_Force_N_Control_i[0:Contact_Force_List_Length]
        Control_i = Contact_Force_N_Control_i[Contact_Force_List_Length:]

        Contact_Force_List.append(Contact_Force_i.tolist())
        Control_List.append(Control_i.tolist())

    # Contact_Force_List_Trans = map(list, zip(*Contact_Force_List))
    # Control_List_Trans = map(list, zip(*Control_List))
    # ListOfLists_Print(Contact_Force_List_Trans)
    # ListOfLists_Print(Control_List_Trans)
    # ipdb.set_trace()
    # Seed zip into a single list
    # ipdb.set_trace()
    Seed_Guess_List = []
    Seed_Guess_List.append(duration)                # The first element is the time duration
    for i in range(0, grids):
        Seed_Guess_List = Seed_Guess_List + Pos_List[i] + Vel_List[i]

    for i in range(0, grids):
        Seed_Guess_List = Seed_Guess_List + Control_List[i]

    for i in range(0, grids):
        Seed_Guess_List = Seed_Guess_List + Contact_Force_List[i]
    DOF = len(Goal_Config)
    return Seed_Guess_List, DOF, Contact_Force_List_Length

def Parabolic_Para_All_Evaluation(parabolic_parameter_list, path_s_list, duration):

    Pos_List = [];          Vel_List = [];              Acc_List = []
    for parabolic_parameter_list_i in parabolic_parameter_list:
        Pos_List_i, Vel_List_i, Acc_List_i = Parabolic_Para_Evaluation(parabolic_parameter_list_i, path_s_list, duration)
        Pos_List.append(Pos_List_i)
        Vel_List.append(Vel_List_i)
        Acc_List.append(Acc_List_i)

    Pos_List_New = map(list, zip(*Pos_List))
    Vel_List_New = map(list, zip(*Vel_List))
    Acc_List_New = map(list, zip(*Acc_List))
    return Pos_List_New, Vel_List_New, Acc_List_New

def Parabolic_Para_Evaluation(parabolic_parameter_list_i, path_s_list, duration):
    # This function is used to evaluate the position, velocity and acceleration given the parabolic parameter list and the s_list
    # The parabolic_parameter_list_i:[a_i. b_i. c_i]

    Pos_list = [];    Vel_list = [];    Acc_list = []
    a_i = parabolic_parameter_list_i[0]
    b_i = parabolic_parameter_list_i[1]
    c_i = parabolic_parameter_list_i[2]

    for i in range(0, len(path_s_list)):
        s_i = path_s_list[i]
        Pos_i = a_i * s_i * s_i + b_i * s_i + c_i
        Vel_i = (2 * a_i * s_i * s_i + b_i)/duration
        Acc_i = (2 * a_i)/(duration * duration)

        Pos_list.append(Pos_i)
        Vel_list.append(Vel_i)
        Acc_list.append(Acc_i)
    return Pos_list, Vel_list, Acc_list

def Seed_Guess_Unzip(seed_guess_list, DOF, contact_force_len, grids):
    # This function is used to get the robot state, control, contact force out of the seed_guess

    # The first element is the time duration
    duration = seed_guess_list[0]

    # Calculate the state, control and contact force number
    State_Number = 2 * DOF
    Total_State_Number = State_Number * grids
    Total_Control_Number = (DOF - 6) * grids
    Total_Contact_Force_Number = contact_force_len * grids

    # Get the state list
    Start_Index = 1
    State_List_All = seed_guess_list[Start_Index:(Start_Index + Total_State_Number)]
    Start_Index = Start_Index + Total_State_Number

    # Get the control list
    Control_List_All = seed_guess_list[Start_Index:(Start_Index + Total_Control_Number)]
    Start_Index = Start_Index + Total_Control_Number

    # Get the contact force list
    Contact_Force_List = seed_guess_list[Start_Index:]


    State_List_Array = np.reshape(State_List_All, (grids, 2 * DOF))
    Control_List_Array = np.reshape(Control_List_All, (grids, DOF-6))
    Contact_Force_List_Array = np.reshape(Contact_Force_List, (grids, contact_force_len))

    return duration, State_List_Array, Control_List_Array, Contact_Force_List_Array
