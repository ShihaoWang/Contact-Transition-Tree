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
from Spline_Fun import *
from Contact_Force_Fun import *

# This files contains all the files needed for the Contact Transition Tree (CTT) optimization

class Nodes_Optimization_Inner_Opt_Prob(probFun):
    # This is the class type used for the initial condition validation optimization problem
    def __init__(self, world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, control_force_len, seed_guess_list):

        self.world = world;
        self.treenode_parent = treenode_parent;                         self.treenode_child = treenode_child;
        self.contact_link_dictionary = contact_link_dictionary;         self.terr_model = terr_model
        self.robot_option = robot_option;                               self.grids = grids
        self.DOF = DOF;                                                 self.control_force_len = control_force_len
        self.seed_guess_list = seed_guess_list

        # Here the derivative is set to be finite difference method
        self.grad = False
        nx = len(seed_guess_list)
        y_val, y_type = Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, control_force_len, seed_guess_list)
        nc = len(y_type)

        probFun.__init__(self, nx, nc)

    def __callf__(self, x, y):
        # This function is the main function to calculate the objective and constraint function value
        y_val, y_type = Nodes_Optimization_ObjNConstraint(self.world, self.treenode_parent, self.treenode_child, self.contact_link_dictionary, self.terr_model, self.robot_option, self.grids, self.DOF, self.control_force_len, self.seed_guess_list)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

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
        ipdb.set_trace()
        Nodes_Optimization_Inner_Opt(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option,  Duration_i, Grids_Number, Duration_min, Duration_max)

def Nodes_Optimization_Inner_Opt(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids, duration_min, duration_max):
    # The inner optimization of this contact transition tree
    Seed_Guess_List, DOF, Control_Force_Len = Seed_Guess_Gene(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids)

    # Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, Seed_Guess_List)

    xlb = [];                                               xub = []
    # The seed_guess_list's format determines the bounds on the optimized variables
    xlb.append(duration_min);                               xub.append(duration_max)
    """
            Robot State bounds
    """
    # The next part is the grids number of robot state
    qmin, qmax = world.robot(0).getJointLimits();           dqmax_val = world.robot(0).getVelocityLimits()
    for i in range(0, grids):
        xlb = List_Append_fn(xlb, qmin);                    xub = List_Append_fn(xub, qmax)
        for j in range(0, len(dqmax_val)):
            xlb.append(-dqmax_val[j]);                      xub.append(dqmax_val[j])
    """
            Control bounds
    """
    tau_bound = world.robot(0).getTorqueLimits ()
    for i in range(0, grids):
        for j in range(0, len(tau_bound) - 6):
            xlb.append(-tau_bound[j + 6]);                  xub.append(tau_bound[j + 6])
    """
            Contact Force bounds: not bounded
    """

    # Optimization problem setup
    Nodes_Optimization_Inner = Nodes_Optimization_Inner_Opt_Prob(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, Seed_Guess_List)
    Nodes_Optimization_Inner.xlb = xlb;                     Nodes_Optimization_Inner.xub = xub

    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, Seed_Guess_List)
    lb, ub = Constraint_Bounds(y_type)
    Nodes_Optimization_Inner.lb = lb;                       Nodes_Optimization_Inner.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 1;
    # cfg.printFile = "result.txt"
    cfg.majorIterLimit = 300
    cfg.minorIterLimit = 50000
    # ipdb.set_trace()

    slv = solver(Nodes_Optimization_Inner, cfg)
    # rst = slv.solveRand()
    rst = slv.solveGuess(np.array(Seed_Guess_List))

    return rst

def Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, control_force_len, seed_guess_list):
    # This is the core optimization function of this contact transition tree project
    sim_robot = world.robot(0);
    # Constraint value and type
    y_val = [];                                                 y_type = []

    Duration, State_List_Array, Control_List_Array, Contact_Force_List_Array = Seed_Guess_Unzip(seed_guess_list, DOF, control_force_len, grids)

    """
        Objective function: Kinetic energy at the terminal time
    """
    Final_State = State_List_Array[grids - 1]
    Robot_State_Update(sim_robot, Final_State)
    Final_KE = Kinetic_Energy_fn(sim_robot, Final_State)
    y_val.append(Final_KE);                                     y_type.append(1)

    """
        Constraint 1:
                    Initial matching constraint condition
    """
    Parent_State = treenode_parent["State"];                    State_0 = State_List_Array[0]
    Initial_Match_Constraint = State_0 - Parent_State
    List_Obj_Update(Initial_Match_Constraint, 0, y_val, y_type)

    """
        Constraint 2:
                    Dynamics constraints at collocation point (middle point)
    """
    # Here T stands for the time duration between two next grids
    T = Duration/(grids - 1.0)

    contact_link_list = contact_link_dictionary.keys()
    contact_type, contact_link_cmp_dict, contact_link_itc_dict = TreeNode_Status_CMP(treenode_parent, treenode_child)

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
        # ipdb.set_trace()
        # Collocation point
        State_Mid, Acc_Mid = State_Mid_from_CubicSpline(DOF, T, State_Front, Acc_Front, State_Back, Acc_Back)
        Robot_State_Update(sim_robot, State_Mid)
        D_q_Mid, CG_q_qdot_Mid, Jac_Full_Mid, Jac_Full_Trans_Mid = Dynamics_Matrices(sim_robot, contact_link_dictionary)
        Control_Mid_Act = 0.5 * Control_Front + 0.5 * Control_Back;         Control_Mid = Contact_Force_2_Full_List(DOF, Control_Mid_Act)
        Contact_Force_Mid = 0.5 * Contact_Force_Front + 0.5 * Contact_Force_Back

        Dynamics_LHS_i = List_Sum_fn(List_Mat_Multi_List_Vec_fn(D_q_Mid, Acc_Mid),CG_q_qdot_Mid)
        Dynamics_RHS_i = List_Sum_fn(np.matmul(Jac_Full_Trans_Mid, Contact_Force_Mid),CG_q_qdot_Mid)
        Dynamics_Constraint_i = List_Minus_fn(Dynamics_LHS_i, Dynamics_RHS_i)
        List_Obj_Update(Dynamics_Constraint_i, 0, y_val, y_type)

        if (i == grids-1) and (contact_type == 0):
            # Inertial shaping finally stabilizes the robot
            List_Obj_Update(Acc_Back, 0, y_val, y_type)

    """
        Constraint 3:
                    Contact Maintenance!
                    The previous unchanged active contact constraint have to be satisfied
    """
                                   # list of dictionaries according to contact_link_list
    Ref_Contact_PosNVel = treenode_parent["Contact_PosNVel"]
    for i in range(0, grids):
        # Robot state dynamics at the each grid i
        At_i_State = State_List_Array[i];
        Robot_State_Update(sim_robot, At_i_State)
        At_i_Contact_Link_PosNVel = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)
        for j in range(0, len(contact_link_itc_dict)):
            contact_link_i_number = contact_link_list[j]
            At_i_Ref_Link_j_Dict = Ref_Contact_PosNVel[j]
            At_i_Opt_Link_i_Dict = At_i_Contact_Link_PosNVel[j]
            for k in contact_link_itc_dict[contact_link_i_number]:
                # This is the order of Ref_Contact_PosNVel/RefContact_PosNVel list
                # If not empty, there exists active pos/vel
                # ipdb.set_trace()

                At_i_Ref_Link_j_Pos_k = Ref_Contact_PosNVel[j]["Pos"][k]
                At_i_Ref_Link_j_Vel_k = Ref_Contact_PosNVel[j]["Vel"][k]

                At_i_Opt_Link_j_Pos_k = At_i_Contact_Link_PosNVel[j]["Pos"][k]
                At_i_Opt_Link_j_Vel_k = At_i_Contact_Link_PosNVel[j]["Vel"][k]

                At_i_Ref_Opt_Link_j_Pos_k_Constraint = List_Minus_fn(At_i_Ref_Link_j_Pos_k, At_i_Opt_Link_j_Pos_k)
                At_i_Ref_Opt_Link_j_Vel_k_Constraint = List_Minus_fn(At_i_Ref_Link_j_Vel_k, At_i_Opt_Link_j_Vel_k)

                List_Obj_Update(At_i_Ref_Opt_Link_j_Pos_k_Constraint, 0, y_val, y_type)
                List_Obj_Update(At_i_Ref_Opt_Link_j_Vel_k_Constraint, 0, y_val, y_type)
    """
        Constraint 4:
                    Terminal contact addition if there exists
    """
    if contact_type == 1:
        Final_State = State_List_Array[grids-1];
        Robot_State_Update(sim_robot, Final_State)
        Final_Contact_Link_PosNVel = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)
        for i in range(0, len(contact_link_list)):
            Final_Contact_PosNVel_Link_i = Final_Contact_Link_PosNVel[i]
            Final_Contact_PosNVel_Link_i_Pos = Final_Contact_PosNVel_Link_i["Pos"]
            for j in contact_link_cmp_dict[contact_link_list[i]]:
                Final_Contact_PosNVel_Link_i_Pos_j = Final_Contact_PosNVel_Link_i_Pos[j]
                Final_Contact_PosNVel_Link_i_Pos_j_Dist, _, _ = Robot_Link_2_All_Terr_Dist(Final_Contact_PosNVel_Link_i_Pos_j, terr_model)
                y_val.append(Final_Contact_PosNVel_Link_i_Pos_j_Dist)
                y_type.append(0)

    """
        Constraint 5:
                    Contact Force: Complementarity and Feasibility
                                   1. The net force should lie within the friction cone.
                                   2. The inactive contact points should not have nonzero forces there.
    """
    Transien_Contact_Status = treenode_parent["Contact_Status"]
    Terminal_Contact_Status = treenode_child["Contact_Status"]
    # ipdb.set_trace()
    for i in range(0, grids):
        At_i_Contact_Force = Contact_Force_List_Array[i]            # This contact force corresponds to the contat point Jacobian matrix
        # Here the contact force is a column vector
        # The order of the contact force obeys the order of the contact link list
        if i<grids-1:
            Contact_Force_Constraint(At_i_Contact_Force, contact_link_list, Transien_Contact_Status, terr_model, y_val, y_type)
        else:
            Contact_Force_Constraint(At_i_Contact_Force, contact_link_list, Terminal_Contact_Status, terr_model, y_val, y_type)

    """
        Constraint 6:
                    Kinetic energy threshold constraint
                    Inertial shaping to minimize the terminal kinetic energy to be within a threshold
    """
    if contact_type == 0:
        KE_Threshold = 0.1
        y_val.append(KE_Threshold - Final_KE);                                     y_type.append(1)

    """
        Constraint 7: self-collision constraint to be added in future!
    """
    # ipdb.set_trace()
    return y_val, y_type
