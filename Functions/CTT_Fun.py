import sys, os, copy
from random import randint
import numpy as np
import math

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
        y_val, y_type = Nodes_Optimization_ObjNConstraint(self.world, self.treenode_parent, self.treenode_child, self.contact_link_dictionary, self.terr_model, self.robot_option, self.grids, self.DOF, self.control_force_len, x)
        for i in range(0,len(y_val)):
            y[i] = y_val[i]

    def __callg__(self, x, y, G, row, col, rec, needg):
        # This function will be used if the analytic gradient is provided
        pass

def Nodes_Optimization_fn(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, Grids_Number):
    """
    This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
    """
    # The outer operation runs over the enumeration of the duration guess
    Duration_min = 0.1
    Duration_max = 2.5
    Duration_Grids_number = 41
    Duration_array = np.linspace(Duration_min, Duration_max, Duration_Grids_number)
    Duration_list = Duration_array.tolist()

    Opt_Flag = False
    for Duration_i in Duration_list:
        # Inner operation
        Opt_Soln, Opt_Flag, Final_State_List = Nodes_Optimization_Inner_Opt(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option,  Duration_i, Grids_Number, Duration_min, Duration_max)
        if Opt_Flag == True:
            return Opt_Soln, Opt_Flag, Final_State_List
    return Opt_Soln, Opt_Flag, Final_State_List

def Nodes_Optimization_Inner_Opt(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids, duration_min, duration_max):
    # The inner optimization of this contact transition tree
    Seed_Guess_List, DOF, Control_Force_Len = Seed_Guess_Gene(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, duration, grids)
    # Single_Robot_Motion_Plot(world, DOF, Control_Force_Len, contact_link_dictionary, terr_model, robot_option, grids, Seed_Guess_List)
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
    tau_bound = world.robot(0).getTorqueLimits()
    for i in range(0, grids):
        for j in range(0, len(tau_bound) - 6):
            xlb.append(-tau_bound[j + 6]);                  xub.append(tau_bound[j + 6])
    """
            Contact Force bounds: not bounded
    """
    force_bound = np.inf
    for i in range(0, grids):
        for j in range(0, Control_Force_Len):
            xlb.append(-force_bound);                       xub.append(force_bound)

    # Optimization problem setup
    Nodes_Optimization_Inner = Nodes_Optimization_Inner_Opt_Prob(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, Seed_Guess_List)
    Nodes_Optimization_Inner.xlb = xlb;                     Nodes_Optimization_Inner.xub = xub

    # This self structure is different from the previous self structure defined in the optimization problem
    y_val, y_type = Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, Seed_Guess_List)

    lb, ub = Constraint_Bounds(y_type)
    Nodes_Optimization_Inner.lb = lb;                       Nodes_Optimization_Inner.ub = ub
    cfg = snoptConfig()
    cfg.printLevel = 1;
    cfg.printFile = "result.txt"
    cfg.majorIterLimit = 3000
    cfg.minorIterLimit = 50000000
    cfg.iterLimit = 50000000

    cfg.addIntOption("Major print level", 1)
    cfg.addIntOption("Minor print level", 1)

    slv = solver(Nodes_Optimization_Inner, cfg)
    # rst = slv.solveRand()
    rst = slv.solveGuess(np.array(Seed_Guess_List))

    y_val, y_type = Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, Control_Force_Len, rst.sol)

    constraint_vio = Optimal_Solution_Validation(y_val, y_type)

    if constraint_vio<= 0.1:
        # Here 0.1 is threshold for constraint validity
        sol_valid_flag = True
    else:
        sol_valid_flag = False

    _, State_List_Array, _, _ = Seed_Guess_Unzip(rst.sol, DOF, Control_Force_Len, grids)

    # Get the final robot state out for potential impact mapping

    return rst.sol, sol_valid_flag, State_List_Array[grids - 1].tolist()

def Nodes_Optimization_ObjNConstraint(world, treenode_parent, treenode_child, contact_link_dictionary, terr_model, robot_option, grids, DOF, control_force_len, seed_guess_list):
    # This is the core optimization function of this contact transition tree project
    sim_robot = world.robot(0);
    contact_link_list = contact_link_dictionary.keys()
    contact_type, contact_link_cmp_dict, contact_link_itc_dict = TreeNode_Status_CMP(treenode_parent, treenode_child)

    # Constraint value and type
    y_val = [];                                                 y_type = []

    T, State_List_Array, Control_List_Array, Contact_Force_List_Array = Seed_Guess_Unzip(seed_guess_list, DOF, control_force_len, grids)
    """
        Objective function: Kinetic energy at the terminal time
    """
    if contact_type == 1:
        # Impact mapping kinetic energy
        Final_State_Config, Final_State_Velocity = Impact_Mapping_fn(world, State_List_Array[grids - 1], treenode_child, contact_link_dictionary)
        Final_State = np.append(Final_State_Config, Final_State_Velocity)
        # Final_State = Final_State_Config + Final_State_Velocity
    else:
        Final_State = State_List_Array[grids - 1]
    Robot_State_Update(sim_robot, Final_State)
    Final_KE = sim_robot.getKineticEnergy()
    y_val.append(Final_KE);                                     y_type.append(1)

    """
        Constraint 1:
                    Initial matching constraint condition
    """
    Parent_State = treenode_parent["State"];                    State_0 = State_List_Array[0]
    Initial_Match_Result = State_0 - Parent_State
    for Initial_Match_i in Initial_Match_Result:
        y_val.append(Initial_Match_i * Initial_Match_i)
        y_type.append(0)

    """
        Constraint 2:
                    Dynamics constraints at collocation point (middle point)
    """
    # Here T stands for the time duration between two next grids
    deltaT = T/(grids - 1.0)
    for i in range(0, grids - 1):
        # Robot state dynamics at the front side
        State_Front = State_List_Array[i]
        Robot_State_Update(sim_robot, State_Front)
        Contact_Force_Front = Contact_Force_List_Array[i]
        Control_Front = Control_List_Array[i]
        Acc_Front = AccFromTorque(sim_robot, contact_link_dictionary, Contact_Force_Front, Control_Front)

        # Robot state dynamics at the back side
        State_Back = State_List_Array[i + 1];
        Robot_State_Update(sim_robot, State_Back)
        Contact_Force_Back = Contact_Force_List_Array[i+1]
        Control_Back = Control_List_Array[i+1]
        Acc_Back = AccFromTorque(sim_robot, contact_link_dictionary, Contact_Force_Back, Control_Back)

        # Collocation point
        State_Mid, Acc_Mid = State_Mid_from_CubicSpline(DOF, deltaT, State_Front, Acc_Front, State_Back, Acc_Back)
        Robot_State_Update(sim_robot, State_Mid)

        # Here the inverse dynamics will be used to accelerate this process.
        Contact_Force_Mid = 0.5 * Contact_Force_Front + 0.5 * Contact_Force_Back
        Control_Mid = 0.5 * Control_Front + 0.5 * Control_Back;
        Dynamics_Constraint_i = CollocationDynamicsConstraint(sim_robot, Acc_Mid, contact_link_dictionary, Contact_Force_Mid, Control_Mid)

        List_Obj_Update(Dynamics_Constraint_i, 0, y_val, y_type)

        if (i == grids-1) and ((contact_type == 0) or (contact_type == -1)):
            # Inertial shaping finally stabilizes the robot
            List_Obj_Update(Acc_Back, 0, y_val, y_type)

    """
        Constraint 3:
                    Contact Maintenance!
                    The previous unchanged active contact constraint have to be satisfied
    """
    Ref_Contact_PosNVel = copy.deepcopy(treenode_parent["Contact_PosNVel"])
    for i in range(0, grids):
        # Robot state dynamics at the each grid i
        At_i_State = State_List_Array[i]
        # print At_i_State
        Robot_State_Update(sim_robot, At_i_State)
        At_i_Contact_Link_PosNVel = Contact_Link_PosNVel(sim_robot, contact_link_dictionary, -1)
        for j in range(0, len(contact_link_itc_dict)):
            contact_link_i_number = contact_link_list[j]
            At_i_Ref_Link_j_Dict = Ref_Contact_PosNVel[j]
            At_i_Opt_Link_i_Dict = At_i_Contact_Link_PosNVel[j]
            for k in contact_link_itc_dict[contact_link_i_number]:
                # This is the order of Ref_Contact_PosNVel/RefContact_PosNVel list
                # If not empty, there exists active pos/vel

                At_i_Ref_Link_j_Pos_k = Ref_Contact_PosNVel[j]["Pos"][k]
                At_i_Ref_Link_j_Vel_k = Ref_Contact_PosNVel[j]["Vel"][k]

                At_i_Opt_Link_j_Pos_k = At_i_Contact_Link_PosNVel[j]["Pos"][k]
                At_i_Opt_Link_j_Vel_k = At_i_Contact_Link_PosNVel[j]["Vel"][k]

                At_i_Ref_Opt_Link_j_Pos_k_Constraint = vectorops.sub(At_i_Ref_Link_j_Pos_k, At_i_Opt_Link_j_Pos_k)
                At_i_Ref_Opt_Link_j_Vel_k_Constraint = vectorops.sub(At_i_Ref_Link_j_Vel_k, At_i_Opt_Link_j_Vel_k)

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
                    * For inertial-shaping, both treenode_parent["Contact_Status"] and treenode_child["Contact_Status"] will work.
                    * For contact retract, now I consider it to be first-grid in contact but then the rest is inertia shaping!
                    * For contact addition, the last grid the robot will make a new contact.
    """
    Parent_Contact_Status = copy.deepcopy(treenode_parent["Contact_Status"])
    Child_Contact_Status = copy.deepcopy(treenode_child["Contact_Status"])
    if contact_type == 0:
        for i in range(0, grids):
            Contact_Force_i = Contact_Force_List_Array[i]
            Contact_Force_Constraint(Contact_Force_i, contact_link_list, Parent_Contact_Status, terr_model, y_val, y_type)
    else:
        if(contact_type == -1):
            for i in range(0, grids):
                Contact_Force_i = Contact_Force_List_Array[i]
                if(i ==0):
                    Contact_Force_Constraint(Contact_Force_i, contact_link_list, Parent_Contact_Status, terr_model, y_val, y_type)
                else:
                    Contact_Force_Constraint(Contact_Force_i, contact_link_list, Child_Contact_Status, terr_model, y_val, y_type)
        else:
            for i in range(0, grids):
                Contact_Force_i = Contact_Force_List_Array[i]
                if(i ==grids-1):
                    Contact_Force_Constraint(Contact_Force_i, contact_link_list, Child_Contact_Status, terr_model, y_val, y_type)
                else:
                    Contact_Force_Constraint(Contact_Force_i, contact_link_list, Parent_Contact_Status, terr_model, y_val, y_type)
    """
        Constraint 6:
                    Kinetic energy threshold constraint
                    Inertial shaping to minimize the terminal kinetic energy to be within a threshold
    """

    if ((contact_type == 0) or (contact_type == -1)):
        KE_Threshold = 0.1
        y_val.append(KE_Threshold - Final_KE);
        y_type.append(1)

    """
        Constraint 7: self-collision constraint to be added in future!
    """
    return y_val, y_type
