#!/usr/bin/python

# This is used to serve as a library such that some functions do not need to be rewritten over and ove again
import sys, os
from random import randint
import numpy as np
import math, ipdb

# This part is for the operation on list

def Cross_Product(a, b):
    # Only defined for 3D vectors
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c

def Dot_Product(list_a, list_b):
    # This is defined for list of arbitrary length
    if type(list_a) is list:
        n = len(list_a)
    else:
        # <type 'numpy.ndarray'>
        n = list_a.size
    list_a_dot_b = 0
    for i in range(0, n):
        list_a_dot_b = list_a_dot_b + list_a[i] * list_b[i]
    return list_a_dot_b

def List_Sum_fn(list_a, list_b):
    # This function is used to add two list up by element wise
    n = len(list_a)
    list_a_sum_b = []
    for i in range(0, n):
        list_a_sum_b.append(list_a[i] + list_b[i])
    return list_a_sum_b

def List_Minus_fn(list_a, list_b):
    # This function outpu the operation of list_a - list_b
    n = len(list_a)
    list_a_minus_b = []
    for i in range(0, n):
        list_a_minus_b.append(list_a[i] - list_b[i])
    return list_a_minus_b

def List_Norm_fn(given_list):
    # This function is used to calculate the norm of a given list
    n = len(given_list)
    given_list_sum_square = 0
    for i in range(0, n):
        given_list_sum_square = given_list_sum_square + given_list[i] * given_list[i]
    return math.sqrt(given_list_sum_square)

def List_Append_fn(list_a, list_b):
    # This function is used to append the list_b to the end of list_a
    list_a_ref = list_a[:]
    for i in range(0, len(list_b)):
        list_a_ref.append(list_b[i])
    return list_a_ref

def List_Flatten_fn(l):
    flat_list = [item for sublist in l for item in sublist]
    return flat_list

# This part of functions is used for file operation: read, write
def Read_Txt_fn(file_name):
    Empty_List = []
    with open(file_name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        for Txt_File_Str_i in Txt_File_Str:
            Txt_File_Str_i.replace("\'","")
            Empty_List.append(float(Txt_File_Str_i))
    return Empty_List

def Write_Txt_fn(file_name, list_array, path_name):
    file_path_name = path_name + file_name
    Txt_File = open(file_path_name, 'w')
    for list_array_i in list_array:
        print>>Txt_File, list_array_i
    Txt_File.close()

def Random_Color_Generator_fn():
    ret = []
    unit = 1.0/255.0
    r = randint(0,255)
    g = randint(0,255)
    b = randint(0,255)

    ret.append(r * unit)
    ret.append(g * unit)
    ret.append(b * unit)

    return ret

# Some dynamics functions
def Kinetic_Energy_fn(sim_robot, x):
    State_Tot_Number = len(x)
    Configuration_Number = State_Tot_Number/2
    D_q = sim_robot.getMassMatrix()
    D_q = np.array(D_q)
    qdot = np.array(x[Configuration_Number:])
    KE_val = np.dot(np.transpose(qdot), np.dot(D_q, qdot))
    return KE_val

def State_Loader_fn(*args):
    if len(args) == 2:
        # In this case, the robot is only given the configuration file
        config_file_path = args[1] + args[0]
        DOF, Config_Init = Configuration_Loader_fn(config_file_path)
        # Then the Velocty_Init is set to be a zero value list
        Velocity_Init = []
        for i in range(0,DOF):
            Velocity_Init.append(0)
    else:
        if len(args) == 3:
            config_file_path = args[2] + args[0]
            velocity_file_path = args[2] + args[1]
            Config_Init = Read_Txt_fn(config_file_path)
            Velocity_Init = Read_Txt_fn(velocity_file_path)
            DOF = len(Config_Init)
        else:
            raise RuntimeError("Input name should be either one config file or two txt files!")
    return DOF, Config_Init, Velocity_Init

def State_Writer_fn(*args):
    if len(args)==3:
        # In this case,the configuration is written to a file with .config
        Configuration_Writer_fn(args[0], args[1], args[2])
    else:
        if len(args) == 5:
            Write_Txt_fn(args[2], args[0], args[4])
            Write_Txt_fn(args[3], args[1], args[4])
        else:
            raise RuntimeError("Inputs should be either one config with string name or config and velocity with string names!")

def Configuration_Writer_fn(*args):
    # This function is used to write the given configuration into a file with .config
    file_list = args[0]
    file_name = args[1]
    path_name = args[2]
    file_path_name = path_name + file_name
    file_object  = open(file_path_name, 'w')
    DOF = len(file_list)
    file_object.write(str(DOF))
    file_object.write('\t')
    for i in range(0,DOF):
        file_object.write(str(file_list[i]))
        file_object.write(' ')
    file_object.close()
    return

def Configuration_Loader_fn(Config_Name):
    # This function is only used to load in the initial configuraiton
    # The initial file will be in the .config format
    with open(Config_Name,'r') as robot_angle_file:
        robotstate_angle_i = robot_angle_file.readlines()
    config_temp = [x.replace('\t',' ') for x in robotstate_angle_i]
    config_temp = [x.replace('\n','') for x in config_temp]
    config_temp = [float(i) for i in config_temp[0].split()]

    DOF = int(config_temp[0])
    # Config_Init = np.array(config_temp[1:])
    Config_Init = config_temp[1:]
    return DOF, Config_Init

def Robot_Link_Attributes_fn(sim_robot, link_index):
    # This function is used to output the necessary attributes for the comparison of dynamics
    # 4 things are to be compared
    # 1. The COM position
    # 2. The Velocity at COM
    # 3. The Angular velocity of this robot link
    # 4. The moment of Inertia of the link_index is robot

    robot_link_i = sim_robot.link(link_index)
    robot_link_i_mass_structure = robot_link_i.getMass()
    robot_link_i_mass_structure_COM_Local = robot_link_i_mass_structure.getCom()
    robot_link_i_mass_structure_Inertia = robot_link_i_mass_structure.getInertia ()     # Not sure about the meaning of the getInertia
    robot_link_i_mass_structure_COM_Global = robot_link_i.getWorldPosition(robot_link_i_mass_structure_COM_Local)
    robot_link_i_mass_structure_COM_Velo = robot_link_i.getPointVelocity(robot_link_i_mass_structure_COM_Local)
    robot_link_i_AngVel_Global = robot_link_i.getAngularVelocity()    # One thing about this Angular Velocity is that it is expressed in the global coordinate

    robot_link_i_axis = robot_link_i.getAxis()

    robot_link_i_Transform = robot_link_i.getTransform()
    robot_link_i_Rotation_Matrix = List_2_Mat_fn(robot_link_i_Transform[0])
    robot_link_i_Rotation_Matrix_Trans = List_Mat_Trans_fn(robot_link_i_Rotation_Matrix)

    # robot_link_i_AngVel_Local = List_Mat_Multi_List_Vec_fn(robot_link_i_Rotation_Matrix_Trans, robot_link_i_AngVel_Global)    # One thing about this Angular Velocity is that it is expressed in the global coordinate

    print ""
    print "For Link " + str(link_index) + " name (" + robot_link_i.getName() + ") "
    print "Global Ref position is: [" + List_to_String_fn(robot_link_i.getWorldPosition([0,0,0])) + "]"
    print "Local COM position is: [" + List_to_String_fn(robot_link_i_mass_structure_COM_Local) + "]"
    print "Global COM position is: [" + List_to_String_fn(robot_link_i_mass_structure_COM_Global) + "]"
    # print "Local COM position is: [" + List_to_String_fn(robot_link_i_mass_structure_COM_Local) + "]"
    # print "Local axis in Global frame is " + List_to_String_fn(robot_link_i_axis)
    print "Velocity at COM is : [" + List_to_String_fn(robot_link_i_mass_structure_COM_Velo) + "]"
    # print "Local Angular velocity of this link is : [" + List_to_String_fn(robot_link_i_AngVel_Local) + "]"
    print "Global Angular velocity of this link is : [" + List_to_String_fn(robot_link_i_AngVel_Global) + "]"
    print "Moment of Inertia of this link is : [" + List_to_String_fn(robot_link_i_mass_structure_Inertia) + "]"
    print "Rotation matrix of this link is : "
    Print_List_Mat_fn(robot_link_i_Rotation_Matrix_Trans)
    # ipdb.set_trace()

def Print_List_Mat_fn(list_matrix):
    # This function is used to print the list of list into an elegent format for visualization
    for col_vec in list_matrix:
        col_vec_string = " "
        for col_vec_ele in col_vec:
            col_vec_string = col_vec_string + str(col_vec_ele) + " "
        print col_vec_string

def List_2_Mat_fn(list_temp):
    # This function is used to convert a squared vector into a squared matrix
    n = len(list_temp)
    Dim = int(math.sqrt(n))
    Mat = []
    Mat_i = []

    for i in range(0, Dim):
        Start_Ind = i * Dim
        End_Ind = (i + 1) * Dim
        Mat_i = list_temp[Start_Ind:End_Ind]
        Mat.append(Mat_i)
    return Mat

def List_Mat_Trans_fn(list_matrix):
    # This function gives out the transposed matrix
    trans_list_matrix = [[list_matrix[j][i] for j in range(len(list_matrix))] for i in range(len(list_matrix[0]))]
    return trans_list_matrix

def List_Mat_Multi_List_Vec_fn(list_matrix, list_vec):
    # This function is used to multiply the list matrix with a list vector
    # According to the default manner in Klampt, the matrix is defined in a column way so each three elements stand for a column instead of a row
    Dim = len(list_vec)
    Dim_row = len(list_matrix)
    Dim_col = len(list_matrix[0])
    if Dim_row != Dim_col:
        raise RuntimeError("Input is not a square matrix!\n")
    Res_List_Vec = []
    for i in range(0, Dim):
        Res_List_i = 0
        for j in range(0, Dim):
            Res_List_i = Res_List_i + list_matrix[j][i] * list_vec[j]
        Res_List_Vec.append(Res_List_i)
    return Res_List_Vec

def List_to_String_fn(list_temp):
    # This function is used to converted the List into string to be printed
    list_string = ""
    for i in range(0, len(list_temp)):
        list_string_i = str(list_temp[i])
        if i == len(list_temp)-1:
            list_string = list_string + list_string_i
        else:
            list_string = list_string + list_string_i + ", "
    return list_string

def String_List_to_Number_List(str_list, string_name):
    # This function is used to convert the string list to a certain type of list
    if string_name =="float":
        res_list = [float(i) for i in str_list]
    else:
        res_list = [int(i) for i in str_list]
    return res_list

# Related to the operation of the robot contact link

def Contact_Link_Reader(File_Name, Path_Name):
    # This function is used to read-in the formation of the certain link and its associated contact points
    # The output of this function is a dictionary of several keys with multiple list values
    # File_Name = "./User_File/Contact_Link.txt"
    Contact_Link_Dictionary = dict()
    # The format of this function should be an integet with a list of contact points
    Link_Number_i = -1
    File_Path_Name = Path_Name +File_Name
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            if "Link" in Txt_File_Str_i:
                # This indicates a contact link number
                # Then the next few points would be the local coordinates of the contact extremities
                Txt_File_Str_i = Txt_File_Str_i.translate(None, 'Link')     # This step is to get the link number out
                Link_Number_i = int(Txt_File_Str_i)
                Contact_Link_Dictionary[Link_Number_i] = []
                continue
            if Link_Number_i>=0:
                # Here the Txt_File_Str_i is a string of 'a, b, c' format
                Txt_File_Str_i = Txt_File_Str_i.split(",")
                Txt_File_Flt_i = String_List_to_Number_List(Txt_File_Str_i,"float")
                Contact_Link_Dictionary[Link_Number_i].append(Txt_File_Flt_i)
    return Contact_Link_Dictionary

def Contact_Status_Reader(File_Name, Path_Name):
    # This function is used to read-in the formation of the certain link and its associated contact points
    # The output of this function is a dictionary of several keys with multiple list values
    # File_Name = "./User_File/Contact_Link.txt"
    Contact_Status_Dictionary = dict()
    # The format of this function should be an integet with a list of contact points
    Link_Number_i = -1
    File_Path_Name = Path_Name + File_Name
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            if "Link" in Txt_File_Str_i:
                # This indicates a contact link number
                # Then the next few points would be the local coordinates of the contact extremities
                Txt_File_Str_i = Txt_File_Str_i.translate(None, 'Link')     # This step is to get the link number out
                Link_Number_i = int(Txt_File_Str_i)
                Contact_Status_Dictionary[Link_Number_i] = []
                continue
            if Link_Number_i>=0:
                Contact_Status_Dictionary[Link_Number_i].append(int(Txt_File_Str_i))
    return Contact_Status_Dictionary

def Robot_State_Update(sim_robot, state_list):
    # This function is used to update the robot configuration and velocity
    DOF = len(state_list)/2
    config_list = state_list[0:DOF]
    velocity_list = state_list[DOF:]
    sim_robot.setConfig(config_list)
    sim_robot.setVelocity(velocity_list)

def Contact_Link_PosNVel(sim_robot, contact_link_dictionary, contact_index_i):
    # This function is used to retrieve the contact link point positions and velocities
    # If contact_index_i == -1 then it will give all the contact points and velocities

    Contact_Link_PosNVel_List = []
    if contact_index_i == -1:
        contact_link_index_list = contact_link_dictionary.keys()            # here the  keys is a  list of link numbers
        for i in range(0, len(contact_link_index_list)):
            link_local_extremities_i = contact_link_dictionary[contact_link_index_list[i]]
            Contact_Link_PosNVel_i = Contact_Link_i_PosNVel(sim_robot, contact_link_index_list[i], link_local_extremities_i)
            Contact_Link_PosNVel_List.append(Contact_Link_PosNVel_i)
    else:
        link_local_extremities_i = contact_link_dictionary[contact_index_i]
        Contact_Link_PosNVel_i = Contact_Link_i_PosNVel(sim_robot, contact_index_i, link_local_extremities_i)
        Contact_Link_PosNVel_List.append(Contact_Link_PosNVel_i)
    return Contact_Link_PosNVel_List

def Contact_Link_i_PosNVel(sim_robot, link_index, link_local_extremities):
    # This function is used to retrieve the robot global position and velocity of local extremities
    robot_link = sim_robot.link(link_index)
    robot_link_local_extremities_num = len(link_local_extremities)

    Contact_Link_PosNVel_i = dict()
    Contact_Link_PosNVel_i["Link"] = link_index
    Contact_Link_PosNVel_i["Pos"] = []
    Contact_Link_PosNVel_i["Vel"] = []

    for i in range(0, robot_link_local_extremities_num):
        link_local_extremity_i = link_local_extremities[i]
        link_local_extremity_i_Pos = robot_link.getWorldPosition(link_local_extremity_i)
        link_local_extremity_i_Vel = robot_link.getPointVelocity(link_local_extremity_i)
        Contact_Link_PosNVel_i["Pos"].append(link_local_extremity_i_Pos)
        Contact_Link_PosNVel_i["Vel"].append(link_local_extremity_i_Vel)
    return Contact_Link_PosNVel_i

def List_Obj_Update(list_value, constraint_type, y_val, y_type):
    # This function is used to update the list according to the function value
    for i in range(0, len(list_value)):
        y_val.append(list_value[i])
        y_type.append(constraint_type)

def ListOfLists_Print(list_of_lists):
    # This function is used to print out the element of list of lists
    NumberOfLists = len(list_of_lists)
    for i in range(0, NumberOfLists):
        print list_of_lists[i]

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
