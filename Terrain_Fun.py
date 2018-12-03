import sys, os
import numpy as np
import math, ipdb
import ipdb
from OwnLib import *
# import xml.etree.ElementTree as ET
# from Main import *
from compiler.ast import flatten
from collections import Counter
import pickle

# This file containts the functions of the terrain definition and related functions

def Robot_Link_2_All_Terr_Dist(robot_link_point, terr_model):
    # This function computes the distance from the robot link to all terr_model
    terr_model_num = len(terr_model)
    Dist_2_Terr_Face_List = []
    Dist_2_Terr_Edge_List = []

    for i in range(0, terr_model_num):
        terr_model_i = terr_model[i]
        Dist_2_Terr_Face_i, Dist_2_Terr_Edge_i = Robot_Link_2_Terr_Dist(robot_link_point, terr_model_i)
        Dist_2_Terr_Face_List.append(Dist_2_Terr_Face_i)
        Dist_2_Terr_Edge_List.append(Dist_2_Terr_Edge_i)

    Dist_2_Face_Index = Dist_2_Terr_Face_List.index(min(Dist_2_Terr_Face_List))
    return Dist_2_Terr_Face_List[Dist_2_Face_Index], Dist_2_Terr_Edge_List[Dist_2_Face_Index]

def Robot_Link_2_Terr_Dist(robot_link_point, terr_model_i):

    # This function is used to calculate the relative distance between the robot given link to the environment obstacles
    # Due to the new method to define the terrain, the terrain is defined by a spatial plane

    # The input to this function: robot_link_point the global coordinate of a point on the robot link

    Dist = []
    terr_model_i_face_num = len(terr_model_i["normals"])        # Each terrain model can have more than one faces
    for i in range(0, terr_model_i_face_num):
        Terr_Normal_i = terr_model_i["normals"][i]              # Here Terr_Normal_i is a 3 by 1 list
        Terr2Point_i = [0, 0, 0]                                # Initialization
        Terr_Model_i_Vertex = terr_model_i["faces"][i][0]
        Terr2Point_i = List_Minus_fn(robot_link_point, Terr_Model_i_Vertex)
        Dist_i = Dot_Product(Terr2Point_i, Terr_Normal_i)
        Dist.append(Dist_i)
    Dist_2_Terr_Face = max(Dist)
    # In case of Many subfaces
    new_list = [k for k in Dist if k == Dist_2_Terr_Face]

    Face_Index = Dist.index(max(Dist))

    # The next job is to calculate the distance to edge to enforce the strict inclusion constraint
    Terr_Normal_i = terr_model_i["normals"][Face_Index]
    Terr_Faces_i = terr_model_i["faces"][Face_Index]

    Dist_2_Terr_Edge = Distance_2_Face_Edge(robot_link_point, Terr_Faces_i, Terr_Normal_i)

    return Dist_2_Terr_Face, Dist_2_Terr_Edge

def Distance_2_Face_Edge(robot_link_point, terr_face_vertices, terr_normal):
    # This function is used to calculate the signed distance of a given point to the edge of a triangle mesh
    # The three edge points are in the counter-clockwise way

    # The projected point on this face should be used for computation.
    point_left = terr_face_vertices[0]

    Vector_left2point = List_Minus_fn(robot_link_point, point_left)
    Projection_Length = Dot_Product(Vector_left2point, terr_normal)
    Projection_Vector = [Projection_Length*x for x in terr_normal]

    robot_link_point_projected = List_Minus_fn(robot_link_point, Projection_Vector)
    Dist_2_Edge = []

    for i in range(0, len(terr_face_vertices)):
        point_left = terr_face_vertices[i]
        if i == len(terr_face_vertices)-1:
            point_right = terr_face_vertices[0]
        else:
            point_right = terr_face_vertices[i+1]

        Dist_i = Distance2Edge(robot_link_point_projected, point_left, point_right, terr_normal)
        Dist_2_Edge.append(Dist_i)
    return min(Dist_2_Edge)

def Distance2Edge(robot_link_point, edge_a, edge_b, terr_normal):
    # This function is used to calculate the distance between a point and an edge
    # Here the edges are made in an counter-clock wise fashion
    edge_a2point = List_Minus_fn(robot_link_point, edge_a)
    edge_b2point = List_Minus_fn(robot_link_point, edge_b)
    double_area_vec = Cross_Product(edge_a2point, edge_b2point)
    double_area = List_Norm_fn(double_area_vec)
    edge_a2b = List_Minus_fn(edge_b, edge_a)
    edge_a2b_norm = List_Norm_fn(edge_a2b)
    dist2edge = double_area/(1.0 * edge_a2b_norm)
    dist2edge_sign = Dot_Product(double_area_vec, terr_normal)
    if dist2edge_sign>=0:
        return dist2edge
    else:
        return -dist2edge

def Terr_Model_Cal(world):
    # This function is used to calculate the necessary terrain model
    Terr_No = world.numTerrains()           # Number of Terrains
    Terr_Model = []                         # A list of dictionaries
    # ipdb.set_trace()
    for i in range(0, Terr_No):
        terr_i = world.terrain(i)
        terr_model_i = Terr_List_Analysis(terr_i)
        Terr_Model.append(terr_model_i)
    return Terr_Model

def Terr_List_Analysis(terr_i):
    # This function is used to analyze the terrain model
    # Here I have to be careful since the terr_model can be of two types: Point Could, Triangle Mesh or Geometric Primitive

    # The output of this function is a list of dictionary with key (faces, surface normal) and value (several 3D point, a 3D vector)
    terr_model_i = dict()
    terr_model_i["faces"] = []
    terr_model_i["normals"] = []

    terr_model_geometry = terr_i.geometry()
    terr_type = terr_model_geometry.type()
    if terr_type == 'TriangleMesh':
        terr_model_TriMesh = terr_model_geometry.getTriangleMesh()
        terr_model_vertices = terr_model_TriMesh.vertices               # This is the pure vertex coordinate
        terr_model_indices = terr_model_TriMesh.indices                 # This is the plane index formation
        # Here according to Swig's way to conduct the conversion between the Python and the low-level C++
        # Terr_model_vertices is of the std::vector<double> type, as a result the operation on the std::vector<double> can be directly applied to this object
        terr_model_vertices_list = []
        for i in range(0, terr_model_vertices.size()/3):                # Here the number 3 indicates that what we have here is a 3-d point
            terr_model_vertices_list.append([terr_model_vertices[3 * i], terr_model_vertices[3 * i + 1], terr_model_vertices[3 * i + 2]])

        # Here the reference points have already been put into one reference list
        OFF_Num_Vertices_In_Face = 3

        terr_model_face_number = terr_model_indices.size()/OFF_Num_Vertices_In_Face            # This is the number of faces in the terrain model according to the default way to define the OFF file

        for i in range(0, terr_model_face_number):
            # For each pair of indices, there is a certain set of surface normal to be computed
            terr_model_face_index_list_i = [terr_model_indices[OFF_Num_Vertices_In_Face * i], terr_model_indices[OFF_Num_Vertices_In_Face * i + 1], terr_model_indices[OFF_Num_Vertices_In_Face * i + 2]]
            terr_model_face_i = [terr_model_vertices_list[j] for j in terr_model_face_index_list_i]           # Each face_i is a set of 3 spatial points according to the OFF file's defnition of spatial terrain

            terr_model_normal_i_flag, terr_model_normal_i = Three_Points_2_Normal(terr_model_face_i)
            if terr_model_normal_i_flag == True:
                terr_model_i["faces"].append(terr_model_face_i)
                terr_model_i["normals"].append(terr_model_normal_i)

        # However, due to the existence of the redundant faces (2 triangle faces are 1 square face), a simplification is conducted to reduce the size of faces and normals
        terr_model_i_ref = terr_model_i.copy()          # Only the value
        list_index_abs = ListofList_Categorization(terr_model_i_ref["normals"])
        # Then the next job is to construct the new terrain with these information
        terr_model_i = Faces_Unification(list_index_abs, terr_model_i)
        return terr_model_i
    else:
        if terr_type == 'PointCloud':
            terr_model_PtCloud = terr_model_geometry.getPointCloud()
            raise RuntimeError('The preferred terrain type is Triangle Mesh!\n However, the given terrain file is Point Clould!\n')
        else:
            raise RuntimeError('The preferred terrain type is Triangle Mesh!\n However, the given terrain file is Geometric Primitive\n')

def Faces_Unification(list_index_abs_set, terr_model_i):
    # This function is used to categorize the triangle mesh into polygon shape to increase the efficiency
    terr_model_ref = dict()
    terr_model_ref["faces"] = []
    terr_model_ref["normals"] = []
    Potential_Poly_Num = len(list_index_abs_set)

    for i in range(0, Potential_Poly_Num):
        list_index_pair_i = list_index_abs_set[i]
        if type(list_index_pair_i) is list:
            # In this case, there are triangle mesh to be unified
            new_faces_indices_i = Faces_Unification_Inner(list_index_pair_i, terr_model_i)
            face_normal_i = terr_model_i["normals"][list_index_pair_i[0]]
            Terr_Model_Ref_Update(new_faces_indices_i, face_normal_i, terr_model_ref)
        else:
            # In this case, list_index_pair_i should be an integer value
            new_face_indices = terr_model_i["faces"][list_index_pair_i]
            face_normal_i = terr_model_i["normals"][list_index_pair_i]
            terr_model_ref["faces"].append(new_face_indices)
            terr_model_ref["normals"].append(face_normal_i)
    return terr_model_ref

def Terr_Model_Ref_Update(face_verticies, face_normal, terr_model_ref):
    # This function is used to update the terrain model given the normal and new face indicies
    for face_vertices_i in face_verticies:
        terr_model_ref["faces"].append(face_vertices_i)
        terr_model_ref["normals"].append(face_normal)

def Faces_Unification_Inner(some_list, terr_model_i):
    # This function is used to unifiy the triangle faces into square or higher segmental face if there exists
    face_number = len(some_list)
    some_list_ref = some_list[:]

    new_some_list = []

    while len(some_list_ref)>0:
        # The iteration is conducted if there is still unclassified element in some_list_ref
        list_face_index = some_list_ref.pop()               # This is the element to be unified with
        face_vertices_i = terr_model_i["faces"][list_face_index]
        Face_Verticies_Shifted =  Face_Unification_Sub(face_vertices_i, some_list_ref, terr_model_i)            # After this step,
        new_some_list.append(Face_Verticies_Shifted)
    return new_some_list

def Face_Unification_Sub(Face_Vertices, triangle_vertex_indices, terr_model_i):
    # This function is used to compare the face vertices with several triangle vertices
    triangle_vertex_indices_ref = triangle_vertex_indices
    for i in range(0, len(triangle_vertex_indices_ref)):
        Triangle_Vertices_Index_i = triangle_vertex_indices_ref[i]
        Triangle_Vertices_i = terr_model_i["faces"][Triangle_Vertices_Index_i]
        Subsub_Res, Face_Vertices_Shifted = Face_Unification_SubSub(Face_Vertices, Triangle_Vertices_i)
        if Subsub_Res == True:
            # An iterative process can be conducted here
            triangle_vertex_indices_ref.pop(i)
            Face_Vertices_Shifted = Face_Unification_Sub(Face_Vertices_Shifted, triangle_vertex_indices_ref, terr_model_i)
            return Face_Vertices_Shifted
    return Face_Vertices

def Face_Unification_SubSub(Face_Vertices, Triangle_Vertices):
    # This function is used to unifiy the Triangle_Vertices into Face_Vertices if there exists a shared boundary edge
    # Compare the same element in two lists
    Face_Vertices_Tuple = [tuple(t) for t in Face_Vertices]
    Triangle_Vertices_Tuple = [tuple(t) for t in Triangle_Vertices]

    shared_element_tuples = list(set(Face_Vertices_Tuple).intersection(Triangle_Vertices_Tuple))
    shared_element_lists = [list(elem) for elem in shared_element_tuples]

    cmp_rest = False
    if len(shared_element_lists)<2:
        return False, Face_Vertices
    else:
        # This means that there are at least two overlapping elements between these two lists
        # Therefore, it is needed to unify them together to reduce the unnecessary effort in figuring out the triangle mesh
        # The Face_Vertices is a reference list which remains to be unchanged. However, we would like to first change the order of Triangle_Vertices to match the shared elements in Face_Vertices

        Face_Shift_Res, Face_Vertices_Shifted = Face_List_Reording(Face_Vertices, shared_element_lists)
        Tria_Shift_Res, Triangle_Vertices_Shifted = Face_List_Reording(Triangle_Vertices, shared_element_lists)

        # Then the job is very straightforward, which is to import the 1st coordinate from the triangle vertices into the Face_Vertices to finish this unification
        Face_Vertices_Shifted.append(Face_Vertices_Shifted[-1])
        Face_Vertices_Shifted[-2] = Triangle_Vertices_Shifted[0]
        return True, Face_Vertices_Shifted

def Face_List_Reording(face_list, shared_elements):
    # This function is used to reorder the list to make sure that the shared elements will be at the last two element of the whole list
    # We have already known that the face vertices are ordered in a counter-clockwise order
    # The shared elements are in the face_list, the main idea is to shift the face_list's order such that the shared elements rest in the last two spots.
    Element_Number = len(face_list)
    i = 0
    face_list_ref = face_list[:]
    # print face_list_ref
    while i < Element_Number:
        List_Shift2Left(face_list_ref, 1)       # One to the left at each steps
        # print face_list_ref
        face_cmp_elements = face_list_ref[-2:]                  # Get the last two elements out
        # print face_cmp_elements
        # print shared_elements
        face_shared_elements = [x for x in face_cmp_elements if x in shared_elements]
        # print face_shared_elements
        if len(face_shared_elements) == 2:
            # In this case, the list should output a successful value
            return True, face_list_ref
        else:
            # In this case, the list should be further shifted
            i = i + 1
    return False, face_list_ref

def List_Shift2Left(list_i, n):
    """
    Shifts the lst over by n indices
    """
    if n < 0:
        raise ValueError('n must be a positive integer')
    if n > 0:
        list_i_first_element = list_i.pop(0)
        list_i.append(list_i_first_element)
        List_Shift2Left(list_i, n-1)  # repea

def ListofList_Categorization(list_of_list):
    # This function is used to categorize the a list of list into the a sorted index
    list_number = len(list_of_list)
    list_index = []
    for list_i in list_of_list:
        list_index_i = [i for i, x in enumerate(list_of_list) if x == list_i]               #### This is the most inefficient step and waits to be accelerated
        list_index.append(list_index_i)

    list_index_tuples = [tuple(t) for t in list_index]
    list_index_abs = list(set(list_index_tuples))
    list_index_abs = [list(elem) for elem in list_index_abs]

    return list_index_abs

def Three_Points_2_Normal(points_list):
    # This function computes the surface normal using the three spatial points
    # This surface normal should be computed with a counter-clock wise direction
    if len(points_list) != 3:
        raise RuntimeError('The length of the input point list is not equal to 3!\n Please change the point list length to calculate the surface normal vector')
    else:
        # In this case, these three points should be in the same plane.
        # The good thing about the OFF file is that the points_list is already defined into counter-clock way to facilitate the computation of surface normal
        normal_flag, normal_i = Terr_Normal_Computation_fn(points_list)
        return normal_flag, normal_i

def Terr_Normal_Computation_fn(terr_faces_i):
    # This function is used to calculate the surface normal given the terr_faces_i
    terr_i = terr_faces_i
    Point_A = terr_i[0]
    Point_B = terr_i[1]
    Point_C = terr_i[2]

    Point_A2B = [Point_B[0] - Point_A[0], Point_B[1] - Point_A[1], Point_B[2] - Point_A[2]]
    Point_B2C = [Point_C[0] - Point_B[0], Point_C[1] - Point_B[1], Point_C[2] - Point_B[2]]
    Normal_Flag = True

    normal_i = Cross_Product(Point_A2B, Point_B2C)

    normal_length_i = math.sqrt(normal_i[0] * normal_i[0] + normal_i[1] * normal_i[1] + normal_i[2] * normal_i[2])
    if normal_length_i ==0:
        return False, Point_A
    else:
        normal_i[0] = normal_i[0]/normal_length_i
        normal_i[1] = normal_i[1]/normal_length_i
        normal_i[2] = normal_i[2]/normal_length_i
        return True, normal_i
