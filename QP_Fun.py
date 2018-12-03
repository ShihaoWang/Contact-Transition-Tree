import sys, os
from random import randint
import numpy as np
import math, ipdb
from OwnLib import *

def Constraint_Force_Pyramid(beta_array):
    # This function is used to generate the pyramid approximation of the contact forces
    # beta_list is the magnitude of the corresponding components along each axis

    n = beta_array.size

    # Then a three-dimensional pyramid is generated for approximation

    # The first point is chosen to be

    Pyramid_Unit = np.zeros((3, n))
    Base_Unit = np.array([mu/np.sqrt(1 + mu * mu), 0, 1/np.sqrt(1 + mu * mu)])
    New_Unit_List = [0, 0, Base_Unit[2]]
    Pyramid_Unit[:,0] = Base_Unit
    Support_Angle = 2*pi/n
    New_Angle = Support_Angle
    for i in range(0, n - 1):
        # The only difference in those support vectors are the x and y components

        New_x = Base_Unit[0] * np.cos(New_Angle)
        New_y = Base_Unit[0] * np.sin(New_Angle)


        New_Unit_List[0] = New_x
        New_Unit_List[1] = New_y

        New_Unit = np.array(New_Unit_List)

        New_Angle = New_Angle + Support_Angle
        Pyramid_Unit[:,i+1] = New_Unit

    # At this step Pyramid_Unit is a 3 by n matrix with each column denoting the supporting vector in that direction
    Contact_Force = np.dot(Pyramid_Unit, beta_array)
    return Contact_Force, Pyramid_Unit

def QP_Controller(self):
    # This function takes in the self structure and computes the QP problem to stabilize the robot

    sim_robot = self.sim_robot

    sigma = self.sigma_t # In this case, the robot sigma is the initial sigma

    # These are four constrains to be added to this problem

    # The first constraint is the dynamics constraint
    D_q = sim_robot.getMassMatrix()
    D_q_sp = sparse.csc_matrix(np.array(D_q))
    Cons1_A = D_q_sp

    Cons1_B = sparse.csc_matrix((self.State_Number, self.Beta_Number * self.Cont_Force_Number))

    C_q_qdot = sim_robot.getCoriolisForces()
    C_q_qdot_sp = np.array(C_q_qdot)

    G_q = sim_robot.getGravityForces([ 0, 0, -9.8])
    G_q_sp = np.array(G_q)

    Tau_list = []
    for i in range(0, self.State_Number-6):
        Tau_list.append(-1)
    Tau_Data_List, Tau_Row_List, Tau_Col_List = Diagnal2DRC(Tau_list)

    B_up_sp = sparse.csc_matrix((6, self.State_Number-6))
    B_down_sp = sparse.csc_matrix((Tau_Data_List,(Tau_Row_List,Tau_Col_List)), shape = (self.State_Number-6,self.State_Number-6))
    Neg_B_sp = sparse.vstack([-B_up_sp, -B_down_sp]).tocsc()

    # So the first layer of the constraint equation is: D*q'' + C(q,q') = J^T * lamda + B * u

    Constraint_Pyramid = self.Constraint_Pyramid

    Jac_q = Contact_Jacobian_Matrix(sim_robot)
    Jac_q_array = np.array(Jac_q)
    Jac_q_T_array = np.transpose(Jac_q_array)
    Jac_q_T_Matrix = np.asmatrix(Jac_q_T_array)
    Neg_Jac_q_T_sp = sparse.csc_matrix(-Jac_q_T_array)

    Cons1_C = Neg_Jac_q_T_sp
    Cons1_D = Neg_B_sp

    Cons1 = sparse.hstack((Cons1_A, Cons1_B, Cons1_C, Cons1_D))

    l1 = - G_q_sp - C_q_qdot_sp
    u1 = - G_q_sp - C_q_qdot_sp

    # The second cobnstraint is the constraint force to beta coefficient

    Cons2_A = sparse.csc_matrix((3 * self.Cont_Force_Number, self.State_Number))

    Constraint_Pyramid = self.Constraint_Pyramid
    Constraint_Pyramid_sp = sparse.csc_matrix(Constraint_Pyramid)

    Cons2_B = Constraint_Pyramid_sp.copy()

    for i in range(0, self.Cont_Force_Number-1):
        Cons2_B = sparse.block_diag((Cons2_B, Constraint_Pyramid_sp), format='csc')

    Cons2_C = sparse.eye(3 * self.Cont_Force_Number)
    Cons2_D = sparse.csc_matrix((3* self.Cont_Force_Number, self.State_Number - 6))

    Cons2 = sparse.hstack((Cons2_A, Cons2_B, -Cons2_C, Cons2_D)).tocsc()
    l2 = np.zeros(3* self.Cont_Force_Number)
    u2 = np.zeros(3* self.Cont_Force_Number)

    # The thrid constraint is the beta feasibility
    Cons3_A = sparse.csc_matrix((self.Beta_Number * self.Cont_Force_Number, self.State_Number))
    Cons3_B = sparse.eye(self.Beta_Number * self.Cont_Force_Number)
    Cons3_C = sparse.csc_matrix((self.Beta_Number * self.Cont_Force_Number, 3 * self.Cont_Force_Number))
    Cons3_D = sparse.csc_matrix((self.Beta_Number * self.Cont_Force_Number, self.State_Number - 6))

    Cons3 = sparse.hstack((Cons3_A, Cons3_B, Cons3_C, Cons3_D)).tocsc()
    l3 = np.zeros(self.Beta_Number * self.Cont_Force_Number)
    u3 = np.inf*np.ones(self.Beta_Number * self.Cont_Force_Number)

    # The fourth constraint is the constraint force active constraint
    Cons4_A = sparse.csc_matrix((3 * self.Cont_Force_Number, self.State_Number))
    Cons4_B = sparse.csc_matrix((3 * self.Cont_Force_Number, self.Beta_Number * self.Cont_Force_Number))
    Cons4_C = Constraint_Force_Selection_Matrix(sigma)
    Cons4_D = sparse.csc_matrix((3 * self.Cont_Force_Number, self.State_Number - 6))
    Cons4 = sparse.hstack((Cons4_A, Cons4_B, Cons4_C, Cons4_D)).tocsc()

    l4 = np.zeros(3 * self.Cont_Force_Number)
    u4 = np.zeros(3 * self.Cont_Force_Number)

    # The five constraint is the torque limit constraint
    Cons5 = sparse.hstack((sparse.csc_matrix(( self.State_Number - 6, self.m - self.State_Number + 6)), sparse.eye(self.State_Number - 6))).tocsc()

    l5 = self.low_torque
    u5 = self.hgh_torque

    # The six constraint is the acceleration limit
    Cons6 = sparse.hstack((sparse.eye(self.State_Number), sparse.csc_matrix((self.State_Number, self.m - self.State_Number)))).tocsc()
    l6 = self.low_acc
    u6 = self.hgh_acc

    # The seventh constraint is the acceleration at the end effector
    Cons7_A_left = Jacobian_Selection_Matrix(sigma).todense()
    Cons7_A_right = np.array(Contact_Jacobian_Matrix(sim_robot))
    Cons7_A = np.dot(Cons7_A_left, Cons7_A_right)
    Cons7_B = sparse.csc_matrix((3 * self.Cont_Force_Number, self.m - self.State_Number))
    Cons7 = sparse.hstack((Cons7_A, Cons7_B)).tocsc()
    l7 = np.zeros(3 * self.Cont_Force_Number)
    u7 = np.zeros(3 * self.Cont_Force_Number)

    # The seventh constraint is at the acceleration which means that we would like the end effectors to have zero acceleration

    Cons = sparse.vstack((Cons1, Cons2, Cons3, Cons4, Cons5, Cons6, Cons7)).tocsc()
    l = np.hstack([l1, l2, l3, l4, l5, l6, l7])
    u = np.hstack([u1, u2, u3, u4, u5, u6, u7])

    # Now it is the test of the second objective function
    M_left = np.eye(36)
    M_right = np.zeros((36, 120))

    M = np.hstack((M_left, M_right))    # Here M is the selection matrix

    P = np.dot(np.transpose(M), M)
    P_sp = sparse.csc_matrix(P)


    K = Damping_Gain_Matrix(self)
    qdot = sim_robot.getVelocity()
    qdot_array = np.transpose(np.matrix(qdot))
    # ipdb.set_trace()

    q_right = np.matmul(K, qdot_array)

    q = np.matmul(np.transpose(M), q_right)

    # q = np.zeros(self.m)

    prob = osqp.OSQP()

    ipdb.set_trace()

    # Setup workspace
    prob.setup(P_sp, q, Cons, l, u)

    # Solve problem
    res = prob.solve()

    soln = res.x

    qddot = soln[0:self.State_Number]

    soln = soln[self.State_Number:]

    beta = soln[0:self.Beta_Number * self.Cont_Force_Number]
    soln = soln[self.Beta_Number * self.Cont_Force_Number:]

    lamda = soln[0:3 * self.Cont_Force_Number]
    soln = soln[3 * self.Cont_Force_Number:]

    tau = soln[:]

    if res.info.status == 'solved':
        flag = 1
        #  In this case, the QP problem is solved successfully

        # Constraint validation
        # 1. Dynamics constraint

        D_q_qddot = np.dot(D_q_sp.todense(), qddot)
        C_q_qdot_G_q = C_q_qdot_sp + G_q_sp
        Jac_T_Lamda = np.dot(Jac_q_T_Matrix, lamda)
        B_Tau = -np.dot(Neg_B_sp.todense(), tau)
        dyn_val = D_q_qddot + C_q_qdot_G_q - Jac_T_Lamda - B_Tau

        # 2. Beta coefficient to lamda
        lamda_beta = []
        for i in range(0, self.Cont_Force_Number):
            beta_i = beta[(i*self.Beta_Number):((i+1)*self.Beta_Number)]
            lamda_beta_i = np.dot(Constraint_Pyramid, beta_i)
            if i == 0:
                lamda_beta = lamda_beta_i.copy()
            else:
                lamda_beta = np.append(lamda_beta, lamda_beta_i)
    else:
        flag = 0
        #  In this case, the QP problem is not solved correctly


    return flag, qddot, beta, lamda, tau
