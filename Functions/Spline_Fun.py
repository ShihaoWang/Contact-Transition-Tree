import sys, os
from random import randint
import numpy as np
import math, ipdb

# This function is used to save the spline functions

def State_Mid_from_CubicSpline(DOF, T, state_front, Acc_front, state_back, Acc_back):
    # This function is used to get the middle state from the State_Front, State_Back, Acc_Front and Acc_Back
    # Here T is the time duration between these two sequential grids
    Pos_front = state_front[0:DOF];     Vel_front = state_front[DOF:]
    Pos_back = state_back[0:DOF];       Vel_back = state_back[DOF:]


    State_mid = [0] * (2 * DOF);        Acc_mid = [0] * DOF
    # ipdb.set_trace()

    for i in range(0, DOF):
        Pos_front_i = Pos_front[i];     Vel_front_i = Vel_front[i];     Acc_front_i = Acc_front[i]
        Pos_back_i = Pos_back[i];       Vel_back_i = Vel_back[i];       Acc_back_i = Acc_back[i]

        Pos_mid_i = CubicSpline_Evaluation_fn(0, T, Pos_front_i, Pos_back_i, Vel_front_i,  Vel_back_i, 0.5)
        Vel_mid_i = CubicSpline_Evaluation_fn(0, T, Vel_front_i, Vel_back_i, Acc_front_i,  Acc_back_i, 0.5)
        Acc_mid_i = CubicSpline_Evaluation_fn(1, T, Vel_front_i, Vel_back_i, Acc_front_i,  Acc_back_i, 0.5)
        State_mid[i] = Pos_mid_i;       State_mid[i + DOF] = Vel_mid_i
        Acc_mid[i] = Acc_mid_i
    return State_mid, Acc_mid


def CubicSpline_Coeff_fn(T, x_init, x_end, xdot_init, xdot_end):
    """
        This function is used to calcualte the coefficients for the cubic spline
	    The cubic spline is expressed to be : y(s) = a*s^3 + b*s^2 + c*s + d
    """
    a = 2 * x_init - 2 * x_end + T * xdot_end + T * xdot_init
    b = 3 * x_end - 3 * x_init - T * xdot_end - 2 * T  * xdot_init
    c = T * xdot_init
    d = x_init
    return a, b, c, d

def CubicSpline_Evaluation_fn(order, T, x_init, x_end, xdot_init, xdot_end, s):
    # Here s is the path parameter not the time!
    a, b, c, d = CubicSpline_Coeff_fn(T, x_init, x_end, xdot_init, xdot_end)
    if order == 0:
        # zero order
        y = a * s * s * s + b * s * s + c * s + d
    else:
        # 1st order
        y = (3 * a * s * s + 2 * b * s + c)/T
    return y
