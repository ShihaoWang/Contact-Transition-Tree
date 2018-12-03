#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include "snopt.hh"
#include "Unified_Header.h"

#include <fstream>
#include <cmath>
#include <dlib/matrix.h>
#include "snoptProblem.hh"
#include <algorithm>

// This function is used to address the 3-dimensional fall mitigation project

// Pre-define of the bounds
double Inf = 1.1e20;										double PI = 3.1415926535897932384626;									   							double pi = PI;
// There are three types of variables in this optimization problem: robot state, control torques and contact forces
double rIxlow = -Inf;                  	double rIxupp = Inf;					double rIylow = -Inf;                  	double rIyupp = Inf;
double thetalow = -pi;                 	double thetaupp = pi;					double q1low = -2.3562;                	double q1upp = 0.733;
double q2low = 0.0;                    	double q2upp = 2.618;					double q3low = -1.3;                   	double q3upp = 0.733;
double q4low = -2.3562;                	double q4upp = 0.733;					double q5low = 0.0;                    	double q5upp = 2.618;
double q6low = -1.3;                   	double q6upp = 0.733;					double q7low = -3.14;                  	double q7upp = 1.047;
double q8low = -2.391;                 	double q8upp = 0.0;						double q9low = -3.14;                  	double q9upp = 1.047;
double q10low = -2.391;                	double q10upp = 0.0;
double AngRateMag = 3.0;			   				double AngRateLow = -AngRateMag;     	double AngRateHgh = AngRateMag;

double rIxdotlow = -Inf;               	double rIxdotupp = Inf;					double rIydotlow = -Inf;               	double rIydotupp = Inf;
double thetadotlow = -Inf;             	double thetadotupp = Inf;				double q1dotlow = AngRateLow;          	double q1dotupp = AngRateHgh;
double q2dotlow = AngRateLow;          	double q2dotupp = AngRateHgh;			double q3dotlow = AngRateLow;          	double q3dotupp = AngRateHgh;
double q4dotlow = AngRateLow;          	double q4dotupp = AngRateHgh;			double q5dotlow = AngRateLow;          	double q5dotupp = AngRateHgh;
double q6dotlow = AngRateLow;          	double q6dotupp = AngRateHgh;			double q7dotlow = AngRateLow;          	double q7dotupp = AngRateHgh;
double q8dotlow = AngRateLow;          	double q8dotupp = AngRateHgh;			double q9dotlow = AngRateLow;          	double q9dotupp = AngRateHgh;
double q10dotlow = AngRateLow;         	double q10dotupp = AngRateHgh;

double tau1_max = 100;             			double tau2_max = 100;					double tau3_max = 100;					double tau4_max = 100;
double tau5_max = 100;             			double tau6_max = 100;					double tau7_max = 60;              		double tau8_max = 50;
double tau9_max = 60;             			double tau10_max = 50;

double Acc_max = 10;										// This value is for the depth 3 left hand contact
// actually in the new method to formulate the optimization problem, the acceleration is left unbounded.
// double Acc_max = 10;

dlib::matrix<double> xlow_vec;							dlib::matrix<double> xupp_vec;
dlib::matrix<double> ctrl_low_vec;					dlib::matrix<double> ctrl_upp_vec;
dlib::matrix<double>  Envi_Map;							dlib::matrix<double> Envi_Map_Normal, Envi_Map_Tange; // The Normal and tangential vector of the plane

/**
 * Some global values are defined
 * Description
 */

double mini = 0.025;			int Grids = 8;			double mu = 0.35;

int Variable_Num = 48 * Grids + 1;

double Time_Seed; 															// This value will be adaptively changed to formulate an optimal solution
std::vector<Tree_Node_Ptr> All_Nodes;						// All nodes are here!
std::vector<Tree_Node_Ptr> Children_Nodes;			// All children nodes!
std::vector<Tree_Node_Ptr> Frontier_Nodes;			// Only Frontier ndoes!
std::vector<double> Frontier_Nodes_Cost;		    // The kinetic energy of each nodes

const int NumOfState = 22;

void D_q_C_q_qdot_Body(dlib::matrix<double> D_q, dlib::matrix<double> C_q_qdot)
{
	// This function is used to get out the share of the Body link to the whole dynamics motion
	// So besides the D(q), the corresponding C(q, qdot) will also be returned

	// Currently I am not pretty sure that whether it is a good idea to make all the needed variables into global
	// Matrix initialization

	D_q_t = dlib::zeros_matrix<double>(NumOfState, NumOfState);
	C_q_qdot_t = dlib::zeros_matrix<double>(NumOfState, 1);

	// Evaluation the mass matrix
	D_q_t(0,0) = 7.474E-1;
	D_q_t(1,1) = 7.474E-1;
	D_q_t(2,2) = 7.474E-1;
	D_q_t(3,3) = sin(beta_t)*(-1.0/4.0E2);
	D_q_t(3,4) = cos(beta_t)*sin(gamma_t)*(1.0/8.0E2);
	D_q_t(3,5) = cos(beta_t)*cos(gamma_t)*5.0E-4+1.0/8.0E2;
	D_q_t(4,3) = cos(beta_t)*sin(gamma_t)*(1.0/8.0E2);
	D_q_t(4,4) = cos(gamma_t)*(1.0/4.0E2);
	D_q_t(4,5) = sin(gamma_t)*(-5.0E-4);
	D_q_t(5,3) = cos(beta_t)*cos(gamma_t)*5.0E-4+1.0/8.0E2;
	D_q_t(5,4) = sin(gamma_t)*(-5.0E-4);

	C_q_qdot_t(2,0) = 7.331994;
	C_q_qdot_t(3,0) = -betadot_t*(alphadot_t*cos(beta_t)*(1.0/4.0E2)+gammadot_t*cos(gamma_t)*sin(beta_t)*5.0E-4+betadot_t*sin(beta_t)*sin(gamma_t)*(1.0/8.0E2))+gammadot_t*(betadot_t*cos(beta_t)*cos(gamma_t)*(1.0/8.0E2)-gammadot_t*cos(beta_t)*sin(gamma_t)*5.0E-4);
	C_q_qdot_t(4,0) = -gammadot_t*(gammadot_t*cos(gamma_t)*5.0E-4+betadot_t*sin(gamma_t)*(1.0/4.0E2)-alphadot_t*cos(beta_t)*cos(gamma_t)*(1.0/8.0E2))+(alphadot_t*alphadot_t)*cos(beta_t)*(1.0/8.0E2)+alphadot_t*gammadot_t*cos(gamma_t)*sin(beta_t)*5.0E-4;
	C_q_qdot_t(5,0) = betadot_t*(betadot_t*sin(gamma_t)*(1.0/8.0E2)-alphadot_t*cos(beta_t)*cos(gamma_t)*(1.0/8.0E2))-alphadot_t*betadot_t*cos(gamma_t)*sin(beta_t)*5.0E-4;
}

void D_q_C_q_qdot_LeftShoulder(dlib::matrix<double> D_q, dlib::matrix<double> C_q_qdot)
{
	D_q_t = dlib::zeros_matrix<double>(NumOfState, NumOfState);
	C_q_qdot_t = dlib::zeros_matrix<double>(NumOfState, 1);

	// Evaluation the mass matrix
	D_q_t[0][0] = 8.7E-3;
	D_q_t[0][3] = sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6808047E-5+cos(alpha_t)*cos(gamma_t)*1.4330727E-5+cos(beta_t)*sin(alpha_t)*5.45726814E-4+cos(alpha_t)*sin(gamma_t)*2.7743169E-4-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.7743169E-4+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5;
	D_q_t[0][4] = cos(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*1.131907068074867E-20;
	D_q_t[0][5] = sin(alpha_t)*sin(gamma_t)*(-1.4330727E-5)+cos(gamma_t)*sin(alpha_t)*2.7743169E-4-cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.4330727E-5-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*2.7743169E-4;
	D_q_t[1][1] = 8.7E-3;
	D_q_t[1][3] = sin(alpha_t)*sin(gamma_t)*2.7743169E-4-cos(alpha_t)*cos(beta_t)*5.45726814E-4-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6808047E-5+cos(gamma_t)*sin(alpha_t)*1.4330727E-5+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.7743169E-4-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5;
	D_q_t[1][4] = sin(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*1.131907068074867E-20;
	D_q_t[1][5] = cos(alpha_t)*cos(gamma_t)*(-2.7743169E-4)+cos(alpha_t)*sin(gamma_t)*1.4330727E-5-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.4330727E-5-sin(alpha_t)*sin(beta_t)*sin(gamma_t)*2.7743169E-4;
	D_q_t[2][2] = 8.7E-3;
	D_q_t[2][4] = cos(beta_t)*5.45726814E-4+cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6808047E-5-cos(gamma_t)*sin(beta_t)*2.7743169E-4+sin(beta_t)*sin(gamma_t)*1.4330727E-5;
	D_q_t[2][5] = cos(beta_t)*(cos(gamma_t)*1.266069221068963E15+sin(gamma_t)*2.451011198930424E16)*(-1.131907068074867E-20);
	D_q_t[3][0] = sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6808047E-5+cos(alpha_t)*cos(gamma_t)*1.4330727E-5+cos(beta_t)*sin(alpha_t)*5.45726814E-4+cos(alpha_t)*sin(gamma_t)*2.7743169E-4-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.7743169E-4+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5;
	D_q_t[3][1] = sin(alpha_t)*sin(gamma_t)*2.7743169E-4-cos(alpha_t)*cos(beta_t)*5.45726814E-4-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6808047E-5+cos(gamma_t)*sin(alpha_t)*1.4330727E-5+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.7743169E-4-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5;
	D_q_t[3][3] = pow(sin(alpha_t)*sin(gamma_t)*2.451011198930424E16-cos(alpha_t)*cos(beta_t)*4.821304057480314E16-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15+cos(gamma_t)*sin(alpha_t)*1.266069221068963E15+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.451011198930424E16-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.266069221068963E15,2.0)*1.472659322710162E-38+pow(sin(beta_t+3.141592653589793*(1.0/1.8E1)),2.0)*3.6E-6+pow(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15+cos(alpha_t)*cos(gamma_t)*1.266069221068963E15+cos(beta_t)*sin(alpha_t)*4.821304057480314E16+cos(alpha_t)*sin(gamma_t)*2.451011198930424E16-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.451011198930424E16+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.266069221068963E15,2.0)*1.472659322710162E-38+pow(cos(beta_t+3.141592653589793*(1.0/1.8E1)),2.0)*pow(cos(gamma_t-q1_t),2.0)*3.6E-6+pow(sin(gamma_t-q1_t),2.0)*pow(cos(beta_t+3.141592653589793*(1.0/1.8E1)),2.0)*7.698E-7;
	D_q_t[3][4] = sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*(-2.8302E-6)+sin(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*(sin(alpha_t)*sin(gamma_t)*3.18887E-2-cos(alpha_t)*cos(beta_t)*6.272722E-2-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*4.23081E-3+cos(gamma_t)*sin(alpha_t)*1.64721E-3+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*3.18887E-2-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*1.131907068074867E-20+cos(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*4.23081E-3+cos(alpha_t)*cos(gamma_t)*1.64721E-3+cos(beta_t)*sin(alpha_t)*6.272722E-2+cos(alpha_t)*sin(gamma_t)*3.18887E-2-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*3.18887E-2+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*1.131907068074867E-20;
	D_q_t[3][5] = sin(beta_t)*(-8.870541649724669E-6)-cos(beta_t)*sin(gamma_t)*8.9892666528894E-7-cos(3.141592653589793*(1.0/1.8E1))*sin(beta_t)*3.6E-6-sin(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*3.6E-6+cos(beta_t)*cos(gamma_t)*1.74025186536018E-5+cos(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*cos(gamma_t)*1.1737607683689E-6-cos(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*sin(gamma_t)*6.063058309887E-8-sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*sin(beta_t)*1.1737607683689E-6+sin(3.141592653589793*(1.0/1.8E1))*sin(beta_t)*sin(gamma_t)*6.063058309887E-8;
	D_q_t[3][6] = sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.6E-6;
	D_q_t[4][0] = cos(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*1.131907068074867E-20;
	D_q_t[4][1] = sin(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*1.131907068074867E-20;
	D_q_t[4][2] = cos(beta_t)*5.45726814E-4+cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6808047E-5-cos(gamma_t)*sin(beta_t)*2.7743169E-4+sin(beta_t)*sin(gamma_t)*1.4330727E-5;
	D_q_t[4][3] = sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*(-2.8302E-6)+sin(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*(sin(alpha_t)*sin(gamma_t)*3.18887E-2-cos(alpha_t)*cos(beta_t)*6.272722E-2-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*4.23081E-3+cos(gamma_t)*sin(alpha_t)*1.64721E-3+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*3.18887E-2-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*1.131907068074867E-20+cos(alpha_t)*(sin(beta_t)*4.821304057480314E16+sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.251861220603796E15-cos(beta_t)*sin(gamma_t)*1.266069221068963E15+cos(beta_t)*cos(gamma_t)*2.451011198930424E16)*(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*4.23081E-3+cos(alpha_t)*cos(gamma_t)*1.64721E-3+cos(beta_t)*sin(alpha_t)*6.272722E-2+cos(alpha_t)*sin(gamma_t)*3.18887E-2-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*3.18887E-2+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*1.131907068074867E-20;
	D_q_t[4][4] = cos(3.141592653589793*(1.0/1.8E1))*4.61773292387868E-6-sin(gamma_t*2.0)*4.569882540849E-7+pow(cos(gamma_t),2.0)*1.165353021608133E-5+pow(cos(q1_t),2.0)*2.8302E-6+sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*2.3475215367378E-6-sin(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.2126116619774E-7-pow(cos(gamma_t),2.0)*pow(cos(q1_t),2.0)*5.6604E-6-cos(gamma_t)*cos(q1_t)*sin(gamma_t)*sin(q1_t)*5.6604E-6+3.518105949182682E-5;
	D_q_t[4][5] = (cos(3.141592653589793*(1.0/1.8E1))*8.12965305150949E14+1.205326014370078E16)*(cos(gamma_t)*1.266069221068963E15+sin(gamma_t)*2.451011198930424E16)*(-5.890637290840647E-38);
	D_q_t[5][0] = sin(alpha_t)*sin(gamma_t)*(-1.4330727E-5)+cos(gamma_t)*sin(alpha_t)*2.7743169E-4-cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.4330727E-5-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*2.7743169E-4;
	D_q_t[5][1] = cos(alpha_t)*cos(gamma_t)*(-2.7743169E-4)+cos(alpha_t)*sin(gamma_t)*1.4330727E-5-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.4330727E-5-sin(alpha_t)*sin(beta_t)*sin(gamma_t)*2.7743169E-4;
	D_q_t[5][2] = cos(beta_t)*(cos(gamma_t)*1.266069221068963E15+sin(gamma_t)*2.451011198930424E16)*(-1.131907068074867E-20);
	D_q_t[5][3] = sin(beta_t)*(-8.870541649724669E-6)-cos(beta_t)*sin(gamma_t)*8.9892666528894E-7-cos(3.141592653589793*(1.0/1.8E1))*sin(beta_t)*3.6E-6-sin(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*3.6E-6+cos(beta_t)*cos(gamma_t)*1.74025186536018E-5+cos(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*cos(gamma_t)*1.1737607683689E-6-cos(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*sin(gamma_t)*6.063058309887E-8-sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*sin(beta_t)*1.1737607683689E-6+sin(3.141592653589793*(1.0/1.8E1))*sin(beta_t)*sin(gamma_t)*6.063058309887E-8;
	D_q_t[5][4] = (cos(3.141592653589793*(1.0/1.8E1))*8.12965305150949E14+1.205326014370078E16)*(cos(gamma_t)*1.266069221068963E15+sin(gamma_t)*2.451011198930424E16)*(-5.890637290840647E-38);
	D_q_t[5][5] = pow(cos(alpha_t)*cos(gamma_t)*2.451011198930424E16-cos(alpha_t)*sin(gamma_t)*1.266069221068963E15+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.266069221068963E15+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*2.451011198930424E16,2.0)*1.472659322710162E-38+pow(sin(alpha_t)*sin(gamma_t)*1.266069221068963E15-cos(gamma_t)*sin(alpha_t)*2.451011198930424E16+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.266069221068963E15+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*2.451011198930424E16,2.0)*1.472659322710162E-38+pow(cos(beta_t),2.0)*pow(cos(gamma_t)*1.266069221068963E15+sin(gamma_t)*2.451011198930424E16,2.0)*1.472659322710162E-38+3.6E-6;
	D_q_t[5][6] = -3.6E-6;
	D_q_t[6][3] = sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.6E-6;
	D_q_t[6][5] = -3.6E-6;
	D_q_t[6][6] = 3.6E-6;

	// Evaluating the Coriolis matrix

	C_q_qdot_t[0][0] = (alphadot_t*alphadot_t)*cos(alpha_t)*cos(beta_t)*5.45726814E-4+(betadot_t*betadot_t)*cos(alpha_t)*cos(beta_t)*5.45726814E-4-(alphadot_t*alphadot_t)*cos(gamma_t)*sin(alpha_t)*1.4330727E-5-(gammadot_t*gammadot_t)*cos(gamma_t)*sin(alpha_t)*1.4330727E-5-(alphadot_t*alphadot_t)*sin(alpha_t)*sin(gamma_t)*2.7743169E-4-(gammadot_t*gammadot_t)*sin(alpha_t)*sin(gamma_t)*2.7743169E-4-(alphadot_t*alphadot_t)*cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.7743169E-4-(betadot_t*betadot_t)*cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.7743169E-4-(gammadot_t*gammadot_t)*cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.7743169E-4+(alphadot_t*alphadot_t)*cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5+(betadot_t*betadot_t)*cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5+(gammadot_t*gammadot_t)*cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5+alphadot_t*gammadot_t*cos(alpha_t)*cos(gamma_t)*5.5486338E-4+(alphadot_t*alphadot_t)*cos(3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(beta_t)*3.6808047E-5+(betadot_t*betadot_t)*cos(3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(beta_t)*3.6808047E-5-alphadot_t*gammadot_t*cos(alpha_t)*sin(gamma_t)*2.8661454E-5-alphadot_t*betadot_t*sin(alpha_t)*sin(beta_t)*1.091453628E-3-(alphadot_t*alphadot_t)*sin(3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*sin(beta_t)*3.6808047E-5-(betadot_t*betadot_t)*sin(3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*sin(beta_t)*3.6808047E-5-betadot_t*gammadot_t*cos(alpha_t)*cos(beta_t)*cos(gamma_t)*2.8661454E-5-alphadot_t*betadot_t*cos(beta_t)*cos(gamma_t)*sin(alpha_t)*5.5486338E-4-betadot_t*gammadot_t*cos(alpha_t)*cos(beta_t)*sin(gamma_t)*5.5486338E-4+alphadot_t*betadot_t*cos(beta_t)*sin(alpha_t)*sin(gamma_t)*2.8661454E-5+alphadot_t*gammadot_t*cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.8661454E-5+alphadot_t*gammadot_t*sin(alpha_t)*sin(beta_t)*sin(gamma_t)*5.5486338E-4-alphadot_t*betadot_t*cos(3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*sin(beta_t)*7.3616094E-5-alphadot_t*betadot_t*sin(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*sin(alpha_t)*7.3616094E-5;
	C_q_qdot_t[1][0] = (alphadot_t*alphadot_t)*cos(alpha_t)*cos(gamma_t)*1.4330727E-5+(gammadot_t*gammadot_t)*cos(alpha_t)*cos(gamma_t)*1.4330727E-5+(alphadot_t*alphadot_t)*cos(beta_t)*sin(alpha_t)*5.45726814E-4+(betadot_t*betadot_t)*cos(beta_t)*sin(alpha_t)*5.45726814E-4+(alphadot_t*alphadot_t)*cos(alpha_t)*sin(gamma_t)*2.7743169E-4+(gammadot_t*gammadot_t)*cos(alpha_t)*sin(gamma_t)*2.7743169E-4-(alphadot_t*alphadot_t)*cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.7743169E-4-(betadot_t*betadot_t)*cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.7743169E-4-(gammadot_t*gammadot_t)*cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.7743169E-4+(alphadot_t*alphadot_t)*sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5+(betadot_t*betadot_t)*sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5+(gammadot_t*gammadot_t)*sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.4330727E-5+alphadot_t*betadot_t*cos(alpha_t)*sin(beta_t)*1.091453628E-3+alphadot_t*gammadot_t*cos(gamma_t)*sin(alpha_t)*5.5486338E-4+(alphadot_t*alphadot_t)*cos(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*sin(alpha_t)*3.6808047E-5+(betadot_t*betadot_t)*cos(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*sin(alpha_t)*3.6808047E-5-alphadot_t*gammadot_t*sin(alpha_t)*sin(gamma_t)*2.8661454E-5-(alphadot_t*alphadot_t)*sin(3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*sin(beta_t)*3.6808047E-5-(betadot_t*betadot_t)*sin(3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*sin(beta_t)*3.6808047E-5+alphadot_t*betadot_t*cos(alpha_t)*cos(beta_t)*cos(gamma_t)*5.5486338E-4-alphadot_t*betadot_t*cos(alpha_t)*cos(beta_t)*sin(gamma_t)*2.8661454E-5-alphadot_t*gammadot_t*cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.8661454E-5-betadot_t*gammadot_t*cos(beta_t)*cos(gamma_t)*sin(alpha_t)*2.8661454E-5-alphadot_t*gammadot_t*cos(alpha_t)*sin(beta_t)*sin(gamma_t)*5.5486338E-4-betadot_t*gammadot_t*cos(beta_t)*sin(alpha_t)*sin(gamma_t)*5.5486338E-4+alphadot_t*betadot_t*cos(3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*sin(beta_t)*7.3616094E-5+alphadot_t*betadot_t*sin(3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(beta_t)*7.3616094E-5;
	C_q_qdot_t[2][0] = (betadot_t*betadot_t)*sin(beta_t)*(-5.45726814E-4)-(betadot_t*betadot_t)*cos(3.141592653589793*(1.0/1.8E1))*sin(beta_t)*3.6808047E-5-(betadot_t*betadot_t)*sin(3.141592653589793*(1.0/1.8E1))*cos(beta_t)*3.6808047E-5-(betadot_t*betadot_t)*cos(beta_t)*cos(gamma_t)*2.7743169E-4-(gammadot_t*gammadot_t)*cos(beta_t)*cos(gamma_t)*2.7743169E-4+(betadot_t*betadot_t)*cos(beta_t)*sin(gamma_t)*1.4330727E-5+(gammadot_t*gammadot_t)*cos(beta_t)*sin(gamma_t)*1.4330727E-5+betadot_t*gammadot_t*cos(gamma_t)*sin(beta_t)*2.8661454E-5+betadot_t*gammadot_t*sin(beta_t)*sin(gamma_t)*5.5486338E-4+8.534699999999999E-2;
	C_q_qdot_t[3][0] = gammaddot*sin(beta_t)*8.870541649724669E-6-gammaddot*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.6E-6+betadot_t*q1dot_t*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6E-6-gammaddot*cos(beta_t)*cos(gamma_t)*1.74025186536018E-5+gammaddot*cos(beta_t)*sin(gamma_t)*8.9892666528894E-7+alphadot_t*betadot_t*sin(beta_t*2.0+3.141592653589793*(1.0/9.0))*3.6E-6+gammaddot*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6-gammaddot*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.8E-6*sqrt(-1.0)+gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6*sqrt(-1.0)-gammaddot*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)+gammaddot*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.8E-6+gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6-betadot_t*gammadot_t*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6E-6+(betadot_t*betadot_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.1737607683689E-6-(gammadot_t*gammadot_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.1737607683689E-6-gammaddot*pow(sin(alpha_t),2.0)*sin(beta_t)*pow(sin(gamma_t),2.0)*8.870541649724669E-6+(betadot_t*betadot_t)*sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.8302E-6-(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*pow(cos(gamma_t),2.0)*sin(beta_t)*4.569882540849E-7+(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*sin(beta_t)*pow(sin(gamma_t),2.0)*4.569882540849E-7-(betadot_t*betadot_t)*pow(cos(gamma_t),2.0)*pow(sin(alpha_t),2.0)*sin(beta_t)*4.569882540849E-7-betadot_t*gammadot_t*cos(beta_t+3.141592653589793*(1.0/1.8E1))*pow(cos(gamma_t-q1_t),2.0)*2.8302E-6+betadot_t*q1dot_t*cos(beta_t+3.141592653589793*(1.0/1.8E1))*pow(cos(gamma_t-q1_t),2.0)*2.8302E-6+(betadot_t*betadot_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*pow(sin(gamma_t),2.0)*4.569882540849E-7+betadot_t*gammadot_t*pow(sin(gamma_t-q1_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.8302E-6+gammaddot*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*5.8688038418445E-7*sqrt(-1.0)-gammaddot*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*5.8688038418445E-7-gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*5.8688038418445E-7-betadot_t*q1dot_t*pow(sin(gamma_t-q1_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.8302E-6+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*pow(cos(gamma_t),2.0)*9.139765081698E-7+gammaddot*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*1.74025186536018E-5-gammaddot*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)+gammaddot*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*5.8688038418445E-7+gammaddot*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*3.0315291549435E-8-gammaddot*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*5.8688038418445E-7*sqrt(-1.0)+gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*3.0315291549435E-8+gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*5.8688038418445E-7*sqrt(-1.0)-alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*pow(sin(gamma_t),2.0)*9.139765081698E-7+alphadot_t*gammadot_t*pow(cos(gamma_t),2.0)*pow(sin(alpha_t),2.0)*9.139765081698E-7+gammaddot*pow(cos(alpha_t),2.0)*cos(gamma_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.1737607683689E-6+gammaddot*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*1.74025186536018E-5-gammaddot*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(gamma_t)*8.9892666528894E-7-gammaddot*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*3.0315291549435E-8+gammaddot*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-alphadot_t*gammadot_t*pow(sin(alpha_t),2.0)*pow(sin(gamma_t),2.0)*9.139765081698E-7+gammaddot*cos(gamma_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.1737607683689E-6-gammaddot*pow(cos(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*6.063058309887E-8-gammaddot*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(gamma_t)*8.9892666528894E-7-gammaddot*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*6.063058309887E-8+(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*8.9892666528894E-7-(gammadot_t*gammadot_t)*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*8.9892666528894E-7+(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*cos(gamma_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*6.063058309887E-8+(betadot_t*betadot_t)*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*8.9892666528894E-7+(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(gamma_t)*1.74025186536018E-5-(gammadot_t*gammadot_t)*pow(cos(alpha_t),2.0)*cos(gamma_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*6.063058309887E-8-gammaddot*pow(cos(alpha_t),2.0)*pow(cos(gamma_t),2.0)*sin(beta_t)*8.870541649724669E-6-(gammadot_t*gammadot_t)*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*8.9892666528894E-7-(gammadot_t*gammadot_t)*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(gamma_t)*1.74025186536018E-5+(betadot_t*betadot_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*6.063058309887E-8+(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.1737607683689E-6+(betadot_t*betadot_t)*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(gamma_t)*1.74025186536018E-5-(gammadot_t*gammadot_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*6.063058309887E-8-(gammadot_t*gammadot_t)*pow(cos(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.1737607683689E-6-gammaddot*pow(cos(alpha_t),2.0)*sin(beta_t)*pow(sin(gamma_t),2.0)*8.870541649724669E-6-gammaddot*pow(cos(gamma_t),2.0)*pow(sin(alpha_t),2.0)*sin(beta_t)*8.870541649724669E-6-(gammadot_t*gammadot_t)*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(gamma_t)*1.74025186536018E-5-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.1145570665614E-7+alphadot_t*gammadot_t*cos(gamma_t)*pow(sin(alpha_t),2.0)*sin(gamma_t)*1.764666043216266E-5-alphadot_t*betadot_t*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*4.61773292387868E-6-alphadot_t*betadot_t*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t)*4.61773292387868E-6-alphadot_t*betadot_t*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.1145570665614E-7-alphadot_t*betadot_t*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t+3.141592653589793*(1.0/1.8E1))*pow(cos(gamma_t-q1_t),2.0)*7.2E-6-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*pow(cos(beta_t),2.0)*cos(gamma_t)*3.48050373072036E-5-betadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*pow(cos(gamma_t),2.0)*4.721143364334E-8-alphadot_t*betadot_t*pow(sin(gamma_t-q1_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t+3.141592653589793*(1.0/1.8E1))*1.5396E-6+alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(gamma_t)*pow(sin(beta_t),2.0)*3.48050373072036E-5+alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*pow(cos(beta_t),2.0)*sin(gamma_t)*1.79785333057788E-6-alphadot_t*betadot_t*pow(cos(beta_t),2.0)*cos(gamma_t)*pow(sin(alpha_t),2.0)*3.48050373072036E-5-betadot_t*gammadot_t*cos(beta_t)*pow(cos(gamma_t),2.0)*pow(sin(alpha_t),2.0)*4.721143364334E-8-betadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*pow(sin(gamma_t),2.0)*1.7693871865806E-5+alphadot_t*betadot_t*cos(gamma_t)*pow(sin(alpha_t),2.0)*pow(sin(beta_t),2.0)*3.48050373072036E-5-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*pow(sin(beta_t),2.0)*sin(gamma_t)*1.79785333057788E-6+alphadot_t*betadot_t*pow(cos(beta_t),2.0)*pow(sin(alpha_t),2.0)*sin(gamma_t)*1.79785333057788E-6-betadot_t*gammadot_t*cos(beta_t)*pow(sin(alpha_t),2.0)*pow(sin(gamma_t),2.0)*1.7693871865806E-5-alphadot_t*gammadot_t*sin(gamma_t-q1_t)*pow(cos(beta_t+3.141592653589793*(1.0/1.8E1)),2.0)*cos(gamma_t-q1_t)*5.6604E-6-alphadot_t*betadot_t*pow(sin(alpha_t),2.0)*pow(sin(beta_t),2.0)*sin(gamma_t)*1.79785333057788E-6+alphadot_t*q1dot_t*sin(gamma_t-q1_t)*pow(cos(beta_t+3.141592653589793*(1.0/1.8E1)),2.0)*cos(gamma_t-q1_t)*5.6604E-6-alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*pow(cos(gamma_t),2.0)*pow(sin(beta_t),2.0)*9.139765081698E-7+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*pow(sin(beta_t),2.0)*pow(sin(gamma_t),2.0)*9.139765081698E-7-alphadot_t*gammadot_t*pow(cos(gamma_t),2.0)*pow(sin(alpha_t),2.0)*pow(sin(beta_t),2.0)*9.139765081698E-7-(betadot_t*betadot_t)*pow(cos(alpha_t),2.0)*cos(gamma_t)*sin(beta_t)*sin(gamma_t)*8.823330216081329E-6+alphadot_t*gammadot_t*pow(sin(alpha_t),2.0)*pow(sin(beta_t),2.0)*pow(sin(gamma_t),2.0)*9.139765081698E-7-(betadot_t*betadot_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*sin(gamma_t)*8.823330216081329E-6-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(beta_t)*6.846385184335416E-5+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(gamma_t)*sin(gamma_t)*1.764666043216266E-5-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*4.61773292387868E-6-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t)*4.61773292387868E-6-alphadot_t*betadot_t*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*6.846385184335416E-5+alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*pow(cos(gamma_t),2.0)*sin(beta_t)*1.7693871865806E-5+alphadot_t*betadot_t*cos(beta_t)*pow(cos(gamma_t),2.0)*pow(sin(alpha_t),2.0)*sin(beta_t)*1.7693871865806E-5+alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(beta_t)*pow(sin(gamma_t),2.0)*4.721143364334E-8-alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(gamma_t)*pow(sin(beta_t),2.0)*sin(gamma_t)*1.764666043216266E-5+alphadot_t*betadot_t*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*pow(sin(gamma_t),2.0)*4.721143364334E-8-alphadot_t*gammadot_t*cos(gamma_t)*pow(sin(alpha_t),2.0)*pow(sin(beta_t),2.0)*sin(gamma_t)*1.764666043216266E-5-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.3475215367378E-6+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*sin(beta_t)*1.79785333057788E-6-betadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*sin(gamma_t)*1.8279530163396E-6-alphadot_t*betadot_t*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.3475215367378E-6+alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.2126116619774E-7+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(gamma_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t)*1.2126116619774E-7+alphadot_t*gammadot_t*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*1.79785333057788E-6+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*sin(beta_t)*sin(gamma_t)*3.48050373072036E-5-betadot_t*gammadot_t*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*sin(gamma_t)*1.8279530163396E-6+alphadot_t*betadot_t*cos(beta_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*1.2126116619774E-7+alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(gamma_t)*sin(beta_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.3475215367378E-6+alphadot_t*gammadot_t*cos(gamma_t)*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t)*1.2126116619774E-7+alphadot_t*gammadot_t*pow(cos(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t)*sin(gamma_t)*2.3475215367378E-6+alphadot_t*gammadot_t*cos(beta_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*sin(gamma_t)*3.48050373072036E-5+alphadot_t*betadot_t*cos(gamma_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.3475215367378E-6-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*sin(beta_t)*sin(gamma_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*1.2126116619774E-7+alphadot_t*gammadot_t*pow(sin(alpha_t),2.0)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*sin(beta_t)*sin(gamma_t)*2.3475215367378E-6-alphadot_t*betadot_t*pow(sin(alpha_t),2.0)*sin(beta_t)*sin(gamma_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*1.2126116619774E-7-alphadot_t*betadot_t*pow(cos(alpha_t),2.0)*cos(beta_t)*cos(gamma_t)*sin(beta_t)*sin(gamma_t)*1.8279530163396E-6-alphadot_t*betadot_t*cos(beta_t)*cos(gamma_t)*pow(sin(alpha_t),2.0)*sin(beta_t)*sin(gamma_t)*1.8279530163396E-6;
	C_q_qdot_t[4][0] = cos(beta_t)*5.35358004534E-3-cos(gamma_t)*sin(beta_t)*2.7216048789E-3+(gammadot_t*gammadot_t)*sin(gamma_t)*8.9892666528894E-7+sin(beta_t)*sin(gamma_t)*1.4058443187E-4+exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.80543470535E-4-exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.80543470535E-4+exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.80543470535E-4*sqrt(-1.0)+exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.80543470535E-4*sqrt(-1.0)+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*1.489832754840737E-5-betaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*2.30886646193934E-6+betaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*2.30886646193934E-6-betaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*2.30886646193934E-6-(gammadot_t*gammadot_t)*cos(gamma_t)*1.74025186536018E-5-(gammadot_t*gammadot_t)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*5.8688038418445E-7+(gammadot_t*gammadot_t)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*5.8688038418445E-7+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(gamma_t)*1.74025186536018E-5+(gammadot_t*gammadot_t)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*3.0315291549435E-8-(gammadot_t*gammadot_t)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*3.0315291549435E-8-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*sin(gamma_t)*8.9892666528894E-7-betadot_t*gammadot_t*cos(gamma_t*2.0)*9.139765081698E-7+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.148430366679825E-7*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.148430366679825E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*1.15443323096967E-6*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*1.15443323096967E-6*sqrt(-1.0)-betadot_t*gammadot_t*sin(gamma_t*2.0)*8.823330216081329E-6-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.148430366679825E-7+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.148430366679825E-7+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*1.15443323096967E-6-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*1.15443323096967E-6-(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*sin(beta_t*2.0)*2.205832554020332E-6+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*sin(gamma_t*2.0)*2.2849412704245E-7-betaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)+betaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)-betaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)+gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8-gammaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8+gammaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8+betaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)-betaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)+betaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)+gammaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7-gammaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7+gammaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7+alphadot_t*gammadot_t*cos(beta_t)*8.870541649724669E-6+alphadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.8E-6-alphadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6+betadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*6.063058309887E-8*sqrt(-1.0)+betadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*6.063058309887E-8*sqrt(-1.0)-alphadot_t*q1dot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.8E-6+alphadot_t*q1dot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6+alphadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)+alphadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)+alphadot_t*gammadot_t*cos(gamma_t*2.0)*cos(beta_t)*8.823330216081329E-6+betadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*1.1737607683689E-6*sqrt(-1.0)+betadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*1.1737607683689E-6*sqrt(-1.0)-alphadot_t*q1dot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)-alphadot_t*q1dot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*5.8688038418445E-7-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*5.8688038418445E-7-alphadot_t*gammadot_t*sin(gamma_t*2.0)*cos(beta_t)*9.139765081698E-7-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*3.0315291549435E-8+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*5.8688038418445E-7*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*3.0315291549435E-8+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*5.8688038418445E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-betadot_t*gammadot_t*cos(gamma_t*2.0)*sin(q1_t*2.0)*2.8302E-6+betadot_t*gammadot_t*cos(q1_t*2.0)*sin(gamma_t*2.0)*2.8302E-6+betadot_t*q1dot_t*cos(gamma_t*2.0)*sin(q1_t*2.0)*2.8302E-6-betadot_t*q1dot_t*cos(q1_t*2.0)*sin(gamma_t*2.0)*2.8302E-6+alphadot_t*gammadot_t*cos(gamma_t)*sin(beta_t)*3.48050373072036E-5-alphadot_t*gammadot_t*sin(beta_t)*sin(gamma_t)*1.79785333057788E-6-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*cos(q1_t*2.0)*sin(beta_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7-(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*cos(q1_t*2.0)*sin(beta_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7-alphadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)-alphadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)+alphadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)+alphadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*1.1737607683689E-6+alphadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)-alphadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*1.1737607683689E-6+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7-alphadot_t*gammadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*6.063058309887E-8+alphadot_t*gammadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*6.063058309887E-8-alphadot_t*gammadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.4151E-6+alphadot_t*gammadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.4151E-6+alphadot_t*q1dot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.4151E-6-alphadot_t*q1dot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.4151E-6-alphadot_t*gammadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)-alphadot_t*gammadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*q1dot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*q1dot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)-alphadot_t*gammadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.4151E-6+alphadot_t*gammadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.4151E-6+alphadot_t*q1dot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.4151E-6-alphadot_t*q1dot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.4151E-6-alphadot_t*gammadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)-alphadot_t*gammadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*q1dot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*q1dot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0);
	C_q_qdot_t[5][0] = cos(beta_t)*sin(gamma_t)*(-2.7216048789E-3)-(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*2.2849412704245E-7+(betadot_t*betadot_t)*cos(gamma_t*2.0)*4.569882540849E-7-(alphadot_t*alphadot_t)*sin(gamma_t*2.0)*2.205832554020332E-6+(betadot_t*betadot_t)*sin(gamma_t*2.0)*4.411665108040665E-6-cos(beta_t)*cos(gamma_t)*1.4058443187E-4-(alphadot_t*alphadot_t)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*1.51576457747175E-8*sqrt(-1.0)-(alphadot_t*alphadot_t)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*1.51576457747175E-8*sqrt(-1.0)-(betadot_t*betadot_t)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-(betadot_t*betadot_t)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-(alphadot_t*alphadot_t)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*2.93440192092225E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*2.93440192092225E-7*sqrt(-1.0)-(betadot_t*betadot_t)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7*sqrt(-1.0)-(betadot_t*betadot_t)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*cos(gamma_t)*4.4946333264447E-7-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*sin(gamma_t)*8.7012593268009E-6-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(gamma_t*2.0)*2.2849412704245E-7-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*sin(gamma_t*2.0)*2.205832554020332E-6-(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*sin(q1_t*2.0)*7.0755E-7+(alphadot_t*alphadot_t)*cos(q1_t*2.0)*sin(gamma_t*2.0)*7.0755E-7+(betadot_t*betadot_t)*cos(gamma_t*2.0)*sin(q1_t*2.0)*1.4151E-6-(betadot_t*betadot_t)*cos(q1_t*2.0)*sin(gamma_t*2.0)*1.4151E-6+alphaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6*sqrt(-1.0)-alphaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6*sqrt(-1.0)+alphaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6*sqrt(-1.0)+betaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8-betaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8+betaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*3.0315291549435E-8+alphaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6-alphaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6+alphaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6+betaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7-betaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7+betaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*5.8688038418445E-7-alphadot_t*betadot_t*cos(beta_t)*8.870541649724669E-6-alphadot_t*betadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)-alphadot_t*betadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.8E-6*sqrt(-1.0)-alphadot_t*betadot_t*cos(gamma_t*2.0)*cos(beta_t)*8.823330216081329E-6+alphadot_t*betadot_t*sin(gamma_t*2.0)*cos(beta_t)*9.139765081698E-7+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*1.51576457747175E-8*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*1.51576457747175E-8*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*2.93440192092225E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*1.51576457747175E-8+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*2.93440192092225E-7*sqrt(-1.0)+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*1.51576457747175E-8-(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(gamma_t)*2.93440192092225E-7+(alphadot_t*alphadot_t)*sin(beta_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(gamma_t)*2.93440192092225E-7-alphaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*5.8688038418445E-7+alphaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*5.8688038418445E-7-alphaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*5.8688038418445E-7+alphaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*3.0315291549435E-8+alphaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*5.8688038418445E-7*sqrt(-1.0)-alphaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*3.0315291549435E-8-alphaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*5.8688038418445E-7*sqrt(-1.0)+alphaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*3.0315291549435E-8+alphaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*5.8688038418445E-7*sqrt(-1.0)-alphaddot*exp(3.141592653589793*2.777777777777778E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)+alphaddot*exp(3.141592653589793*6.111111111111111E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-alphaddot*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*3.0315291549435E-8*sqrt(-1.0)-alphadot_t*betadot_t*cos(gamma_t)*sin(beta_t)*3.48050373072036E-5+alphadot_t*betadot_t*sin(beta_t)*sin(gamma_t)*1.79785333057788E-6-alphadot_t*betadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.8E-6+alphadot_t*betadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.8E-6-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(q1_t*2.0)*sin(gamma_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7+(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7-(alphadot_t*alphadot_t)*cos(beta_t*2.0)*cos(q1_t*2.0)*sin(gamma_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7+alphadot_t*betadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)+alphadot_t*betadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*cos(gamma_t)*1.1737607683689E-6*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*sin(beta_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(q1_t*2.0)*sin(beta_t*2.0)*sin(gamma_t*2.0)*exp(3.141592653589793*1.111111111111111E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)-(alphadot_t*alphadot_t)*cos(gamma_t*2.0)*sin(beta_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)+(alphadot_t*alphadot_t)*cos(q1_t*2.0)*sin(beta_t*2.0)*sin(gamma_t*2.0)*exp(3.141592653589793*8.888888888888889E-1*sqrt(-1.0))*3.53775E-7*sqrt(-1.0)-alphadot_t*betadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)-alphadot_t*betadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*1.1737607683689E-6-alphadot_t*betadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*sin(gamma_t)*6.063058309887E-8*sqrt(-1.0)+alphadot_t*betadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(gamma_t)*sin(beta_t)*1.1737607683689E-6+alphadot_t*betadot_t*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*6.063058309887E-8-alphadot_t*betadot_t*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*sin(gamma_t)*6.063058309887E-8+alphadot_t*betadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.4151E-6-alphadot_t*betadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.4151E-6+alphadot_t*betadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*betadot_t*cos(gamma_t*2.0)*cos(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*betadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*cos(beta_t)*1.4151E-6-alphadot_t*betadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*cos(beta_t)*1.4151E-6+alphadot_t*betadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*5.555555555555556E-2*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0)+alphadot_t*betadot_t*sin(gamma_t*2.0)*sin(q1_t*2.0)*exp(3.141592653589793*9.444444444444444E-1*sqrt(-1.0))*sin(beta_t)*1.4151E-6*sqrt(-1.0);
	C_q_qdot_t[6][0] = (alphadot_t*alphadot_t)*sin(gamma_t*2.0-q1_t*2.0)*(-7.0755E-7)+(betadot_t*betadot_t)*sin(gamma_t*2.0-q1_t*2.0)*1.4151E-6+(alphadot_t*alphadot_t)*sin(beta_t*2.0-gamma_t*2.0+q1_t*2.0+3.141592653589793*(1.0/9.0))*3.53775E-7-(alphadot_t*alphadot_t)*sin(beta_t*2.0+gamma_t*2.0-q1_t*2.0+3.141592653589793*(1.0/9.0))*3.53775E-7-alphadot_t*betadot_t*cos(beta_t-gamma_t*2.0+q1_t*2.0+3.141592653589793*(1.0/1.8E1))*1.4151E-6-alphadot_t*betadot_t*cos(beta_t+gamma_t*2.0-q1_t*2.0+3.141592653589793*(1.0/1.8E1))*1.4151E-6+alphadot_t*betadot_t*cos(beta_t+3.141592653589793*(1.0/1.8E1))*3.6E-6;

}

void D_q_C_q_qdot_LeftArm(dlib::matrix<double> D_q, dlib::matrix<double> C_q_qdot)
{
	D_q_t = dlib::zeros_matrix<double>(NumOfState, NumOfState);
	C_q_qdot_t = dlib::zeros_matrix<double>(NumOfState, 1);

	// Evaluation the mass matrix

	D_q_t[0][0] = 2.7E1/2.0E2;
	D_q_t[0][3] = sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.63341395E-3+cos(alpha_t)*cos(gamma_t-q1_t)*6.301529999999981E-5+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*4.85381025E-3+cos(alpha_t)*cos(gamma_t)*2.2237335E-4+cos(beta_t)*sin(alpha_t)*8.4681747E-3+cos(alpha_t)*sin(gamma_t)*4.3049745E-3-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.47625E-3+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*3.5392653E-3-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*4.3049745E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*2.2237335E-4;
	D_q_t[0][4] = cos(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*4.683753385137379E-22;
	D_q_t[0][5] = sin(alpha_t)*sin(gamma_t)*(-2.2237335E-4)-sin(gamma_t-q1_t)*sin(alpha_t)*6.301529999999981E-5+cos(gamma_t)*sin(alpha_t)*4.3049745E-3+cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3-cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.2237335E-4-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*4.3049745E-3;
	D_q_t[0][6] = sin(gamma_t-q1_t)*sin(alpha_t)*6.301529999999981E-5-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3;
	D_q_t[0][7] = cos(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*3.278627369596165E-18;
	D_q_t[1][1] = 2.7E1/2.0E2;
	D_q_t[1][3] = sin(alpha_t)*sin(gamma_t)*4.3049745E-3+sin(alpha_t)*cos(gamma_t-q1_t)*6.301529999999981E-5-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*4.85381025E-3-cos(alpha_t)*cos(beta_t)*8.4681747E-3-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.63341395E-3+cos(gamma_t)*sin(alpha_t)*2.2237335E-4+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.47625E-3-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*3.5392653E-3+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*4.3049745E-3-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*2.2237335E-4;
	D_q_t[1][4] = sin(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*4.683753385137379E-22;
	D_q_t[1][5] = sin(gamma_t-q1_t)*cos(alpha_t)*6.301529999999981E-5-cos(alpha_t)*cos(gamma_t)*4.3049745E-3+cos(alpha_t)*sin(gamma_t)*2.2237335E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.2237335E-4-sin(alpha_t)*sin(beta_t)*sin(gamma_t)*4.3049745E-3;
	D_q_t[1][6] = sin(gamma_t-q1_t)*cos(alpha_t)*(-6.301529999999981E-5)-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3;
	D_q_t[1][7] = sin(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*3.278627369596165E-18;
	D_q_t[2][2] = 2.7E1/2.0E2;
	D_q_t[2][4] = cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*4.85381025E-3+cos(beta_t)*8.4681747E-3+cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.63341395E-3-cos(gamma_t)*sin(beta_t)*4.3049745E-3+sin(beta_t)*sin(gamma_t)*2.2237335E-4-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.47625E-3+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.5392653E-3;
	D_q_t[2][5] = cos(beta_t)*sin(gamma_t)*(-4.3049745E-3)+cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.5392653E-3-cos(beta_t)*cos(gamma_t)*2.2237335E-4;
	D_q_t[2][6] = cos(gamma_t-q1_t)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*9.445590451108341E17-cos(beta_t+3.141592653589793*(1.0/1.8E1))*9.277415232383222E17)*3.747002708109903E-21;
	D_q_t[2][7] = cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*4.85381025E-3+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.5392653E-3;
	D_q_t[3][0] = sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.63341395E-3+cos(alpha_t)*cos(gamma_t-q1_t)*6.301529999999981E-5+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*4.85381025E-3+cos(alpha_t)*cos(gamma_t)*2.2237335E-4+cos(beta_t)*sin(alpha_t)*8.4681747E-3+cos(alpha_t)*sin(gamma_t)*4.3049745E-3-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.47625E-3+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*3.5392653E-3-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*4.3049745E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*2.2237335E-4;
	D_q_t[3][1] = sin(alpha_t)*sin(gamma_t)*4.3049745E-3+sin(alpha_t)*cos(gamma_t-q1_t)*6.301529999999981E-5-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*4.85381025E-3-cos(alpha_t)*cos(beta_t)*8.4681747E-3-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.63341395E-3+cos(gamma_t)*sin(alpha_t)*2.2237335E-4+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.47625E-3-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*3.5392653E-3+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*4.3049745E-3-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*2.2237335E-4;
	D_q_t[3][3] = pow(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18+cos(alpha_t)*cos(gamma_t-q1_t)*1.345401749800956E17+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*1.036307817871507E19+cos(alpha_t)*cos(gamma_t)*4.747759579008611E17+cos(beta_t)*sin(alpha_t)*1.807989021555118E19+cos(alpha_t)*sin(gamma_t)*9.19129199598909E18-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*7.556472360886673E18-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*9.19129199598909E18+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*4.747759579008611E17,2.0)*1.62500339057673E-42+pow(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1)),2.0)*pow(cos(gamma_t-q1_t),2.0)*1.25E-4+pow(sin(gamma_t-q1_t),2.0)*pow(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1)),2.0)*1.25E-4+pow(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1)),2.0)*3.515445E-5+pow(sin(alpha_t)*sin(gamma_t)*9.19129199598909E18+sin(alpha_t)*cos(gamma_t-q1_t)*1.345401749800956E17-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*1.036307817871507E19-cos(alpha_t)*cos(beta_t)*1.807989021555118E19-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18+cos(gamma_t)*sin(alpha_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*7.556472360886673E18+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*9.19129199598909E18-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*4.747759579008611E17,2.0)*1.62500339057673E-42;
	D_q_t[3][4] = (cos(gamma_t)*4.747759579008611E17+sin(gamma_t)*9.19129199598909E18+cos(gamma_t-q1_t)*1.345401749800956E17)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*1.62500339057673E-42;
	D_q_t[3][5] = sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*(-3.515445E-5)-(sin(alpha_t)*sin(gamma_t)*1.64721E-3+sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(gamma_t)*sin(alpha_t)*3.18887E-2-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.64721E-3+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*3.595415E-2+cos(alpha_t)*cos(gamma_t)*1.64721E-3+cos(beta_t)*sin(alpha_t)*6.272722E-2+cos(alpha_t)*sin(gamma_t)*3.18887E-2-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*2.621678E-2-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*3.18887E-2+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2)-(sin(gamma_t-q1_t)*cos(alpha_t)*(-4.667799999999986E-4)+cos(alpha_t)*cos(gamma_t)*3.18887E-2-cos(alpha_t)*sin(gamma_t)*1.64721E-3-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.64721E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(sin(alpha_t)*sin(gamma_t)*3.18887E-2+sin(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*3.595415E-2-cos(alpha_t)*cos(beta_t)*6.272722E-2-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(gamma_t)*sin(alpha_t)*1.64721E-3+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*3.18887E-2-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2);
	D_q_t[3][6] = sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.515445E-5-(sin(gamma_t-q1_t)*cos(alpha_t)*4.667799999999986E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(alpha_t)*sin(gamma_t)*3.18887E-2+sin(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*3.595415E-2-cos(alpha_t)*cos(beta_t)*6.272722E-2-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(gamma_t)*sin(alpha_t)*1.64721E-3+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*3.18887E-2-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2)+(sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*3.595415E-2+cos(alpha_t)*cos(gamma_t)*1.64721E-3+cos(beta_t)*sin(alpha_t)*6.272722E-2+cos(alpha_t)*sin(gamma_t)*3.18887E-2-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*2.621678E-2-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*3.18887E-2+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2);
	D_q_t[3][7] = (sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(cos(gamma_t)*4.747759579008611E17+sin(gamma_t)*9.19129199598909E18+cos(gamma_t-q1_t)*1.345401749800956E17)*1.137502373403711E-38;
	D_q_t[4][0] = cos(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*4.683753385137379E-22;
	D_q_t[4][1] = sin(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*4.683753385137379E-22;
	D_q_t[4][2] = cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*4.85381025E-3+cos(beta_t)*8.4681747E-3+cos(beta_t+3.141592653589793*(1.0/1.8E1))*2.63341395E-3-cos(gamma_t)*sin(beta_t)*4.3049745E-3+sin(beta_t)*sin(gamma_t)*2.2237335E-4-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*3.47625E-3+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.5392653E-3;
	D_q_t[4][3] = (cos(gamma_t)*4.747759579008611E17+sin(gamma_t)*9.19129199598909E18+cos(gamma_t-q1_t)*1.345401749800956E17)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*1.62500339057673E-42;
	D_q_t[4][4] = pow(cos(gamma_t-q1_t),2.0)*1.25E-4+pow(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+cos(beta_t)*1.807989021555118E19+cos(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(gamma_t)*sin(beta_t)*9.19129199598909E18+sin(beta_t)*sin(gamma_t)*4.747759579008611E17-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18,2.0)*1.62500339057673E-42+pow(sin(gamma_t-q1_t),2.0)*1.25E-4+pow(cos(alpha_t),2.0)*pow(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18,2.0)*1.62500339057673E-42+pow(sin(alpha_t),2.0)*pow(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18,2.0)*1.62500339057673E-42;
	D_q_t[4][5] = cos(gamma_t)*(-1.3948862047587E-5)-sin(gamma_t)*2.7003908255589E-4-sin(gamma_t)*sin(q1_t)*5.944086627349498E-5-cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*4.3377857925795E-6-cos(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*8.3976147427365E-5+sin(3.141592653589793*(1.0/1.8E1))*cos(q1_t)*1.10853093375E-4-sin(3.141592653589793*(1.0/1.8E1))*sin(q1_t)*5.7261137625E-6-cos(gamma_t)*cos(q1_t)*5.944086627349498E-5+cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q1_t)*2.18055498525E-4-cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q2_t)*7.9952447819025E-6-cos(3.141592653589793*(1.0/1.8E1))*cos(q2_t)*sin(gamma_t)*1.54781698919175E-4-cos(3.141592653589793*(1.0/1.8E1))*cos(q1_t)*sin(q2_t)*1.1286256937211E-4-sin(3.141592653589793*(1.0/1.8E1))*cos(q1_t)*cos(q2_t)*1.1286256937211E-4+cos(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*sin(q1_t)*2.18055498525E-4+sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*sin(q2_t)*7.9952447819025E-6+cos(3.141592653589793*(1.0/1.8E1))*sin(q1_t)*sin(q2_t)*5.829913194813E-6+sin(3.141592653589793*(1.0/1.8E1))*cos(q2_t)*sin(q1_t)*5.829913194813E-6+sin(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*sin(q2_t)*1.54781698919175E-4+cos(gamma_t)*cos(q1_t)*cos(q2_t)*5.5945979761419E-5+cos(q2_t)*sin(gamma_t)*sin(q1_t)*5.5945979761419E-5-cos(3.141592653589793*(1.0/1.8E1))*cos(q2_t)*sin(gamma_t)*sin(q1_t)*2.22008273111466E-4+sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q1_t)*sin(q2_t)*2.22008273111466E-4+sin(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*sin(q1_t)*sin(q2_t)*2.22008273111466E-4-cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q1_t)*cos(q2_t)*2.22008273111466E-4;
	D_q_t[4][6] = sin(alpha_t)*(sin(gamma_t-q1_t)*cos(alpha_t)*4.667799999999986E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*(-4.683753385137379E-22)+cos(gamma_t-q1_t)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*9.445590451108341E17-cos(beta_t+3.141592653589793*(1.0/1.8E1))*9.277415232383222E17)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+cos(beta_t)*6.272722E-2+cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2-cos(gamma_t)*sin(beta_t)*3.18887E-2+sin(beta_t)*sin(gamma_t)*1.64721E-3-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*3.747002708109903E-21+cos(alpha_t)*(sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*4.683753385137379E-22;
	D_q_t[4][7] = pow(cos(gamma_t-q1_t),2.0)*1.25E-4+pow(sin(gamma_t-q1_t),2.0)*1.25E-4+(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+cos(beta_t)*6.272722E-2+cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2-cos(gamma_t)*sin(beta_t)*3.18887E-2+sin(beta_t)*sin(gamma_t)*1.64721E-3-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*(2.7E1/2.0E2)+pow(cos(alpha_t),2.0)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*1.137502373403711E-38+pow(sin(alpha_t),2.0)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*1.137502373403711E-38;
	D_q_t[5][0] = sin(alpha_t)*sin(gamma_t)*(-2.2237335E-4)-sin(gamma_t-q1_t)*sin(alpha_t)*6.301529999999981E-5+cos(gamma_t)*sin(alpha_t)*4.3049745E-3+cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3-cos(alpha_t)*cos(gamma_t)*sin(beta_t)*2.2237335E-4-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*4.3049745E-3;
	D_q_t[5][1] = sin(gamma_t-q1_t)*cos(alpha_t)*6.301529999999981E-5-cos(alpha_t)*cos(gamma_t)*4.3049745E-3+cos(alpha_t)*sin(gamma_t)*2.2237335E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*2.2237335E-4-sin(alpha_t)*sin(beta_t)*sin(gamma_t)*4.3049745E-3;
	D_q_t[5][2] = cos(beta_t)*sin(gamma_t)*(-4.3049745E-3)+cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.5392653E-3-cos(beta_t)*cos(gamma_t)*2.2237335E-4;
	D_q_t[5][3] = sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*(-3.515445E-5)-(sin(alpha_t)*sin(gamma_t)*1.64721E-3+sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(gamma_t)*sin(alpha_t)*3.18887E-2-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.64721E-3+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*3.595415E-2+cos(alpha_t)*cos(gamma_t)*1.64721E-3+cos(beta_t)*sin(alpha_t)*6.272722E-2+cos(alpha_t)*sin(gamma_t)*3.18887E-2-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*2.621678E-2-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*3.18887E-2+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2)-(sin(gamma_t-q1_t)*cos(alpha_t)*(-4.667799999999986E-4)+cos(alpha_t)*cos(gamma_t)*3.18887E-2-cos(alpha_t)*sin(gamma_t)*1.64721E-3-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.64721E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(sin(alpha_t)*sin(gamma_t)*3.18887E-2+sin(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*3.595415E-2-cos(alpha_t)*cos(beta_t)*6.272722E-2-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(gamma_t)*sin(alpha_t)*1.64721E-3+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*3.18887E-2-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2);
	D_q_t[5][4] = cos(gamma_t)*(-1.3948862047587E-5)-sin(gamma_t)*2.7003908255589E-4-sin(gamma_t)*sin(q1_t)*5.944086627349498E-5-cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*4.3377857925795E-6-cos(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*8.3976147427365E-5+sin(3.141592653589793*(1.0/1.8E1))*cos(q1_t)*1.10853093375E-4-sin(3.141592653589793*(1.0/1.8E1))*sin(q1_t)*5.7261137625E-6-cos(gamma_t)*cos(q1_t)*5.944086627349498E-5+cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q1_t)*2.18055498525E-4-cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q2_t)*7.9952447819025E-6-cos(3.141592653589793*(1.0/1.8E1))*cos(q2_t)*sin(gamma_t)*1.54781698919175E-4-cos(3.141592653589793*(1.0/1.8E1))*cos(q1_t)*sin(q2_t)*1.1286256937211E-4-sin(3.141592653589793*(1.0/1.8E1))*cos(q1_t)*cos(q2_t)*1.1286256937211E-4+cos(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*sin(q1_t)*2.18055498525E-4+sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*sin(q2_t)*7.9952447819025E-6+cos(3.141592653589793*(1.0/1.8E1))*sin(q1_t)*sin(q2_t)*5.829913194813E-6+sin(3.141592653589793*(1.0/1.8E1))*cos(q2_t)*sin(q1_t)*5.829913194813E-6+sin(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*sin(q2_t)*1.54781698919175E-4+cos(gamma_t)*cos(q1_t)*cos(q2_t)*5.5945979761419E-5+cos(q2_t)*sin(gamma_t)*sin(q1_t)*5.5945979761419E-5-cos(3.141592653589793*(1.0/1.8E1))*cos(q2_t)*sin(gamma_t)*sin(q1_t)*2.22008273111466E-4+sin(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q1_t)*sin(q2_t)*2.22008273111466E-4+sin(3.141592653589793*(1.0/1.8E1))*sin(gamma_t)*sin(q1_t)*sin(q2_t)*2.22008273111466E-4-cos(3.141592653589793*(1.0/1.8E1))*cos(gamma_t)*cos(q1_t)*cos(q2_t)*2.22008273111466E-4;
	D_q_t[5][5] = pow(cos(beta_t)*sin(gamma_t)*9.19129199598909E18-cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*7.421932185906577E18+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*4.747759579008611E17,2.0)*1.62500339057673E-42+pow(sin(alpha_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*sin(alpha_t)*1.345401749800956E17-cos(gamma_t)*sin(alpha_t)*9.19129199598909E18-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*7.421932185906577E18+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*7.556472360886673E18+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*4.747759579008611E17+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*9.19129199598909E18,2.0)*1.62500339057673E-42+pow(sin(gamma_t-q1_t)*cos(alpha_t)*-1.345401749800956E17+cos(alpha_t)*cos(gamma_t)*9.19129199598909E18-cos(alpha_t)*sin(gamma_t)*4.747759579008611E17-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*7.421932185906577E18+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*7.556472360886673E18+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*4.747759579008611E17+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*9.19129199598909E18,2.0)*1.62500339057673E-42+3.515445E-5;
	D_q_t[5][6] = (sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(alpha_t)*sin(gamma_t)*1.64721E-3+sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(gamma_t)*sin(alpha_t)*3.18887E-2-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.64721E-3+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(-2.7E1/2.0E2)+(sin(gamma_t-q1_t)*cos(alpha_t)*4.667799999999986E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(gamma_t-q1_t)*cos(alpha_t)*(-4.667799999999986E-4)+cos(alpha_t)*cos(gamma_t)*3.18887E-2-cos(alpha_t)*sin(gamma_t)*1.64721E-3-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.64721E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(2.7E1/2.0E2)-cos(gamma_t-q1_t)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*9.445590451108341E17-cos(beta_t+3.141592653589793*(1.0/1.8E1))*9.277415232383222E17)*(cos(beta_t)*sin(gamma_t)*3.18887E-2-cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.621678E-2+cos(beta_t)*cos(gamma_t)*1.64721E-3)*3.747002708109903E-21-3.515445E-5;
	D_q_t[5][7] = (cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*(cos(beta_t)*sin(gamma_t)*3.18887E-2-cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.621678E-2+cos(beta_t)*cos(gamma_t)*1.64721E-3)*(-2.7E1/2.0E2)-cos(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(alpha_t)*sin(gamma_t)*1.64721E-3+sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(gamma_t)*sin(alpha_t)*3.18887E-2-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.64721E-3+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*3.278627369596165E-18-sin(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(gamma_t-q1_t)*cos(alpha_t)*(-4.667799999999986E-4)+cos(alpha_t)*cos(gamma_t)*3.18887E-2-cos(alpha_t)*sin(gamma_t)*1.64721E-3-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.64721E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*3.278627369596165E-18;
	D_q_t[6][0] = sin(gamma_t-q1_t)*sin(alpha_t)*6.301529999999981E-5-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3;
	D_q_t[6][1] = sin(gamma_t-q1_t)*cos(alpha_t)*(-6.301529999999981E-5)-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*3.47625E-3+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*3.5392653E-3;
	D_q_t[6][2] = cos(gamma_t-q1_t)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*9.445590451108341E17-cos(beta_t+3.141592653589793*(1.0/1.8E1))*9.277415232383222E17)*3.747002708109903E-21;
	D_q_t[6][3] = sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.515445E-5-(sin(gamma_t-q1_t)*cos(alpha_t)*4.667799999999986E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(alpha_t)*sin(gamma_t)*3.18887E-2+sin(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4-cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*3.595415E-2-cos(alpha_t)*cos(beta_t)*6.272722E-2-cos(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(gamma_t)*sin(alpha_t)*1.64721E-3+sin(gamma_t-q1_t)*cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2-sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*3.18887E-2-cos(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2)+(sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(alpha_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2+cos(alpha_t)*cos(gamma_t-q1_t)*4.667799999999986E-4+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*3.595415E-2+cos(alpha_t)*cos(gamma_t)*1.64721E-3+cos(beta_t)*sin(alpha_t)*6.272722E-2+cos(alpha_t)*sin(gamma_t)*3.18887E-2-sin(gamma_t-q1_t)*sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*2.621678E-2-cos(gamma_t)*sin(alpha_t)*sin(beta_t)*3.18887E-2+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*1.64721E-3)*(2.7E1/2.0E2);
	D_q_t[6][4] = sin(alpha_t)*(sin(gamma_t-q1_t)*cos(alpha_t)*4.667799999999986E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*(-4.683753385137379E-22)+cos(gamma_t-q1_t)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*9.445590451108341E17-cos(beta_t+3.141592653589793*(1.0/1.8E1))*9.277415232383222E17)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+cos(beta_t)*6.272722E-2+cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2-cos(gamma_t)*sin(beta_t)*3.18887E-2+sin(beta_t)*sin(gamma_t)*1.64721E-3-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*3.747002708109903E-21+cos(alpha_t)*(sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*4.683753385137379E-22;
	D_q_t[6][5] = (sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(alpha_t)*sin(gamma_t)*1.64721E-3+sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(gamma_t)*sin(alpha_t)*3.18887E-2-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.64721E-3+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(-2.7E1/2.0E2)+(sin(gamma_t-q1_t)*cos(alpha_t)*4.667799999999986E-4+sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2-sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2)*(sin(gamma_t-q1_t)*cos(alpha_t)*(-4.667799999999986E-4)+cos(alpha_t)*cos(gamma_t)*3.18887E-2-cos(alpha_t)*sin(gamma_t)*1.64721E-3-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.64721E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*(2.7E1/2.0E2)-cos(gamma_t-q1_t)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*9.445590451108341E17-cos(beta_t+3.141592653589793*(1.0/1.8E1))*9.277415232383222E17)*(cos(beta_t)*sin(gamma_t)*3.18887E-2-cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.621678E-2+cos(beta_t)*cos(gamma_t)*1.64721E-3)*3.747002708109903E-21-3.515445E-5;
	D_q_t[6][6] = pow(cos(gamma_t-q1_t),2.0)*1.8227216295E-4-cos(q2_t)*pow(cos(gamma_t-q1_t),2.0)*1.8227216295E-4+3.5183864281734E-5;
	D_q_t[6][7] = cos(gamma_t-q1_t)*(cos(q2_t)*6.867327096399238E32+sin(gamma_t-q1_t)*sin(q2_t)*5.007466555998059E32-6.991813734925732E32)*(-1.820003797445938E-37);
	D_q_t[7][0] = cos(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*3.278627369596165E-18;
	D_q_t[7][1] = sin(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*3.278627369596165E-18;
	D_q_t[7][2] = cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*4.85381025E-3+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.5392653E-3;
	D_q_t[7][3] = (sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(cos(gamma_t)*4.747759579008611E17+sin(gamma_t)*9.19129199598909E18+cos(gamma_t-q1_t)*1.345401749800956E17)*1.137502373403711E-38;
	D_q_t[7][4] = pow(cos(gamma_t-q1_t),2.0)*1.25E-4+pow(sin(gamma_t-q1_t),2.0)*1.25E-4+(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*(cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+cos(beta_t)*6.272722E-2+cos(beta_t+3.141592653589793*(1.0/1.8E1))*1.950677E-2-cos(gamma_t)*sin(beta_t)*3.18887E-2+sin(beta_t)*sin(gamma_t)*1.64721E-3-sin(gamma_t-q1_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*2.575E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*(2.7E1/2.0E2)+pow(cos(alpha_t),2.0)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*1.137502373403711E-38+pow(sin(alpha_t),2.0)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.036307817871507E19+sin(beta_t)*1.807989021555118E19+sin(beta_t+3.141592653589793*(1.0/1.8E1))*5.622443654604926E18-cos(beta_t)*sin(gamma_t)*4.747759579008611E17+sin(gamma_t-q1_t)*cos(beta_t+3.141592653589793*(1.0/1.8E1))*7.421932185906577E18-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*7.556472360886673E18+cos(beta_t)*cos(gamma_t)*9.19129199598909E18)*1.137502373403711E-38;
	D_q_t[7][5] = (cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*3.595415E-2+sin(gamma_t-q1_t)*sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*2.621678E-2)*(cos(beta_t)*sin(gamma_t)*3.18887E-2-cos(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.621678E-2+cos(beta_t)*cos(gamma_t)*1.64721E-3)*(-2.7E1/2.0E2)-cos(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(alpha_t)*sin(gamma_t)*1.64721E-3+sin(gamma_t-q1_t)*sin(alpha_t)*4.667799999999986E-4-cos(gamma_t)*sin(alpha_t)*3.18887E-2-cos(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*cos(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(alpha_t)*cos(gamma_t)*sin(beta_t)*1.64721E-3+cos(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*3.278627369596165E-18-sin(alpha_t)*(sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.480439739816438E15-sin(gamma_t-q1_t)*cos(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*1.079496051555239E15)*(sin(gamma_t-q1_t)*cos(alpha_t)*(-4.667799999999986E-4)+cos(alpha_t)*cos(gamma_t)*3.18887E-2-cos(alpha_t)*sin(gamma_t)*1.64721E-3-sin(alpha_t)*sin(beta_t+3.141592653589793*(1.0/1.8E1))*cos(gamma_t-q1_t)*2.575E-2+sin(beta_t+q2_t+3.141592653589793*(1.0/1.8E1))*sin(alpha_t)*cos(gamma_t-q1_t)*2.621678E-2+cos(gamma_t)*sin(alpha_t)*sin(beta_t)*1.64721E-3+sin(alpha_t)*sin(beta_t)*sin(gamma_t)*3.18887E-2)*3.278627369596165E-18;
	D_q_t[7][6] = cos(gamma_t-q1_t)*(cos(q2_t)*6.867327096399238E32+sin(gamma_t-q1_t)*sin(q2_t)*5.007466555998059E32-6.991813734925732E32)*(-1.820003797445938E-37);
	D_q_t[7][7] = pow(sin(gamma_t-q1_t),2.0)*9.278813973173399E-5+2.995146218000375E-4;

	// Evaluating the Coriolis matrix
	

}

void Envi_Map_Defi()
{
	// 	This function is used to define the environment map for the simulation
	// Whenever this function gets called, it will return the array with the
	// environment obstacle information
	//
	// This is the default flat ground
	// This map is defined in a polyline manner with the first value denoting
	// the line length and the second value denoting the relative angle
	Envi_Map = dlib::ones_matrix<double>(2,4);
	Envi_Map(0,0) = -5.0;					Envi_Map(0,1) = 0.0;				Envi_Map(0,2) = 4.975;				Envi_Map(0,3) = 0.0;
	Envi_Map(1,0) = 4.975;				Envi_Map(1,1) = 0.0;				Envi_Map(1,2) = 4.975;				Envi_Map(1,3) = 3.0;
}
void Envi_Map_Normal_Cal(dlib::matrix<double> &Envi_Map)
{
	// This function is used to calculate the surface normal vector and tangential vector
	int NumOfObs = Envi_Map.nr();
	int Dimension = Envi_Map.nc()/2;
	Envi_Map_Normal = dlib::zeros_matrix<double>(NumOfObs,Dimension);
	Envi_Map_Tange = dlib::zeros_matrix<double>(NumOfObs,Dimension);
	double Envi_Map_Edge_A_x, Envi_Map_Edge_A_y, Envi_Map_Edge_B_x, Envi_Map_Edge_B_y;
	double Slope_Angle;
	for (int i = 0; i < NumOfObs; i++)
	{
		Envi_Map_Edge_A_x = Envi_Map(i, 0);						Envi_Map_Edge_A_y = Envi_Map(i, 1);
		Envi_Map_Edge_B_x = Envi_Map(i, 2);						Envi_Map_Edge_B_y = Envi_Map(i, 3);
		Slope_Angle = atan2(Envi_Map_Edge_B_y - Envi_Map_Edge_A_y, Envi_Map_Edge_B_x - Envi_Map_Edge_A_x);
		Envi_Map_Tange(i,0) = cos(Slope_Angle);					Envi_Map_Tange(i,1) = sin(Slope_Angle);
		Envi_Map_Normal(i,0) = cos(Slope_Angle + pi/2.0);		Envi_Map_Normal(i,1) = sin(Slope_Angle + pi/2.0);
	}
	return;
}
void Add_Node2Tree(Tree_Node &Current_Node)
{
	// This function will add the current node to the All_Nodes vector
	All_Nodes.push_back(&Current_Node);
	Frontier_Nodes.push_back(&Current_Node);
	Frontier_Nodes_Cost.push_back(Current_Node.KE);
}
int Minimum_Index(std::vector<double> &Given_vec)
{
	int index = 0;
	for(int i = 1; i < Given_vec.size(); i++)
	{
		if(Given_vec[i] < Given_vec[index])
			index = i;
	}
	return index;
}
Tree_Node Pop_Node()
{
	// This function will pop the node outfrom the current Frontier according to the kinetic energy
	int Min_Ind = Minimum_Index(Frontier_Nodes_Cost);
	Tree_Node Current_Node = *Frontier_Nodes[Min_Ind];
	Frontier_Nodes.erase(Frontier_Nodes.begin()+Min_Ind);
	Frontier_Nodes_Cost.erase(Frontier_Nodes_Cost.begin()+Min_Ind);
	return Current_Node;
}
std::vector<double> Default_Init(const std::vector<double> &sigma_i)
{
	// This function is used to initialize the whole optimization process
	// First, is to substitute the map info into the Envi_Map matrix
	// Second, is to give the proper bounds to the variables to be optimized
	// Thrid, is to generate a kinematically feasible initial robot state

	Envi_Map_Defi();
	// Then it is to compute the normal and tangential vector of the map
	Envi_Map_Normal_Cal(Envi_Map);
	xlow_vec = dlib::zeros_matrix<double>(26,1);					xupp_vec = dlib::zeros_matrix<double>(26,1);
	ctrl_low_vec = dlib::zeros_matrix<double>(10,1);			ctrl_upp_vec = dlib::matrix<double>(10,1) ;
	xlow_vec(0) = rIxlow; 						xupp_vec(0) = rIxupp;
	xlow_vec(1) = rIylow; 						xupp_vec(1) = rIyupp;
	xlow_vec(2) = thetalow; 					xupp_vec(2) = thetaupp;
	xlow_vec(3) = q1low; 							xupp_vec(3) = q1upp;
	xlow_vec(4) = q2low; 							xupp_vec(4) = q2upp;
	xlow_vec(5) = q3low; 							xupp_vec(5) = q3upp;
	xlow_vec(6) = q4low; 							xupp_vec(6) = q4upp;
	xlow_vec(7) = q5low; 							xupp_vec(7) = q5upp;
	xlow_vec(8) = q6low; 							xupp_vec(8) = q6upp;
	xlow_vec(9) = q7low; 							xupp_vec(9) = q7upp;
	xlow_vec(10) = q8low; 						xupp_vec(10) = q8upp;
	xlow_vec(11) = q9low; 						xupp_vec(11) = q9upp;
	xlow_vec(12) = q10low; 						xupp_vec(12) = q10upp;
	xlow_vec(0+13) = rIxdotlow; 			xupp_vec(0+13) = rIxdotupp;
	xlow_vec(1+13) = rIydotlow; 			xupp_vec(1+13) = rIydotupp;
	xlow_vec(2+13) = thetadotlow; 		xupp_vec(2+13) = thetadotupp;
	xlow_vec(3+13) = q1dotlow; 				xupp_vec(3+13) = q1dotupp;
	xlow_vec(4+13) = q2dotlow; 				xupp_vec(4+13) = q2dotupp;
	xlow_vec(5+13) = q3dotlow; 				xupp_vec(5+13) = q3dotupp;
	xlow_vec(6+13) = q4dotlow; 				xupp_vec(6+13) = q4dotupp;
	xlow_vec(7+13) = q5dotlow; 				xupp_vec(7+13) = q5dotupp;
	xlow_vec(8+13) = q6dotlow; 				xupp_vec(8+13) = q6dotupp;
	xlow_vec(9+13) = q7dotlow; 				xupp_vec(9+13) = q7dotupp;
	xlow_vec(10+13) = q8dotlow; 			xupp_vec(10+13) = q8dotupp;
	xlow_vec(11+13) = q9dotlow; 			xupp_vec(11+13) = q9dotupp;
	xlow_vec(12+13) = q10dotlow; 			xupp_vec(12+13) = q10dotupp;

	ctrl_low_vec(0) = -tau1_max;			ctrl_upp_vec(0) = -ctrl_low_vec(0);
	ctrl_low_vec(1) = -tau2_max;			ctrl_upp_vec(1) = -ctrl_low_vec(1);
	ctrl_low_vec(2) = -tau3_max;			ctrl_upp_vec(2) = -ctrl_low_vec(2);
	ctrl_low_vec(3) = -tau4_max;			ctrl_upp_vec(3) = -ctrl_low_vec(3);
	ctrl_low_vec(4) = -tau5_max;			ctrl_upp_vec(4) = -ctrl_low_vec(4);
	ctrl_low_vec(5) = -tau6_max;			ctrl_upp_vec(5) = -ctrl_low_vec(5);
	ctrl_low_vec(6) = -tau7_max;			ctrl_upp_vec(6) = -ctrl_low_vec(6);
	ctrl_low_vec(7) = -tau8_max;			ctrl_upp_vec(7) = -ctrl_low_vec(7);
	ctrl_low_vec(8) = -tau9_max;			ctrl_upp_vec(8) = -ctrl_low_vec(8);
	ctrl_low_vec(9) = -tau10_max;			ctrl_upp_vec(9) = -ctrl_low_vec(9);

	vector<double> Robot_State_Init;
	ifstream Initial_Robot_State_File;              // This is to read the initial angle and angular velocities
	Initial_Robot_State_File.open("robot_angle_init.txt");
	if(Initial_Robot_State_File.is_open())
	{
		double data_each_line = 0.0;
		while(Initial_Robot_State_File>>data_each_line)
		{
			Robot_State_Init.push_back(data_each_line);
		}
		Initial_Robot_State_File.close();
	}
	else
	{
		printf("Unable to open robot_angle_init.txt file!\n");
	}
	Initial_Robot_State_File.open("robot_velocity_init.txt");
	if(Initial_Robot_State_File.is_open())
	{
		double data_each_line = 0.0;
		while(Initial_Robot_State_File>>data_each_line)
		{
			Robot_State_Init.push_back(data_each_line);
		}
		Initial_Robot_State_File.close();
	}
	else
	{
		printf("Unable to open robot_velocity_init.txt file!\n");
	}

	// Here the robot initial state has been read into the Robot_State_Init vector
	Tree_Node RootNode;
	RootNode.Node_StateNDot = StateVec2StateNDot(Robot_State_Init);
	RootNode.sigma = sigma_i;
	Structure_P.Node_i = RootNode;
	// Robot_Plot_fn(RootNode.Node_StateNDot);
	// // If the default configuration would like to be viewed
	// Robot_StateNDot Robot_StateNDot_init(Robot_State_Init);
	// std::string input_name = "init_given";
	// Robot_Plot_fn(Robot_StateNDot_init,input_name);

	Robot_State_Init = Default_Init_Opt(Robot_State_Init);
	return Robot_State_Init;
}
std::vector<double> Default_Init_Opt(std::vector<double> &Robot_State_Init)
{
		snoptProblem Default_Init_Pr;                     // This is the name of the Optimization problem
		// Allocate and initialize
		std:vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
		integer n = Robot_State_Init.size();
		Default_Init_Pr_ObjNConstraint(Robot_State_Init, ObjNConstraint_Val, ObjNConstraint_Type);
		integer neF = ObjNConstraint_Val.size();     							  // 1 objective function
		integer lenA  =  n * neF;                         // This is the number of nonzero elements in the linear part A    F(x) = f(x)+Ax

		integer *iAfun = new integer[lenA];              integer *jAvar = new integer[lenA];				doublereal *A  = new doublereal[lenA];

		integer lenG   = lenA;							 integer *iGfun = new integer[lenG];				integer *jGvar = new integer[lenG];

		doublereal *x      = new doublereal[n];			doublereal *xlow   = new doublereal[n];				doublereal *xupp   = new doublereal[n];
		doublereal *xmul   = new doublereal[n];			integer    *xstate = new    integer[n];

		doublereal *F      = new doublereal[neF];		doublereal *Flow   = new doublereal[neF];			doublereal *Fupp   = new doublereal[neF];
		doublereal *Fmul   = new doublereal[neF];		integer    *Fstate = new integer[neF];

		integer nxnames = 1;							integer nFnames = 1;						char *xnames = new char[nxnames*8];						char *Fnames = new char[nFnames*8];

		integer    ObjRow = 0;							doublereal ObjAdd = 0;

		// Set the upper and lower bounds.
		for (int i = 0; i < n; i++) {
			xlow[i] = xlow_vec(i);						xupp[i] = xupp_vec(i);						xstate[i] = 0.0;					x[i] = Robot_State_Init[i];  	// Initial guess
		}

		for(int i = 0; i<neF; i++){
			// The lower bound is the same
			Flow[i] = 0.0;
			if(ObjNConstraint_Type[i]>0)	// Inequality constraint
				{	Fupp[i] = Inf;}
			else{
				Fupp[i] = 0.0;}
			}

		// Load the data for ToyProb ...
			Default_Init_Pr.setPrintFile  ( "Default_Init_Pr.out" );
			Default_Init_Pr.setProblemSize( n, neF );
			Default_Init_Pr.setObjective  ( ObjRow, ObjAdd );
			Default_Init_Pr.setA          ( lenA, iAfun, jAvar, A );
			Default_Init_Pr.setG          ( lenG, iGfun, jGvar );
			Default_Init_Pr.setX          ( x, xlow, xupp, xmul, xstate );
			Default_Init_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
			Default_Init_Pr.setXNames     ( xnames, nxnames );
			Default_Init_Pr.setFNames     ( Fnames, nFnames );
			Default_Init_Pr.setProbName   ( "Default_Init_Pr" );
			Default_Init_Pr.setUserFun    ( Default_Init_Pr_);
		// snopta will compute the Jacobian by finite-differences.
		// The user has the option of calling  snJac  to define the
		// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
			Default_Init_Pr.computeJac    ();
			Default_Init_Pr.setIntParameter( "Derivative option", 0 );
			Default_Init_Pr.setIntParameter( "Major print level", 0 );
			Default_Init_Pr.setIntParameter( "Minor print level", 0 );
			integer Cold = 0, Basis = 1, Warm = 2;
			Default_Init_Pr.solve( Cold );

		// Take the value out from x
			for (int i = 0; i < n; i++)
			{
				Robot_State_Init[i] = x[i];
			}
		// Robot_StateNDot Init_Opt_vec(Robot_State_Init);
		// Robot_Plot_fn(Init_Opt_vec);

			delete []iAfun;  delete []jAvar;  delete []A;
			delete []iGfun;  delete []jGvar;

			delete []x;      delete []xlow;   delete []xupp;
			delete []xmul;   delete []xstate;

			delete []F;		 delete []Flow;	  delete []Fupp;
			delete []Fmul;	 delete []Fstate;

			delete []xnames; delete []Fnames;

			return Robot_State_Init;
		}

		std::vector<double> StateNDot2StateVec(const Robot_StateNDot &Robot_StateNDot_i)
		{
			std::vector<double> StateVec(26);
			StateVec[0] = Robot_StateNDot_i.rIx;
			StateVec[1] = Robot_StateNDot_i.rIy;
			StateVec[2] = Robot_StateNDot_i.theta;
			StateVec[3] = Robot_StateNDot_i.q1;
			StateVec[4] = Robot_StateNDot_i.q2;
			StateVec[5] = Robot_StateNDot_i.q3;
			StateVec[6] = Robot_StateNDot_i.q4;
			StateVec[7] = Robot_StateNDot_i.q5;
			StateVec[8] = Robot_StateNDot_i.q6;
			StateVec[9] = Robot_StateNDot_i.q7;
			StateVec[10] = Robot_StateNDot_i.q8;
			StateVec[11] = Robot_StateNDot_i.q9;
			StateVec[12] = Robot_StateNDot_i.q10;

			StateVec[0+13] = Robot_StateNDot_i.rIxdot;
			StateVec[1+13] = Robot_StateNDot_i.rIydot;
			StateVec[2+13] = Robot_StateNDot_i.thetadot;
			StateVec[3+13] = Robot_StateNDot_i.q1dot;
			StateVec[4+13] = Robot_StateNDot_i.q2dot;
			StateVec[5+13] = Robot_StateNDot_i.q3dot;
			StateVec[6+13] = Robot_StateNDot_i.q4dot;
			StateVec[7+13] = Robot_StateNDot_i.q5dot;
			StateVec[8+13] = Robot_StateNDot_i.q6dot;
			StateVec[9+13] = Robot_StateNDot_i.q7dot;
			StateVec[10+13] = Robot_StateNDot_i.q8dot;
			StateVec[11+13] = Robot_StateNDot_i.q9dot;
			StateVec[12+13] = Robot_StateNDot_i.q10dot;
			return StateVec;
		}
		Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec)
		{
			Robot_StateNDot Robot_StateNDot_i;
			Robot_StateNDot_i.rIx = StateVec[0];
			Robot_StateNDot_i.rIy = StateVec[1];
			Robot_StateNDot_i.theta = StateVec[2];
			Robot_StateNDot_i.q1 = StateVec[3];
			Robot_StateNDot_i.q2 = StateVec[4];
			Robot_StateNDot_i.q3 = StateVec[5];
			Robot_StateNDot_i.q4 = StateVec[6];
			Robot_StateNDot_i.q5 = StateVec[7];
			Robot_StateNDot_i.q6 = StateVec[8];
			Robot_StateNDot_i.q7 = StateVec[9];
			Robot_StateNDot_i.q8 = StateVec[10];
			Robot_StateNDot_i.q9 = StateVec[11];
			Robot_StateNDot_i.q10 = StateVec[12];

			Robot_StateNDot_i.rIxdot = StateVec[13];
			Robot_StateNDot_i.rIydot = StateVec[14];
			Robot_StateNDot_i.thetadot = StateVec[15];
			Robot_StateNDot_i.q1dot = StateVec[16];
			Robot_StateNDot_i.q2dot = StateVec[17];
			Robot_StateNDot_i.q3dot = StateVec[18];
			Robot_StateNDot_i.q4dot = StateVec[19];
			Robot_StateNDot_i.q5dot = StateVec[20];
			Robot_StateNDot_i.q6dot = StateVec[21];
			Robot_StateNDot_i.q7dot = StateVec[22];
			Robot_StateNDot_i.q8dot = StateVec[23];
			Robot_StateNDot_i.q9dot = StateVec[24];
			Robot_StateNDot_i.q10dot = StateVec[25];

			return Robot_StateNDot_i;

		}
		int Default_Init_Pr_(integer    *Status, integer *n,    doublereal x[],
			integer    *needF,  integer *neF,  doublereal F[],
			integer    *needG,  integer *neG,  doublereal G[],
			char       *cu,     integer *lencu,
			integer    iu[],    integer *leniu,
			doublereal ru[],    integer *lenru )
		{
			std::vector<double> ObjNConstraint_Val,ObjNConstraint_Type;

	// Initial guess of the robot configurations
			std::vector<double> Robot_State_Init = StateNDot2StateVec(Structure_P.Node_i.Node_StateNDot);
			std::vector<double> Robot_State_Opt;
			for (int i = 0; i < 26; i++)
			{
				Robot_State_Opt.push_back(x[i]);
			}

			Default_Init_Pr_ObjNConstraint(Robot_State_Opt, ObjNConstraint_Val,ObjNConstraint_Type);
			for (int i = 0; i < ObjNConstraint_Val.size(); i++){
				F[i] = ObjNConstraint_Val[i];}
				return 0;
			}
			void Default_Init_Pr_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
			{
				Robot_StateNDot StateNDot_Init_i(Opt_Seed);		dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
				End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);
				std::vector<double> Robostate_ref = StateNDot2StateVec(Structure_P.Node_i.Node_StateNDot);
				std::vector<double> Robostate_offset = Vec_Minus(Opt_Seed, Robostate_ref);
				double Robostate_offset_val = 0.0;
				for (int i = 0; i < Robostate_offset.size()/2; i++){
					Robostate_offset_val = Robostate_offset_val + Robostate_offset[i] * Robostate_offset[i];
				}
				std::vector<double> sigma = Structure_P.Node_i.sigma;
				ObjNConstraint_Val.push_back(Robostate_offset_val);
				ObjNConstraint_Type.push_back(1);

				dlib::matrix<double,6,1> End_Effector_Dist;
				std::vector<int> End_Effector_Obs(6);

				End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

				dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Matrix_result;
				std::vector<double> sigma_temp;
				sigma_temp = Sigma2Pos(sigma, 0);
				Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
				sigma_temp = Sigma2Pos(sigma, 1);
				Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
				sigma_temp = Sigma2Vel(sigma);
				Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);
	// cout<<Eqn_Pos_Matrix<<endl;				cout<<Ineqn_Pos_Matrix<<endl;			cout<<Eqn_Vel_Matrix<<endl;

	// 1. Active constraints have to be satisfied: Position and Velocity
				Matrix_result = Eqn_Pos_Matrix * End_Effector_Dist;
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	// cout<<Matrix_result<<endl;

				Matrix_result = Eqn_Vel_Matrix * End_Effector_Vel;
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	// cout<<Matrix_result<<endl;

	// // 2. Inactive constraints have to be strictly away from the obstacle
	// dlib::matrix<double> ones_vector, temp_matrix;
	// ones_vector = ONES_VECTOR_fn(6);
	// Matrix_result = Ineqn_Pos_Matrix * (End_Effector_Dist - ones_vector * mini);
	// ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

				double rAx = End_Effector_Pos(0);			double rAy = End_Effector_Pos(1);
				double rBx = End_Effector_Pos(2);			double rBy = End_Effector_Pos(3);
				double rCx = End_Effector_Pos(4);			double rCy = End_Effector_Pos(5);
				double rDx = End_Effector_Pos(6);			double rDy = End_Effector_Pos(7);
				double rEx = End_Effector_Pos(8);			double rEy = End_Effector_Pos(9);
				double rFx = End_Effector_Pos(10);			double rFy = End_Effector_Pos(11);
	//
	// ObjNConstraint_Val.push_back((-rAx + rDx - 0.2) * (-rAx + rDx - 0.2));										ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.rIx - 4.3) * (StateNDot_Init_i.rIx - 4.3));					ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q5 - 0.75) * (StateNDot_Init_i.q5 - 0.75));					ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q2 - 0.75) * (StateNDot_Init_i.q2 - 0.75));					ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.rIxdot - 0.35) * (StateNDot_Init_i.rIxdot - 0.35));			ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.theta ) * (StateNDot_Init_i.theta));				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rAx - 0.579) * (rAx - 0.579));				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rBx - 0.1109) * (rBx - 0.1109));				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rAx - 0.3359) * (rAx - 0.3359));				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rEx - 0.7897) * (rEx - 0.7897));			ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rFx - 0.9188) * (rFx - 0.9188));				ObjNConstraint_Type.push_back(0);

				double KE_init = Kinetic_Energy_fn(StateNDot_Init_i);

				ObjNConstraint_Val.push_back(10 - KE_init);			ObjNConstraint_Type.push_back(1);
				ObjNConstraint_Val.push_back(KE_init - 10);			ObjNConstraint_Type.push_back(1);
	// ObjNConstraint_Val.push_back(49.18  - KE_init);			ObjNConstraint_Type.push_back(1);
	// ObjNConstraint_Val.push_back((68.57  - KE_init) * (68.57  - KE_init));			ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(32 - KE_init);			ObjNConstraint_Type.push_back(0);

	//
				std::vector<double> vCOM_init = Ang_Vel_fn(StateNDot_Init_i, "vCOM");
	// std::vector<double> vI_init = Ang_Vel_fn(StateNDot_Init_i, "vI");
	// ObjNConstraint_Val.push_back(vCOM_init[0] - mini);	ObjNConstraint_Type.push_back(1);
	// ObjNConstraint_Val.push_back(Opt_Seed[14] - 0.5);	ObjNConstraint_Type.push_back(1);
	// ObjNConstraint_Val.push_back(Opt_Seed[14] - Opt_Seed[15]);	ObjNConstraint_Type.push_back(1);

	// std::vector<double> vCOM_ref = Ang_Vel_fn(Structure_P.Node_i.Node_StateNDot, "vCOM");
	// ObjNConstraint_Val.push_back((rBx - 3.988) * (rBx - 3.988));	ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rDx - 4.418) * (rDx - 4.418));	ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rEx - 4.966) * (rEx - 4.966));	ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((rEy - 0.6164) * (rEy - 0.6164));	ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back(vCOM_init[0] - 1.5);
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.rIx - 0.25)*(StateNDot_Init_i.q10dot - 0.25));
	// ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q10dot - 3.0)*(StateNDot_Init_i.q10dot - 3.0));
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q1dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q2dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q3dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q4dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q5dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q6dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q7dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q8dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q9dot);				ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q10dot);				ObjNConstraint_Type.push_back(0);

	// std::vector<double> sigma_i = sigma;
	// std::vector<double> sigma_i_child = sigma;
	// std::vector<double> sigma_maint(12);
	//
	// sigma_maint[0] = sigma_i[0] * sigma_i_child[0];		sigma_maint[1] = sigma_i[0] * sigma_i_child[0];
	// sigma_maint[2] = sigma_i[0] * sigma_i_child[0];		sigma_maint[3] = sigma_i[0] * sigma_i_child[0];
	// sigma_maint[4] = sigma_i[1] * sigma_i_child[1];		sigma_maint[5] = sigma_i[1] * sigma_i_child[1];
	// sigma_maint[6] = sigma_i[1] * sigma_i_child[1];		sigma_maint[7] = sigma_i[1] * sigma_i_child[1];
	// sigma_maint[8] = sigma_i[2] * sigma_i_child[2];		sigma_maint[9] = sigma_i[2] * sigma_i_child[2];
	// sigma_maint[10] = sigma_i[3] * sigma_i_child[3];	sigma_maint[11] = sigma_i[3] * sigma_i_child[3];
	//
	// dlib::matrix<double> Maint_Matrix = Diag_Matrix_fn(sigma_maint);
	// dlib::matrix<double> End_Effector_PosDlib = End_Effector_Pos;
	// dlib::matrix<double> End_Effector_Pos_ref = Structure_P.Node_i.End_Effector_Pos;
	//
	// dlib::matrix<double> Matrix_Minus_result = End_Effector_PosDlib - End_Effector_Pos_ref;
	// Matrix_result = Maint_Matrix * Matrix_Minus_result;
	// ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	// // 3. Middle joints have to be strictly away from the obs
	// temp_matrix = Middle_Joint_Obs_Dist_Fn(StateNDot_Init_i);
	// Matrix_result = temp_matrix - ones_vector * mini;
	// ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

	// double KE_Init = Kinetic_Energy_fn(StateNDot_Init_i);
	// ObjNConstraint_Val.push_back(KE_Init - 20);
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back(Opt_Seed[13] - 0.35);
	// ObjNConstraint_Type.push_back(1);

				return;
			}
			dlib::matrix<double> Middle_Joint_Obs_Dist_Fn(Robot_StateNDot &StateNDot_Init_i){
				dlib::matrix<double> Middle_Joint_Obs_Dist_Matrix;
				Middle_Joint_Obs_Dist_Matrix = dlib::ones_matrix<double>(6,1);
				std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
				std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
				std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
				std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
				std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
				std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");
				int Obs_Dist_Index;		double rH_Obs, rK_Obs, rM_Obs, rN_Obs, rI_Obs, rT_Obs;

				Obs_Dist_Fn(rH, rH_Obs, Obs_Dist_Index, "None");
				Obs_Dist_Fn(rK, rK_Obs, Obs_Dist_Index, "None");
				Obs_Dist_Fn(rM, rM_Obs, Obs_Dist_Index, "None");
				Obs_Dist_Fn(rN, rN_Obs, Obs_Dist_Index, "None");
				Obs_Dist_Fn(rI, rI_Obs, Obs_Dist_Index, "None");
				Obs_Dist_Fn(rT, rT_Obs, Obs_Dist_Index, "None");

				Middle_Joint_Obs_Dist_Matrix(0) = rH_Obs;
				Middle_Joint_Obs_Dist_Matrix(1) = rK_Obs;
				Middle_Joint_Obs_Dist_Matrix(2) = rM_Obs;
				Middle_Joint_Obs_Dist_Matrix(3) = rN_Obs;
				Middle_Joint_Obs_Dist_Matrix(4) = rI_Obs;
				Middle_Joint_Obs_Dist_Matrix(5) = rT_Obs;
				return Middle_Joint_Obs_Dist_Matrix;}
				dlib::matrix<double> ONES_VECTOR_fn(int Dim){
					dlib::matrix<double> Ones_vector;
					Ones_vector = dlib::ones_matrix<double>(Dim,1);
					return Ones_vector;}
					void ObjNConstraint_ValNType_Update(dlib::matrix<double> &Matrix_result, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type, int Constraint_Type)
					{
	// This function is used to update the value of ObjNConstraint_Val and Type
						for (int i = 0; i < Matrix_result.nr(); i++)
						{
							ObjNConstraint_Val.push_back(Matrix_result(i));
							ObjNConstraint_Type.push_back(Constraint_Type);
						}
					}
					void Obs_Dist_Fn(std::vector<double> &r_Pos, double &Obs_Dist, int &Obs_Dist_Index, const char* name)
					{
	// 	This function is used to calculate the relative distance between the robot end effector and the nearby environment
						double Temp_Edge_x, Temp_Edge_y, Temp_offset_x, Temp_offset_y, Normal_vector_i_x, Normal_vector_i_y;
						std::vector<double> Obs_Dist_vec;
						for (int i = 0; i < Envi_Map.nr(); i++)
						{
							Temp_Edge_x = Envi_Map(i,0);
							Temp_Edge_y = Envi_Map(i,1);
							Temp_offset_x = r_Pos[0] - Temp_Edge_x;
							Temp_offset_y = r_Pos[1] - Temp_Edge_y;
							Normal_vector_i_x = Envi_Map_Normal(i,0);
							Normal_vector_i_y = Envi_Map_Normal(i,1);
							Obs_Dist_vec.push_back(Temp_offset_x * Normal_vector_i_x + Temp_offset_y * Normal_vector_i_y);
							if(strcmp(name,"floor")==0)
							{
			// In this case, we are talking about the robot foot, it can only make contact with the floor
			// So the iteration will terminate after the first run
								break;
							}
							if(strcmp(name,"wall")==0)
							{
								Obs_Dist_vec[0] = 10;
							}
						}
						Obs_Dist_Index = Minimum_Index(Obs_Dist_vec);
						Obs_Dist = Obs_Dist_vec[Obs_Dist_Index];
						return;
					}
Robot_StateNDot::Robot_StateNDot(){// A default constructor
	rIx = 0;			rIy = 0.7230;			theta = -0.0900;
	q1 = 0.3768;		q2 = 0.0045;			q3 = -0.2913;			q4 = -1.0015;			q5 = 0.1500;
	q6 = 0.2698;		q7 = -0.6600;			q8 = -0.6251;			q9 = 0.6900;			q10 = -0.2951;
	rIxdot = 0.2000;	rIydot = -0.0605;		thetadot = -0.2100;
	q1dot = -0.1239;	q2dot = 1.3108;			q3dot = -0.9768;		q4dot = -1.4999;		q5dot = 2.0000;
	q6dot = -1.2999;	q7dot = 1.0000;			q8dot = -2.0000;		q9dot = -1.5708;		q10dot = -1.5000;
}
Robot_StateNDot::Robot_StateNDot(std::vector<double> &Robot_AngleNRate){	// An evaluated constructor
	rIx = Robot_AngleNRate[0];		rIy = Robot_AngleNRate[1];		theta = Robot_AngleNRate[2];
	q1 = Robot_AngleNRate[3];		q2 = Robot_AngleNRate[4];		q3 = Robot_AngleNRate[5];		q4 = Robot_AngleNRate[6];		q5 = Robot_AngleNRate[7];
	q6 = Robot_AngleNRate[8];		q7 = Robot_AngleNRate[9];		q8 = Robot_AngleNRate[10];		q9 = Robot_AngleNRate[11];		q10 = Robot_AngleNRate[12];

	rIxdot = Robot_AngleNRate[13];	rIydot = Robot_AngleNRate[14];	thetadot = Robot_AngleNRate[15];
	q1dot = Robot_AngleNRate[16];	q2dot = Robot_AngleNRate[17];	q3dot = Robot_AngleNRate[18];	q4dot = Robot_AngleNRate[19];	q5dot = Robot_AngleNRate[20];
	q6dot = Robot_AngleNRate[21];	q7dot = Robot_AngleNRate[22];	q8dot = Robot_AngleNRate[23];	q9dot = Robot_AngleNRate[24];	q10dot = Robot_AngleNRate[25];
}
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i)
{
	// This function is used to plot the robot configuration
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
	std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
	std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");
	std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
	std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

	std::vector<double> x(2), y(2);
	// AB
	x.at(0) = rA[0];	x.at(1) = rB[0];	y.at(0) = rA[1];	y.at(1) = rB[1];	plt::plot(x,y);
	// CD
	x.at(0) = rC[0];	x.at(1) = rD[0];	y.at(0) = rC[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// CJ
	x.at(0) = rC[0];	x.at(1) = rJ[0];	y.at(0) = rC[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// JD
	x.at(0) = rJ[0];	x.at(1) = rD[0];	y.at(0) = rJ[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// KJ
	x.at(0) = rK[0];	x.at(1) = rJ[0];	y.at(0) = rK[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// IK
	x.at(0) = rI[0];	x.at(1) = rK[0];	y.at(0) = rI[1];	y.at(1) = rK[1];	plt::plot(x,y);
	// BG
	x.at(0) = rB[0];	x.at(1) = rG[0];	y.at(0) = rB[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// AG
	x.at(0) = rA[0];	x.at(1) = rG[0];	y.at(0) = rA[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// GH
	x.at(0) = rG[0];	x.at(1) = rH[0];	y.at(0) = rG[1];	y.at(1) = rH[1];	plt::plot(x,y);
	// HI
	x.at(0) = rH[0];	x.at(1) = rI[0];	y.at(0) = rH[1];	y.at(1) = rI[1];	plt::plot(x,y);
	// IT
	x.at(0) = rI[0];	x.at(1) = rT[0];	y.at(0) = rI[1];	y.at(1) = rT[1];	plt::plot(x,y);
	// LM
	x.at(0) = rL[0];	x.at(1) = rM[0];	y.at(0) = rL[1];	y.at(1) = rM[1];	plt::plot(x,y);
	// ME
	x.at(0) = rM[0];	x.at(1) = rE[0];	y.at(0) = rM[1];	y.at(1) = rE[1];	plt::plot(x,y);
	// LN
	x.at(0) = rL[0];	x.at(1) = rN[0];	y.at(0) = rL[1];	y.at(1) = rN[1];	plt::plot(x,y);
	// NF
	x.at(0) = rN[0];	x.at(1) = rF[0];	y.at(0) = rN[1];	y.at(1) = rF[1];	plt::plot(x,y);
	const char* filename = "./init.png";
	std::cout<<"Saving result to "<<filename<<std::endl;
	plt::save(filename);
}
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name)
{
	// This function is used to plot the robot configuration
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
	std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
	std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");
	std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
	std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

	std::vector<double> x(2), y(2);
	// AB
	x.at(0) = rA[0];	x.at(1) = rB[0];	y.at(0) = rA[1];	y.at(1) = rB[1];	plt::plot(x,y);
	// CD
	x.at(0) = rC[0];	x.at(1) = rD[0];	y.at(0) = rC[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// CJ
	x.at(0) = rC[0];	x.at(1) = rJ[0];	y.at(0) = rC[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// JD
	x.at(0) = rJ[0];	x.at(1) = rD[0];	y.at(0) = rJ[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// KJ
	x.at(0) = rK[0];	x.at(1) = rJ[0];	y.at(0) = rK[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// IK
	x.at(0) = rI[0];	x.at(1) = rK[0];	y.at(0) = rI[1];	y.at(1) = rK[1];	plt::plot(x,y);
	// BG
	x.at(0) = rB[0];	x.at(1) = rG[0];	y.at(0) = rB[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// AG
	x.at(0) = rA[0];	x.at(1) = rG[0];	y.at(0) = rA[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// GH
	x.at(0) = rG[0];	x.at(1) = rH[0];	y.at(0) = rG[1];	y.at(1) = rH[1];	plt::plot(x,y);
	// HI
	x.at(0) = rH[0];	x.at(1) = rI[0];	y.at(0) = rH[1];	y.at(1) = rI[1];	plt::plot(x,y);
	// IT
	x.at(0) = rI[0];	x.at(1) = rT[0];	y.at(0) = rI[1];	y.at(1) = rT[1];	plt::plot(x,y);
	// LM
	x.at(0) = rL[0];	x.at(1) = rM[0];	y.at(0) = rL[1];	y.at(1) = rM[1];	plt::plot(x,y);
	// ME
	x.at(0) = rM[0];	x.at(1) = rE[0];	y.at(0) = rM[1];	y.at(1) = rE[1];	plt::plot(x,y);
	// LN
	x.at(0) = rL[0];	x.at(1) = rN[0];	y.at(0) = rL[1];	y.at(1) = rN[1];	plt::plot(x,y);
	// NF
	x.at(0) = rN[0];	x.at(1) = rF[0];	y.at(0) = rN[1];	y.at(1) = rF[1];	plt::plot(x,y);

	std::string pre_filename = "./";
	std::string post_filename = ".png";
	std::string filename = pre_filename + name + post_filename;

	std::cout<<"Saving result to "<<filename<<std::endl;
	plt::save(filename);
}
dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;							double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;						double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;							double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;							double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;							double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;							double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	dlib::matrix<double>  T;									T = dlib::zeros_matrix<double>(13,13);
	T(0,0) = 1.13E2/2.0;
	T(0,2) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q4+q5+theta)*(9.1E1/1.0E2)-cos(q7+q8+theta)*(9.9E1/1.6E2)-cos(q9+q10+theta)*(9.9E1/1.6E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-cos(q4+theta)*(3.9E1/2.0E1)-cos(q7+theta)*(1.9E1/1.6E1)-cos(q9+theta)*(1.9E1/1.6E1)+cos(theta)*(2.97E2/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,3) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,4) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,5) = cos(q1+q2+q3+theta)*(-2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,6) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q4+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,7) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,8) = cos(q4+q5+q6+theta)*(-2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(0,9) = cos(q7+q8+theta)*(-9.9E1/1.6E2)-cos(q7+theta)*(1.9E1/1.6E1);
	T(0,10) = cos(q7+q8+theta)*(-9.9E1/1.6E2);
	T(0,11) = cos(q9+q10+theta)*(-9.9E1/1.6E2)-cos(q9+theta)*(1.9E1/1.6E1);
	T(0,12) = cos(q9+q10+theta)*(-9.9E1/1.6E2);
	T(1,1) = 1.13E2/2.0;
	T(1,2) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q9+q10+theta)*(9.9E1/1.6E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sin(q4+theta)*(3.9E1/2.0E1)+sin(q7+theta)*(1.9E1/1.6E1)+sin(q9+theta)*(1.9E1/1.6E1)-sin(theta)*(2.97E2/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,3) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,4) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,5) = cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,6) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,7) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,8) = cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,9) = sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q7+theta)*(1.9E1/1.6E1);
	T(1,10) = sin(q7+q8+theta)*(9.9E1/1.6E2);
	T(1,11) = sin(q9+q10+theta)*(9.9E1/1.6E2)+sin(q9+theta)*(1.9E1/1.6E1);
	T(1,12) = sin(q9+q10+theta)*(9.9E1/1.6E2);
	T(2,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q4+q5+theta)*(9.1E1/1.0E2)-cos(q7+q8+theta)*(9.9E1/1.6E2)-cos(q9+q10+theta)*(9.9E1/1.6E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-cos(q4+theta)*(3.9E1/2.0E1)-cos(q7+theta)*(1.9E1/1.6E1)-cos(q9+theta)*(1.9E1/1.6E1)+cos(theta)*(2.97E2/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(2,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q9+q10+theta)*(9.9E1/1.6E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sin(q4+theta)*(3.9E1/2.0E1)+sin(q7+theta)*(1.9E1/1.6E1)+sin(q9+theta)*(1.9E1/1.6E1)-sin(theta)*(2.97E2/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(2,2) = cos(q2+q3)*(1.3E1/2.5E2)+cos(q5+q6)*(1.3E1/2.5E2)-cos(q7+q8)*6.80625E-1-cos(q9+q10)*6.80625E-1-sin(q2+q3)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-cos(q7)*(2.09E2/1.6E2)+cos(q8)*(9.9E1/3.2E2)-cos(q9)*(2.09E2/1.6E2)+cos(q10)*(9.9E1/3.2E2)-sin(q3)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/2.5E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/2.5E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.171739583333333E1;
	T(2,3) = cos(q2+q3)*(1.3E1/2.5E2)-sin(q2+q3)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
	T(2,4) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(2,5) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(2,6) = cos(q5+q6)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
	T(2,7) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(2,8) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(2,9) = cos(q7+q8)*(-3.403125E-1)-cos(q7)*(2.09E2/3.2E2)+cos(q8)*(9.9E1/3.2E2)+1.045989583333333;
	T(2,10) = cos(q7+q8)*(-3.403125E-1)+cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
	T(2,11) = cos(q9+q10)*(-3.403125E-1)-cos(q9)*(2.09E2/3.2E2)+cos(q10)*(9.9E1/3.2E2)+1.045989583333333;
	T(2,12) = cos(q9+q10)*(-3.403125E-1)+cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
	T(3,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-cos(q1+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(3,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(3,2) = cos(q2+q3)*(1.3E1/2.5E2)-sin(q2+q3)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
	T(3,3) = cos(q2+q3)*(1.3E1/2.5E2)-sin(q2+q3)*(1.3E1/2.5E2)+cos(q2)*5.915E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
	T(3,4) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(3,5) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(4,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(4,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(4,2) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(4,3) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q2)*2.9575E-1+cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(4,4) = cos(q3)*(1.3E1/2.5E2)-sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+5.29375E-1;
	T(4,5) = cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
	T(5,0) = cos(q1+q2+q3+theta)*(-2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(5,1) = cos(q1+q2+q3+theta)*(2.0/2.5E1)+sin(q1+q2+q3+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(5,2) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(5,3) = cos(q2+q3)*(1.3E1/5.0E2)-sin(q2+q3)*(1.3E1/5.0E2)+cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(5,4) = cos(q3)*(1.3E1/5.0E2)-sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
	T(5,5) = sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
	T(6,0) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-cos(q4+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(6,1) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+theta)*(3.9E1/2.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(6,2) = cos(q5+q6)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
	T(6,6) = cos(q5+q6)*(1.3E1/2.5E2)-sin(q5+q6)*(1.3E1/2.5E2)+cos(q5)*5.915E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+1.409583333333333;
	T(6,7) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(6,8) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(7,0) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(7,1) = sin(q4+q5+theta)*(9.1E1/1.0E2)+cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(7,2) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(7,6) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q5)*2.9575E-1+cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+5.29375E-1;
	T(7,7) = cos(q6)*(1.3E1/2.5E2)-sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+5.29375E-1;
	T(7,8) = cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
	T(8,0) = cos(q4+q5+q6+theta)*(-2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(8,1) = cos(q4+q5+q6+theta)*(2.0/2.5E1)+sin(q4+q5+q6+theta)*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(8,2) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(8,6) = cos(q5+q6)*(1.3E1/5.0E2)-sin(q5+q6)*(1.3E1/5.0E2)+cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+8.575E-2;
	T(8,7) = cos(q6)*(1.3E1/5.0E2)-sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
	T(8,8) = sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+8.575E-2;
	T(9,0) = cos(q7+q8+theta)*(-9.9E1/1.6E2)-cos(q7+theta)*(1.9E1/1.6E1);
	T(9,1) = sin(q7+q8+theta)*(9.9E1/1.6E2)+sin(q7+theta)*(1.9E1/1.6E1);
	T(9,2) = cos(q7+q8)*(-3.403125E-1)-cos(q7)*(2.09E2/3.2E2)+cos(q8)*(9.9E1/3.2E2)+1.045989583333333;
	T(9,9) = cos(q8)*(9.9E1/3.2E2)+1.045989583333333;
	T(9,10) = cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
	T(10,0) = cos(q7+q8+theta)*(-9.9E1/1.6E2);
	T(10,1) = sin(q7+q8+theta)*(9.9E1/1.6E2);
	T(10,2) = cos(q7+q8)*(-3.403125E-1)+cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
	T(10,9) = cos(q8)*(9.9E1/6.4E2)+6.0328125E-1;
	T(10,10) = 6.0328125E-1;
	T(11,0) = cos(q9+q10+theta)*(-9.9E1/1.6E2)-cos(q9+theta)*(1.9E1/1.6E1);
	T(11,1) = sin(q9+q10+theta)*(9.9E1/1.6E2)+sin(q9+theta)*(1.9E1/1.6E1);
	T(11,2) = cos(q9+q10)*(-3.403125E-1)-cos(q9)*(2.09E2/3.2E2)+cos(q10)*(9.9E1/3.2E2)+1.045989583333333;
	T(11,11) = cos(q10)*(9.9E1/3.2E2)+1.045989583333333;
	T(11,12) = cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
	T(12,0) = cos(q9+q10+theta)*(-9.9E1/1.6E2);
	T(12,1) = sin(q9+q10+theta)*(9.9E1/1.6E2);
	T(12,2) = cos(q9+q10)*(-3.403125E-1)+cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
	T(12,11) = cos(q10)*(9.9E1/6.4E2)+6.0328125E-1;
	T(12,12) = 6.0328125E-1;

	return T;
}
dlib::matrix<double> B_q_fn()
{
	dlib::matrix<double> B_q;
	B_q = dlib::zeros_matrix<double>(13,10);
	B_q(3,0) = 1.0;
	B_q(4,1) = 1.0;
	B_q(5,2) = 1.0;
	B_q(6,3) = 1.0;
	B_q(7,4) = 1.0;
	B_q(8,5) = 1.0;
	B_q(9,6) = 1.0;
	B_q(10,7) = 1.0;
	B_q(11,8) = 1.0;
	B_q(12,9) = 1.0;
	return B_q;
}
dlib::matrix<double> C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;					double rIy = Robot_StateNDot_i.rIy;						double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;					double q2 = Robot_StateNDot_i.q2;						double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;					double q5 = Robot_StateNDot_i.q5;						double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;					double q8 = Robot_StateNDot_i.q8;						double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	double rIxdot = Robot_StateNDot_i.rIxdot;			double rIydot = Robot_StateNDot_i.rIydot;				double thetadot = Robot_StateNDot_i.thetadot;
	double q1dot = Robot_StateNDot_i.q1dot;				double q2dot = Robot_StateNDot_i.q2dot;					double q3dot = Robot_StateNDot_i.q3dot;
	double q4dot = Robot_StateNDot_i.q4dot;				double q5dot = Robot_StateNDot_i.q5dot;					double q6dot = Robot_StateNDot_i.q6dot;
	double q7dot = Robot_StateNDot_i.q7dot;				double q8dot = Robot_StateNDot_i.q8dot;					double q9dot = Robot_StateNDot_i.q9dot;
	double q10dot = Robot_StateNDot_i.q10dot;

	dlib::matrix<double> T;
	T = dlib::zeros_matrix<double>(13,1);
	T(0,0) = (thetadot*thetadot)*sin(theta)*(-2.97E2/2.0E1)+(q10dot*q10dot)*sin(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(q2dot*q2dot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(q4dot*q4dot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(q5dot*q5dot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(q7dot*q7dot)*sin(q7+q8+theta)*(9.9E1/1.6E2)+(q8dot*q8dot)*sin(q7+q8+theta)*(9.9E1/1.6E2)+(q9dot*q9dot)*sin(q9+q10+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*sin(q7+q8+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*sin(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*sin(q1+theta)*(3.9E1/2.0E1)+(q4dot*q4dot)*sin(q4+theta)*(3.9E1/2.0E1)+(q7dot*q7dot)*sin(q7+theta)*(1.9E1/1.6E1)+(q9dot*q9dot)*sin(q9+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*sin(q1+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*sin(q4+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*sin(q7+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*sin(q9+theta)*(1.9E1/1.6E1)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+q10dot*q9dot*sin(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*q8dot*sin(q7+q8+theta)*(9.9E1/8.0E1)+q10dot*thetadot*sin(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*thetadot*sin(q7+q8+theta)*(9.9E1/8.0E1)+q8dot*thetadot*sin(q7+q8+theta)*(9.9E1/8.0E1)+q9dot*thetadot*sin(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*q2dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*sin(q1+theta)*(3.9E1/1.0E1)+q4dot*thetadot*sin(q4+theta)*(3.9E1/1.0E1)+q7dot*thetadot*sin(q7+theta)*(1.9E1/8.0)+q9dot*thetadot*sin(q9+theta)*(1.9E1/8.0)+sqrt(4.1E1)*q1dot*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q3dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q6dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1);
	T(1,0) = (thetadot*thetadot)*cos(theta)*(-2.97E2/2.0E1)+(q10dot*q10dot)*cos(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(q2dot*q2dot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(q4dot*q4dot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(q5dot*q5dot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(q7dot*q7dot)*cos(q7+q8+theta)*(9.9E1/1.6E2)+(q8dot*q8dot)*cos(q7+q8+theta)*(9.9E1/1.6E2)+(q9dot*q9dot)*cos(q9+q10+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*cos(q7+q8+theta)*(9.9E1/1.6E2)+(thetadot*thetadot)*cos(q9+q10+theta)*(9.9E1/1.6E2)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)-(q1dot*q1dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q2dot*q2dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q3dot*q3dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q4dot*q4dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(q5dot*q5dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(q6dot*q6dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(thetadot*thetadot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(thetadot*thetadot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*cos(q1+theta)*(3.9E1/2.0E1)+(q4dot*q4dot)*cos(q4+theta)*(3.9E1/2.0E1)+(q7dot*q7dot)*cos(q7+theta)*(1.9E1/1.6E1)+(q9dot*q9dot)*cos(q9+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*cos(q1+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*cos(q4+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*cos(q7+theta)*(1.9E1/1.6E1)+(thetadot*thetadot)*cos(q9+theta)*(1.9E1/1.6E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+q10dot*q9dot*cos(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*q8dot*cos(q7+q8+theta)*(9.9E1/8.0E1)+q10dot*thetadot*cos(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*thetadot*cos(q7+q8+theta)*(9.9E1/8.0E1)+q8dot*thetadot*cos(q7+q8+theta)*(9.9E1/8.0E1)+q9dot*thetadot*cos(q9+q10+theta)*(9.9E1/8.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)-q1dot*q2dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q1dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q2dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q4dot*q5dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q4dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q5dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q1dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q2dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q3dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q4dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q5dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q6dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+theta)*(3.9E1/1.0E1)+q4dot*thetadot*cos(q4+theta)*(3.9E1/1.0E1)+q7dot*thetadot*cos(q7+theta)*(1.9E1/8.0)+q9dot*thetadot*cos(q9+theta)*(1.9E1/8.0)+sqrt(4.1E1)*q1dot*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q3dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q6dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+5.54265E2;
	T(2,0) = sin(q1+q2+theta)*8.9271+sin(q4+q5+theta)*8.9271+sin(q7+q8+theta)*6.0699375+sin(q9+q10+theta)*6.0699375+sin(q1+theta)*1.91295E1+sin(q4+theta)*1.91295E1+sin(q7+theta)*1.1649375E1+sin(q9+theta)*1.1649375E1-sin(theta)*1.456785E2-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q10dot*q10dot)*sin(q10)*(9.9E1/6.4E2)-(q2dot*q2dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*sin(q7)*(2.09E2/3.2E2)-(q8dot*q8dot)*sin(q8)*(9.9E1/6.4E2)+(q9dot*q9dot)*sin(q9)*(2.09E2/3.2E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1-(q2dot*q2dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q5+q6)*(1.3E1/5.0E2)+(q10dot*q10dot)*sin(q9+q10)*3.403125E-1-(q2dot*q2dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*sin(q5+q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*sin(q7+q8)*3.403125E-1+(q8dot*q8dot)*sin(q7+q8)*3.403125E-1+(q9dot*q9dot)*sin(q9+q10)*3.403125E-1-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q10dot*q9dot*sin(q10)*(9.9E1/3.2E2)-q1dot*q2dot*sin(q2)*5.915E-1-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5)*5.915E-1-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q7dot*q8dot*sin(q8)*(9.9E1/3.2E2)-q10dot*thetadot*sin(q10)*(9.9E1/3.2E2)-q2dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)+q7dot*thetadot*sin(q7)*(2.09E2/1.6E2)-q8dot*thetadot*sin(q8)*(9.9E1/3.2E2)+q9dot*thetadot*sin(q9)*(2.09E2/1.6E2)-sqrt(4.1E1)*(q2dot*q2dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q3dot*q3dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q5dot*q5dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q4dot*q5dot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q2dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q5dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q10dot*q9dot*sin(q9+q10)*6.80625E-1-q1dot*q2dot*sin(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)+q7dot*q8dot*sin(q7+q8)*6.80625E-1+q10dot*thetadot*sin(q9+q10)*6.80625E-1-q2dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)+q7dot*thetadot*sin(q7+q8)*6.80625E-1+q8dot*thetadot*sin(q7+q8)*6.80625E-1+q9dot*thetadot*sin(q9+q10)*6.80625E-1-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q2dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q5dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(3,0) = sin(q1+q2+theta)*8.9271+sin(q1+theta)*1.91295E1-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q2dot*q2dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1-(q2dot*q2dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q2dot*q2dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q1dot*q2dot*sin(q2)*5.915E-1-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q2dot*q2dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q3dot*q3dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q2dot*sin(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q2dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(4,0) = sin(q1+q2+theta)*8.9271-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q1dot*q1dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2)*2.9575E-1+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(5,0) = (q1dot*q1dot)*cos(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*cos(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q1dot*q1dot)*sin(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*sin(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q3-8.960553845713439E-1)*6.5E-3+q1dot*q2dot*cos(q3)*(1.3E1/2.5E2)+q1dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q2dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q1dot*q2dot*sin(q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q3)*(1.3E1/2.5E2)+q2dot*thetadot*sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q2dot*q2dot)*sin(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*q2dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q2dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(6,0) = sin(q4+q5+theta)*8.9271+sin(q4+theta)*1.91295E1-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q5dot*q5dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1-(q5dot*q5dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5)*5.915E-1-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q5dot*q5dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-q4dot*q5dot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q5dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(7,0) = sin(q4+q5+theta)*8.9271-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q4dot*q4dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5)*2.9575E-1+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)+(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(8,0) = (q4dot*q4dot)*cos(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*cos(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q4dot*q4dot)*sin(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*sin(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)+(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q6-8.960553845713439E-1)*6.5E-3+q4dot*q5dot*cos(q6)*(1.3E1/2.5E2)+q4dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q5dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q4dot*q5dot*sin(q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q6)*(1.3E1/2.5E2)+q5dot*thetadot*sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q5dot*q5dot)*sin(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*q5dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q5dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(9,0) = sin(q7+q8+theta)*6.0699375+sin(q7+theta)*1.1649375E1-(q8dot*q8dot)*sin(q8)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q7)*(2.09E2/3.2E2)-(thetadot*thetadot)*sin(q7+q8)*3.403125E-1-q7dot*q8dot*sin(q8)*(9.9E1/3.2E2)-q8dot*thetadot*sin(q8)*(9.9E1/3.2E2);
	T(10,0) = sin(q7+q8+theta)*6.0699375+(q7dot*q7dot)*sin(q8)*(9.9E1/6.4E2)+(thetadot*thetadot)*sin(q8)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q7+q8)*3.403125E-1+q7dot*thetadot*sin(q8)*(9.9E1/3.2E2);
	T(11,0) = sin(q9+q10+theta)*6.0699375+sin(q9+theta)*1.1649375E1-(q10dot*q10dot)*sin(q10)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q9)*(2.09E2/3.2E2)-(thetadot*thetadot)*sin(q9+q10)*3.403125E-1-q10dot*q9dot*sin(q10)*(9.9E1/3.2E2)-q10dot*thetadot*sin(q10)*(9.9E1/3.2E2);
	T(12,0) = sin(q9+q10+theta)*6.0699375+(q9dot*q9dot)*sin(q10)*(9.9E1/6.4E2)+(thetadot*thetadot)*sin(q10)*(9.9E1/6.4E2)-(thetadot*thetadot)*sin(q9+q10)*3.403125E-1+q9dot*thetadot*sin(q10)*(9.9E1/3.2E2);

	return T;
}
dlib::matrix<double> Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;
	double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	dlib::matrix<double> T;
	T = dlib::zeros_matrix<double>(12,13);
	T(0,0) = 1.0;
	T(0,2) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(0,3) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(0,4) = sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(0,5) = sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
	T(1,1) = 1.0;
	T(1,2) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(1,3) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(1,4) = cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(1,5) = sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
	T(2,0) = 1.0;
	T(2,2) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
	T(2,3) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
	T(2,4) = sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
	T(2,5) = sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(-1.0/1.0E1);
	T(3,1) = 1.0;
	T(3,2) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
	T(3,3) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
	T(3,4) = cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
	T(3,5) = sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(-1.0/1.0E1);
	T(4,0) = 1.0;
	T(4,2) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(4,6) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(4,7) = sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(4,8) = sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
	T(5,1) = 1.0;
	T(5,2) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(5,6) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(5,7) = cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	T(5,8) = sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
	T(6,0) = 1.0;
	T(6,2) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
	T(6,6) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
	T(6,7) = sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
	T(6,8) = sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(-1.0/1.0E1);
	T(7,1) = 1.0;
	T(7,2) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
	T(7,6) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
	T(7,7) = cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
	T(7,8) = sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(-1.0/1.0E1);
	T(8,0) = 1.0;
	T(8,2) = sin(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
	T(8,9) = sin(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
	T(8,10) = sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
	T(9,1) = 1.0;
	T(9,2) = cos(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
	T(9,9) = cos(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
	T(9,10) = cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
	T(10,0) = 1.0;
	T(10,2) = sin(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
	T(10,11) = sin(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
	T(10,12) = sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
	T(11,1) = 1.0;
	T(11,2) = cos(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
	T(11,11) = cos(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
	T(11,12) = cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
	return T;
}
std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s)
{

	std::vector<double> r_pos(2);
	r_pos[0] = 0.0;
	r_pos[1] = 0.0;

	double rIx = Robot_StateNDot_i.rIx;							double rIy = Robot_StateNDot_i.rIy;							double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;							double q2 = Robot_StateNDot_i.q2;							double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;							double q5 = Robot_StateNDot_i.q5;							double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;							double q8 = Robot_StateNDot_i.q8;							double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	if (strcmp(s,"rA")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"rB")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
		r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"rC")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"rD")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
		r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"rE")==0)
	{
		r_pos[0] = rIx+cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
		r_pos[1] = rIy-sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rF")==0)
	{
		r_pos[0] = rIx+cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
		r_pos[1] = rIy-sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rG")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
		r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rH")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
		r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rI")==0)
	{
		r_pos[0] = rIx;
		r_pos[1] = rIy;
	}
	if (strcmp(s,"rJ")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
		r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rK")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
		r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rL")==0)
	{
		r_pos[0] = rIx+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
		r_pos[1] = rIy-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rM")==0)
	{
		r_pos[0] = rIx+cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
		r_pos[1] = rIy-sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rN")==0)
	{
		r_pos[0] = rIx+cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
		r_pos[1] = rIy-sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rT")==0)
	{
		r_pos[0] = rIx+cos(theta-PI*(1.0/2.0))*(7.0/1.0E1);
		r_pos[1] = rIy-sin(theta-PI*(1.0/2.0))*(7.0/1.0E1);
	}
	if (strcmp(s,"rCOM")==0)
	{
		r_pos[0] = rIx-sin(q1+q2+theta)*1.678966789667897E-2-sin(q4+q5+theta)*1.678966789667897E-2-sin(q7+q8+theta)*8.717712177121771E-3-sin(q9+q10+theta)*8.717712177121771E-3-sin(q1+theta)*3.597785977859779E-2-sin(q4+theta)*3.597785977859779E-2-sin(q7+theta)*1.775830258302583E-2-sin(q9+theta)*1.775830258302583E-2+sin(theta)*2.506457564575646E-1-sqrt(2.0)*sin(q1+q2+q3+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(2.0)*sin(q4+q5+q6+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.690036900369004E-4;
		r_pos[1] = rIy-cos(q1+q2+theta)*1.678966789667897E-2-cos(q4+q5+theta)*1.678966789667897E-2-cos(q7+q8+theta)*8.717712177121771E-3-cos(q9+q10+theta)*8.717712177121771E-3-cos(q1+theta)*3.597785977859779E-2-cos(q4+theta)*3.597785977859779E-2-cos(q7+theta)*1.775830258302583E-2-cos(q9+theta)*1.775830258302583E-2+cos(theta)*2.506457564575646E-1-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(2.0)*cos(q1+q2+q3+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(2.0)*cos(q4+q5+q6+theta+PI*(1.0/4.0))*1.476014760147601E-3;
	}
	return r_pos;
}

std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s)
{
	double PI = 3.1415926535897932384626;
	std::vector<double> T(2);

	T[0] = 0.0;
	T[1] = 0.0;

	// double rIx = Robot_StateNDot_i.rIx;
	// double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;					double q1 = Robot_StateNDot_i.q1;				double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;						double q4 = Robot_StateNDot_i.q4;				double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;						double q7 = Robot_StateNDot_i.q7;				double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;						double q10 = Robot_StateNDot_i.q10;

	double rIxdot = Robot_StateNDot_i.rIxdot;				double rIydot = Robot_StateNDot_i.rIydot;		double thetadot = Robot_StateNDot_i.thetadot;
	double q1dot = Robot_StateNDot_i.q1dot;					double q2dot = Robot_StateNDot_i.q2dot;			double q3dot = Robot_StateNDot_i.q3dot;
	double q4dot = Robot_StateNDot_i.q4dot;					double q5dot = Robot_StateNDot_i.q5dot;			double q6dot = Robot_StateNDot_i.q6dot;
	double q7dot = Robot_StateNDot_i.q7dot;					double q8dot = Robot_StateNDot_i.q8dot;			double q9dot = Robot_StateNDot_i.q9dot;
	double q10dot = Robot_StateNDot_i.q10dot;

	if (strcmp(s,"vA")==0)
	{
		T[0] = rIxdot-q1dot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q2dot*(sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q3dot*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		T[1] = rIydot-q1dot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q2dot*(cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q3dot*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"vB")==0)
	{
		T[0] = rIxdot-q1dot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-thetadot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-q2dot*(sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q3dot*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
		T[1] = rIydot-q1dot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-thetadot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-q2dot*(cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q3dot*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"vC")==0)
	{
		T[0] = rIxdot-q4dot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q5dot*(sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q6dot*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		T[1] = rIydot-q4dot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q5dot*(cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q6dot*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"vD")==0)
	{
		T[0] = rIxdot-q4dot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-thetadot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-q5dot*(sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q6dot*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
		T[1] = rIydot-q4dot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-thetadot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-q5dot*(cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q6dot*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"vE")==0)
	{
		T[0] = rIxdot-q7dot*(sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q8dot*sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1);
		T[1] = rIydot-q7dot*(cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q8dot*cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1);
	}
	if(strcmp(s,"vF")==0)
	{
		T[0] = rIxdot-q9dot*(sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q10dot*sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1);
		T[1] = rIydot-q9dot*(cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q10dot*cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1);
	}
	if(strcmp(s,"vI")==0)
	{
		T[0] = rIxdot;
		T[1] = rIydot;
	}
	if(strcmp(s,"vM")==0)
	{
		T[0] = rIxdot-q7dot*sin(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(sin(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
		T[1] = rIydot-q7dot*cos(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(cos(q7+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
	}
	if(strcmp(s,"vN")==0)
	{
		T[0] = rIxdot-q9dot*sin(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(sin(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
		T[1] = rIydot-q9dot*cos(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)-thetadot*(cos(q9+theta+3.141592653589793*(1.0/2.0))*(1.0/4.0)+cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1));
	}
	if(strcmp(s,"vG")==0)
	{
		T[0] = rIxdot-q1dot*(sin(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(sin(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q2dot*sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
		T[1] = rIydot-q1dot*(cos(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(cos(q1+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q2dot*cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if(strcmp(s,"vJ")==0)
	{
		T[0] = rIxdot-q4dot*(sin(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(sin(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q5dot*sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
		T[1] = rIydot-q4dot*(cos(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-thetadot*(cos(q4+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1))-q5dot*cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if(strcmp(s,"vL")==0)
	{
		T[0] = rIxdot-thetadot*sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
		T[1] = rIydot-thetadot*cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if(strcmp(s,"vCOM")==0)
	{
		T[0] = rIxdot-q1dot*cos(q1+q2+q3+theta)*1.415929203539823E-3-q2dot*cos(q1+q2+q3+theta)*1.415929203539823E-3-q3dot*cos(q1+q2+q3+theta)*1.415929203539823E-3-q4dot*cos(q4+q5+q6+theta)*1.415929203539823E-3-q5dot*cos(q4+q5+q6+theta)*1.415929203539823E-3-q6dot*cos(q4+q5+q6+theta)*1.415929203539823E-3-thetadot*cos(q1+q2+q3+theta)*1.415929203539823E-3-thetadot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q1dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q2dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q3dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q4dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q5dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q6dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+thetadot*sin(q1+q2+q3+theta)*1.415929203539823E-3+thetadot*sin(q4+q5+q6+theta)*1.415929203539823E-3-q1dot*cos(q1+theta)*3.451327433628319E-2-q4dot*cos(q4+theta)*3.451327433628319E-2-q7dot*cos(q7+theta)*(1.9E1/9.04E2)-q9dot*cos(q9+theta)*(1.9E1/9.04E2)-thetadot*cos(q1+theta)*3.451327433628319E-2-thetadot*cos(q4+theta)*3.451327433628319E-2-thetadot*cos(q7+theta)*(1.9E1/9.04E2)-thetadot*cos(q9+theta)*(1.9E1/9.04E2)+thetadot*cos(theta)*2.628318584070796E-1-q10dot*cos(q9+q10+theta)*1.095132743362832E-2-q1dot*cos(q1+q2+theta)*1.610619469026549E-2-q2dot*cos(q1+q2+theta)*1.610619469026549E-2-q4dot*cos(q4+q5+theta)*1.610619469026549E-2-q5dot*cos(q4+q5+theta)*1.610619469026549E-2-q7dot*cos(q7+q8+theta)*1.095132743362832E-2-q8dot*cos(q7+q8+theta)*1.095132743362832E-2-q9dot*cos(q9+q10+theta)*1.095132743362832E-2-thetadot*cos(q1+q2+theta)*1.610619469026549E-2-thetadot*cos(q4+q5+theta)*1.610619469026549E-2-thetadot*cos(q7+q8+theta)*1.095132743362832E-2-thetadot*cos(q9+q10+theta)*1.095132743362832E-2-sqrt(4.1E1)*q1dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q4dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4-sqrt(4.1E1)*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4;
		T[1] = rIydot+q1dot*cos(q1+q2+q3+theta)*1.415929203539823E-3+q2dot*cos(q1+q2+q3+theta)*1.415929203539823E-3+q3dot*cos(q1+q2+q3+theta)*1.415929203539823E-3+q4dot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q5dot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q6dot*cos(q4+q5+q6+theta)*1.415929203539823E-3+thetadot*cos(q1+q2+q3+theta)*1.415929203539823E-3+thetadot*cos(q4+q5+q6+theta)*1.415929203539823E-3+q1dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q2dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q3dot*sin(q1+q2+q3+theta)*1.415929203539823E-3+q4dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q5dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q6dot*sin(q4+q5+q6+theta)*1.415929203539823E-3+thetadot*sin(q1+q2+q3+theta)*1.415929203539823E-3+thetadot*sin(q4+q5+q6+theta)*1.415929203539823E-3+q1dot*sin(q1+theta)*3.451327433628319E-2+q4dot*sin(q4+theta)*3.451327433628319E-2+q7dot*sin(q7+theta)*(1.9E1/9.04E2)+q9dot*sin(q9+theta)*(1.9E1/9.04E2)+thetadot*sin(q1+theta)*3.451327433628319E-2+thetadot*sin(q4+theta)*3.451327433628319E-2+thetadot*sin(q7+theta)*(1.9E1/9.04E2)+thetadot*sin(q9+theta)*(1.9E1/9.04E2)-thetadot*sin(theta)*2.628318584070796E-1+q10dot*sin(q9+q10+theta)*1.095132743362832E-2+q1dot*sin(q1+q2+theta)*1.610619469026549E-2+q2dot*sin(q1+q2+theta)*1.610619469026549E-2+q4dot*sin(q4+q5+theta)*1.610619469026549E-2+q5dot*sin(q4+q5+theta)*1.610619469026549E-2+q7dot*sin(q7+q8+theta)*1.095132743362832E-2+q8dot*sin(q7+q8+theta)*1.095132743362832E-2+q9dot*sin(q9+q10+theta)*1.095132743362832E-2+thetadot*sin(q1+q2+theta)*1.610619469026549E-2+thetadot*sin(q4+q5+theta)*1.610619469026549E-2+thetadot*sin(q7+q8+theta)*1.095132743362832E-2+thetadot*sin(q9+q10+theta)*1.095132743362832E-2+sqrt(4.1E1)*q1dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q4dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.539823008849558E-4+sqrt(4.1E1)*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.539823008849558E-4;
	}
	return T;
}
double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i)
{
	std::vector<double> Robotstate_vector = StateNDot2StateVec(Robot_StateNDot_i);
	dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full;
	Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
	dlib::matrix<double> Robotstatedot_Dlib, Robotstatedot_Dlib_Trans; Robotstatedot_Dlib = dlib::zeros_matrix<double>(13,1);
	for (int i = 0; i < 13; i++) {Robotstatedot_Dlib(i) = Robotstate_vector[13+i];}
		Robotstatedot_Dlib_Trans = dlib::trans(Robotstatedot_Dlib);
	double T;
	T = 0.5 * Robotstatedot_Dlib_Trans * D_q * Robotstatedot_Dlib;
	return T;
}
std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2)
{
	int Vec_Length = vec1.size();

	std::vector<double> vec_minus;

	for (int i = 0; i < Vec_Length; i++)
	{
		vec_minus.push_back(vec1[i] - vec2[i]);
	}
	return vec_minus;
}
void End_Effector_Obs_Dist_Fn(dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,6,1> &End_Effector_Dist, std::vector<int> &End_Effector_Obs)
{
	std::vector<double> r_Pos(2);		double Obs_Dist_i;			int Obs_Dist_Index;
	for (int i = 0; i < 6; i++){
		r_Pos[0] = End_Effector_Pos(2*i);
		r_Pos[1] = End_Effector_Pos(2*i+1);
		if(i<4)
		{
			Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "floor");
		}
		else
		{
			// Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "None");
			Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "floor");
			// Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index, "wall");
		}
		End_Effector_Dist(i) = Obs_Dist_i;
		End_Effector_Obs[i] = Obs_Dist_Index;}
		return;
	}
	void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,12,1> &End_Effector_Vel)
	{
		std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
		std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
		std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
		std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
		std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
		std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");

		End_Effector_Pos(0) = rA[0];
		End_Effector_Pos(1) = rA[1];
		End_Effector_Pos(2) = rB[0];
		End_Effector_Pos(3) = rB[1];
		End_Effector_Pos(4) = rC[0];
		End_Effector_Pos(5) = rC[1];
		End_Effector_Pos(6) = rD[0];
		End_Effector_Pos(7) = rD[1];
		End_Effector_Pos(8) = rE[0];
		End_Effector_Pos(9) = rE[1];
		End_Effector_Pos(10) = rF[0];
		End_Effector_Pos(11) = rF[1];

		std::vector<double> vA = Ang_Vel_fn(StateNDot_Init_i, "vA");
		std::vector<double> vB = Ang_Vel_fn(StateNDot_Init_i, "vB");
		std::vector<double> vC = Ang_Vel_fn(StateNDot_Init_i, "vC");
		std::vector<double> vD = Ang_Vel_fn(StateNDot_Init_i, "vD");
		std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");
		std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");

		End_Effector_Vel(0) = vA[0];
		End_Effector_Vel(1) = vA[1];
		End_Effector_Vel(2) = vB[0];
		End_Effector_Vel(3) = vB[1];
		End_Effector_Vel(4) = vC[0];
		End_Effector_Vel(5) = vC[1];
		End_Effector_Vel(6) = vD[0];
		End_Effector_Vel(7) = vD[1];
		End_Effector_Vel(8) = vE[0];
		End_Effector_Vel(9) = vE[1];
		End_Effector_Vel(10) = vF[0];
		End_Effector_Vel(11) = vF[1];
		return;
	}
	void Node_UpdateNCon(Tree_Node &Node_i, Robot_StateNDot &Node_StateNDot_i, std::vector<double> &sigma)
{	// This function is used to substitute the attribute and add it to the All_Nodes
	Node_i.Node_StateNDot = Node_StateNDot_i;
	Node_i.sigma = sigma;
	Node_i.KE = Kinetic_Energy_fn(Node_StateNDot_i);
	Node_i.Node_Index = All_Nodes.size();
	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
	End_Effector_PosNVel(Node_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
	Node_i.End_Effector_Pos = End_Effector_Pos;
	Node_i.End_Effector_Vel = End_Effector_Vel;
	All_Nodes.push_back(&Node_i);
	Frontier_Nodes.push_back(&Node_i);
	Frontier_Nodes_Cost.push_back(Node_i.KE);
}
void Nodes_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	// This function is the main optimization constraint function
	double T_tot, T;										dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;
	// Opt_Seed = Opt_Soln_Load();
	Opt_Seed_Unzip(Opt_Seed, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);
	T = T_tot/(Grids - 1)*1.0;
	Tree_Node Node_i, Node_i_child;							Node_i = Structure_P.Node_i;						Node_i_child = Structure_P.Node_i_child;
	std::vector<double> Robotstate_Initial_Vec = StateNDot2StateVec(Node_i.Node_StateNDot);
	int Opt_Type_Flag;										int Critical_Grid_Index;
	std::vector<double> sigma, sigma_i, sigma_i_child;		sigma_i = Node_i.sigma;								sigma_i_child = Node_i_child.sigma;
	Sigma_TransNGoal(sigma_i, sigma_i_child, Opt_Type_Flag, Critical_Grid_Index);

	double KE_ref = 0.001;			// A sufficiently trivial value

	// Objective function initialization
	ObjNConstraint_Val.push_back(0);
	ObjNConstraint_Type.push_back(1);

	// 1. The first constraint is to make sure that the initial condition matches the given initial condition
	dlib::matrix<double> StateNDot_Traj_1st_colm, Robotstate_Initial_VecDlib, Matrix_result;
	StateNDot_Traj_1st_colm = dlib::colm(StateNDot_Traj, 0);

	for (int i = 0; i < 26; i++)
	{
		ObjNConstraint_Val.push_back(StateNDot_Traj_1st_colm(i));
		ObjNConstraint_Type.push_back(1);
	}
	// THe constraint will be added term by term
	dlib::matrix<double> Robostate_Dlib_Front, Robostate_Dlib_Back, Ctrl_Front, Ctrl_Back, Contact_Force_Front, Contact_Force_Back, Robotstate_Mid_Acc;

	// 2. The second constraint is Dynamics constraints
	dlib::matrix<double> D_q_Front, B_q_Front, 	C_q_qdot_Front, Jac_Full_Front, Jac_Full_Trans_Front;
	dlib::matrix<double> D_q_Mid, 	B_q_Mid, 	C_q_qdot_Mid, 	Jac_Full_Mid, 	Jac_Full_Trans_Mid;
	dlib::matrix<double> D_q_Back, 	B_q_Back, 	C_q_qdot_Back, 	Jac_Full_Back, 	Jac_Full_Trans_Back;
	dlib::matrix<double> Dynamics_LHS, Dynamics_RHS;
	Robot_StateNDot Robot_StateNDot_Front, Robot_StateNDot_Mid, Robot_StateNDot_Back;
	for (int i = 0; i < Grids-1; i++)
	{
		// Get the robot state, ctrl, and contact force at the front and end of each segment
		Robostate_Dlib_Front = dlib::colm(StateNDot_Traj, i);						Robostate_Dlib_Back = dlib::colm(StateNDot_Traj, i+1);
		Ctrl_Front = dlib::colm(Ctrl_Traj,i);										Ctrl_Back = dlib::colm(Ctrl_Traj,i+1);
		Contact_Force_Front = dlib::colm(Contact_Force_Traj,i);						Contact_Force_Back = dlib::colm(Contact_Force_Traj,i+1);

		// Compute the Dynamics matrices at Front and Back
		Robot_StateNDot_Front = DlibRobotstate2StateNDot(Robostate_Dlib_Front);		Dynamics_Matrices(Robot_StateNDot_Front, D_q_Front, B_q_Front, C_q_qdot_Front, Jac_Full_Front);
		Robot_StateNDot_Back = DlibRobotstate2StateNDot(Robostate_Dlib_Back);		Dynamics_Matrices(Robot_StateNDot_Back, D_q_Back, B_q_Back, C_q_qdot_Back, Jac_Full_Back);

		Robot_StateNDot_MidNAcc(T, Robot_StateNDot_Front, Robot_StateNDot_Back, Ctrl_Front, Ctrl_Back, Contact_Force_Front, Contact_Force_Back, Robot_StateNDot_Mid, Robotstate_Mid_Acc, ObjNConstraint_Val, ObjNConstraint_Type);
		Dynamics_Matrices(Robot_StateNDot_Mid, D_q_Mid, B_q_Mid, C_q_qdot_Mid, Jac_Full_Mid);		Jac_Full_Trans_Mid = dlib::trans(Jac_Full_Mid);
		Dynamics_LHS = D_q_Mid * Robotstate_Mid_Acc + C_q_qdot_Mid;
		Dynamics_RHS = Jac_Full_Trans_Mid * (0.5 * Contact_Force_Front + 0.5 * Contact_Force_Back) + B_q_Mid * (0.5 * Ctrl_Front + 0.5 * Ctrl_Back);
		Matrix_result = Dynamics_LHS - Dynamics_RHS;
		// if((Opt_Type_Flag == -1)&&(i == Grids-2))
		// {
		// 	// In this case, the dynamics does not have to be satisfied
		// }
		// else
		// {
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
		// }
	}

	// 3. Complementarity constraints: Distance and Maintenance!
	//								   3.1 The previous unchanged active constraint have to be satisfied
	//								   3.2 Active constraints have to remain vacant relative distance.

	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;					dlib::matrix<double,6,1> End_Effector_Dist;
	dlib::matrix<double> Robostate_Dlib_i;											Robot_StateNDot Robot_StateNDot_i;
	std::vector<int> End_Effector_Obs(6);
	dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix;			std::vector<double> sigma_temp;
	dlib::matrix<double> End_Effector_Pos_ref;										End_Effector_Pos_ref = Node_i.End_Effector_Pos;
	std::vector<double> sigma_maint(12);											dlib::matrix<double> Matrix_Minus_result, End_Effector_PosDlib;

	//	3.1 The previous unchanged active constraint have to be satisfied
	for (int i = 0; i < Grids; i++)
	{
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);										Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
		if(i<Critical_Grid_Index){	sigma = sigma_i;}			else{			sigma = sigma_i_child;}
		if(Opt_Type_Flag == -1)
		{
			sigma_maint[0] = sigma[0];							sigma_maint[1] = sigma[0];
			sigma_maint[2] = sigma[0];							sigma_maint[3] = sigma[0];
			sigma_maint[4] = sigma[1];							sigma_maint[5] = sigma[1];
			sigma_maint[6] = sigma[1];							sigma_maint[7] = sigma[1];
			sigma_maint[8] = sigma[2];							sigma_maint[9] = sigma[2];
			sigma_maint[10] = sigma[3];							sigma_maint[11] = sigma[3];
		}
		else
		{
			sigma_maint[0] = sigma_i[0] * sigma_i_child[0];		sigma_maint[1] = sigma_i[0] * sigma_i_child[0];
			sigma_maint[2] = sigma_i[0] * sigma_i_child[0];		sigma_maint[3] = sigma_i[0] * sigma_i_child[0];
			sigma_maint[4] = sigma_i[1] * sigma_i_child[1];		sigma_maint[5] = sigma_i[1] * sigma_i_child[1];
			sigma_maint[6] = sigma_i[1] * sigma_i_child[1];		sigma_maint[7] = sigma_i[1] * sigma_i_child[1];
			sigma_maint[8] = sigma_i[2] * sigma_i_child[2];		sigma_maint[9] = sigma_i[2] * sigma_i_child[2];
			sigma_maint[10] = sigma_i[3] * sigma_i_child[3];	sigma_maint[11] = sigma_i[3] * sigma_i_child[3];
		}

		dlib::matrix<double> Maint_Matrix = Diag_Matrix_fn(sigma_maint);
		End_Effector_PosDlib = End_Effector_Pos;
		Matrix_Minus_result = End_Effector_PosDlib - End_Effector_Pos_ref;
		Matrix_result = Maint_Matrix * Matrix_Minus_result;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	}

	//	3.2 Active constraints have to remain vacant relative distance. However, when the motion type is to retract the contact, we have to deal with that differently.
	// The critical for the making/retracting contact is the same: the last grid (grids -1)

	for (int i = 0; i < Grids; i++)
	{
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);							Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
		End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);
		if(i<Critical_Grid_Index){	sigma = sigma_i;}			else{			sigma = sigma_i_child;}
		if(Opt_Type_Flag == -1)	{	sigma = sigma_i;}

		if((Opt_Type_Flag == -1)&&(i == Critical_Grid_Index-1))
		{
			dlib::matrix<double> Upper_Normal_Speed;
			End_Effector_Upper_Vel(Robot_StateNDot_i, Upper_Normal_Speed);

			ObjNConstraint_Val.push_back(sigma_i[0] * (!sigma_i_child[0]) * Robot_StateNDot_i.q2);
			ObjNConstraint_Type.push_back(0);
			ObjNConstraint_Val.push_back(sigma_i[1] * (!sigma_i_child[1]) * Robot_StateNDot_i.q5);
			ObjNConstraint_Type.push_back(0);
			ObjNConstraint_Val.push_back(sigma_i[2] * (!sigma_i_child[2]) * Robot_StateNDot_i.q8);
			ObjNConstraint_Type.push_back(0);
			ObjNConstraint_Val.push_back(sigma_i[3] * (!sigma_i_child[3]) * Robot_StateNDot_i.q10);
			ObjNConstraint_Type.push_back(0);

			// Maintain the COM to be within the SP
			std::vector<double> rCOM_i = Ang_Pos_fn(Robot_StateNDot_i, "rCOM");
			std::vector<double> rA_i = Ang_Pos_fn(Robot_StateNDot_i, "rA");
			std::vector<double> rB_i = Ang_Pos_fn(Robot_StateNDot_i, "rB");
			std::vector<double> rC_i = Ang_Pos_fn(Robot_StateNDot_i, "rC");
			std::vector<double> rD_i = Ang_Pos_fn(Robot_StateNDot_i, "rD");

			ObjNConstraint_Val.push_back(sigma_i[0] * (!sigma_i_child[0])*(rCOM_i[0] - rC_i[0]));
			ObjNConstraint_Type.push_back(1);
			ObjNConstraint_Val.push_back(sigma_i[0] * (!sigma_i_child[0])*(rD_i[0] - rCOM_i[0]));
			ObjNConstraint_Type.push_back(1);

			ObjNConstraint_Val.push_back(sigma_i[1] * (!sigma_i_child[1])*(rCOM_i[0] - rB_i[0]));
			ObjNConstraint_Type.push_back(1);
			ObjNConstraint_Val.push_back(sigma_i[1] * (!sigma_i_child[1])*(rA_i[0] - rCOM_i[0]));
			ObjNConstraint_Type.push_back(1);
		}
		else
		{
			sigma_temp = Sigma2Pos(sigma, 0);			Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
			sigma_temp = Sigma2Pos(sigma, 1);			Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
			sigma_temp = Sigma2Vel(sigma);				Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);

			// 3.1. Active constraints have to be satisfied
			Matrix_result = Eqn_Pos_Matrix * End_Effector_Dist;
			ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
			// cout<<"End Effector Dist:"<<endl;
			// cout<<End_Effector_Dist<<endl;
			// cout<<Eqn_Pos_Matrix<<endl;
			// cout<<Matrix_result<<endl;

			if((Opt_Type_Flag == 1)&&(i == Critical_Grid_Index))
			{
				// In this case, the velocity constraint at the end frame does not have to satisfied
				// Because an impact mapping will exist to shape the post-impact state
				double Impact_Mapping_Val = 1;
				// ObjNConstraint_Val.push_back(Robot_StateNDot_i.q5);
				// ObjNConstraint_Type.push_back(0);
			}
			else
			{
				if((Opt_Type_Flag == 1)&&(i == Critical_Grid_Index-1))
				{
					// ObjNConstraint_Val.push_back( Robot_StateNDot_i.q2);
					// ObjNConstraint_Type.push_back(0);
					// ObjNConstraint_Val.push_back((Robot_StateNDot_i.theta - PI/2.0) * (Robot_StateNDot_i.theta - PI/2.0));
					// ObjNConstraint_Type.push_back(0);
					// ObjNConstraint_Val.push_back(Robot_StateNDot_i.q8);
					// ObjNConstraint_Type.push_back(0);
					// ObjNConstraint_Val.push_back(Robot_StateNDot_i.q10);
					// ObjNConstraint_Type.push_back(0);
				}
				else
				{
					// Matrix_result = Eqn_Vel_Matrix * End_Effector_Vel;
					// ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
				}
			}
			// 3.2. Inactive constraints have to be strictly away from the obstacle
			dlib::matrix<double> ones_vector, temp_matrix;
			ones_vector = ONES_VECTOR_fn(6);
			Matrix_result = Ineqn_Pos_Matrix * (End_Effector_Dist - ones_vector * mini);
			ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
			// cout<<End_Effector_Dist<<endl;
			// cout<<Ineqn_Pos_Matrix<<endl;
			// cout<<Matrix_result<<endl;

			// 3.3 Middle joints have to be strictly away from the obs
			temp_matrix = Middle_Joint_Obs_Dist_Fn(Robot_StateNDot_i);
			ones_vector(5) = 10;
			Matrix_result = temp_matrix - ones_vector * mini;
			ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
		}
	}

	// 4. Complementarity Constraints: Contact Force!
	dlib::matrix<double> Contact_Force_Complem_Matrix, Contact_Force_i;
	for (int i = 0; i < Grids; i++) {
		if(i<Critical_Grid_Index){	sigma = sigma_i;}			else{			sigma = sigma_i_child;}
		std::vector<double> sigma_temp;				Contact_Force_i = dlib::colm(Contact_Force_Traj,i);
		sigma_temp.push_back(!sigma[0]);			sigma_temp.push_back(!sigma[0]);			sigma_temp.push_back(!sigma[0]);				sigma_temp.push_back(!sigma[0]);
		sigma_temp.push_back(!sigma[1]);			sigma_temp.push_back(!sigma[1]);			sigma_temp.push_back(!sigma[1]);				sigma_temp.push_back(!sigma[1]);
		sigma_temp.push_back(!sigma[2]);			sigma_temp.push_back(!sigma[2]);			sigma_temp.push_back(!sigma[3]);				sigma_temp.push_back(!sigma[3]);
		Contact_Force_Complem_Matrix = Diag_Matrix_fn(sigma_temp);
		Matrix_result = Contact_Force_Complem_Matrix * Contact_Force_i;
		// cout<<Contact_Force_Complem_Matrix<<endl;
		// cout<<Contact_Force_i<<endl;
		// cout<<Matrix_result<<endl;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	}

	// 5. Contact force feasibility constraints
	std::vector<double> Normal_Force_Front_vec(4), 	Tange_Force_Front_vec(4);
	std::vector<double> Normal_Force_Back_vec(4), 	Tange_Force_Back_vec(4);
	std::vector<double> Full_Normal_Front(6), 		Full_Normal_Back(6);

	for (int i = 0; i < Grids-1; i++) {
		// 5. Contact force feasibility constraints: 1. Normal force should be positive and
		//											 2. the friction cone constraint has to be satisfied
		//											 3. the horizontal force has to lie on the same side to avoid oscilation
		Contact_Force_Proj(StateNDot_Traj, Contact_Force_Traj, Normal_Force_Front_vec, Tange_Force_Front_vec, i, Full_Normal_Front);
		Contact_Force_Proj(StateNDot_Traj, Contact_Force_Traj, Normal_Force_Back_vec, Tange_Force_Back_vec, i+1, Full_Normal_Back);
		for (int k = 0; k < 4; k++)
		{
			ObjNConstraint_Val.push_back(Normal_Force_Front_vec[k] * Normal_Force_Front_vec[k] * mu * mu - Tange_Force_Front_vec[k] * Tange_Force_Front_vec[k]);
			ObjNConstraint_Type.push_back(1);
			ObjNConstraint_Val.push_back(Tange_Force_Front_vec[k] * Tange_Force_Back_vec[k]);
			ObjNConstraint_Type.push_back(1);
		}
		for (int mm = 0; mm < 6; mm++)
		{
			ObjNConstraint_Val.push_back(Full_Normal_Front[mm]);
			ObjNConstraint_Type.push_back(1);
		}
	}
	for (int k = 0; k < 4; k++)
	{
		ObjNConstraint_Val.push_back(Normal_Force_Back_vec[k] * Normal_Force_Back_vec[k] * mu * mu - Tange_Force_Back_vec[k] * Tange_Force_Back_vec[k]);
		ObjNConstraint_Type.push_back(1);
	}
	for (int mm = 0; mm < 6; mm++)
		{	ObjNConstraint_Val.push_back(Full_Normal_Back[mm]);
			ObjNConstraint_Type.push_back(1);}

			if(Opt_Type_Flag == 0)
			{
				dlib::matrix<double> Robostate_Dlib_j = dlib::colm(StateNDot_Traj, Grids-1);
				Robot_StateNDot Robot_StateNDot_j = DlibRobotstate2StateNDot(Robostate_Dlib_j);
				double KE_j = Kinetic_Energy_fn(Robot_StateNDot_j);
		// ObjNConstraint_Val[0] = KE_j;
		// ObjNConstraint_Val[0] = Robot_StateNDot_End.theta * Robot_StateNDot_End.theta;
		// ObjNConstraint_Val[0] = Objective_Function_Cal(StateNDot_Traj, Opt_Type_Flag, sigma_i_child);
		// ObjNConstraint_Val[0] = Robot_StateNDot_j.theta * Robot_StateNDot_j.theta;
				ObjNConstraint_Val[0] = Traj_Variation(StateNDot_Traj);
		// ObjNConstraint_Val[0] = KE_Variation_fn(StateNDot_Traj);
		// ObjNConstraint_Val[0] = 10 * Traj_Variation(StateNDot_Traj) + Variation_Cal(Contact_Force_Traj);
		// ObjNConstraint_Val[0] = Variation_Cal(Contact_Force_Traj);

		// double KE_End = Kinetic_Energy_End_Frame(StateNDot_Traj);
		// ObjNConstraint_Val.push_back(55 - KE_End);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(-45 + KE_End);
		// ObjNConstraint_Type.push_back(1);

				dlib::matrix<double> Robostate_Dlib_End = dlib::colm(StateNDot_Traj, Grids-1);
				Robot_StateNDot Robot_StateNDot_End = DlibRobotstate2StateNDot(Robostate_Dlib_End);

		// ObjNConstraint_Val[0] = (Robot_StateNDot_End.theta) * (Robot_StateNDot_End.theta);

				std::vector<double> rA_opt = Ang_Pos_fn(Robot_StateNDot_End, "rA");
				std::vector<double> rB_opt = Ang_Pos_fn(Robot_StateNDot_End, "rB");
				std::vector<double> rC_opt = Ang_Pos_fn(Robot_StateNDot_End, "rC");
				std::vector<double> rD_opt = Ang_Pos_fn(Robot_StateNDot_End, "rD");
				std::vector<double> rE_opt = Ang_Pos_fn(Robot_StateNDot_End, "rE");
				std::vector<double> rF_opt = Ang_Pos_fn(Robot_StateNDot_End, "rF");
				std::vector<double> rT_opt = Ang_Pos_fn(Robot_StateNDot_End, "rT");
				std::vector<double> vCOM_opt = Ang_Vel_fn(Robot_StateNDot_End, "vCOM");
		// ObjNConstraint_Val[0] = 0*(Robot_StateNDot_End.theta-PI/6) * (Robot_StateNDot_End.theta-PI/6);

		// ObjNConstraint_Val.push_back(vCOM_opt[1]);
		// ObjNConstraint_Type.push_back(1);

		// ObjNConstraint_Val.push_back((Robot_StateNDot_End.q1 + PI/2) * (Robot_StateNDot_End.q1 + PI/2));
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back((Robot_StateNDot_End.q2-PI/2) * (Robot_StateNDot_End.q2-PI/2));
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back(-Robot_StateNDot_End.q1 + Robot_StateNDot_End.q4 - mini);
		// ObjNConstraint_Type.push_back(1);

		//
		// ObjNConstraint_Val.push_back(rB_opt[0] - rC_opt[0] - mini);
		// ObjNConstraint_Type.push_back(1);


		// One more constraint to be added is the acceleration to be trivial value
				dlib::matrix<double> End_State_Traj = dlib::colm(StateNDot_Traj,Grids-1);			Robot_StateNDot End_StateNDot = DlibRobotstate2StateNDot(End_State_Traj);
				dlib::matrix<double> End_Ctrl_Traj = dlib::colm(Ctrl_Traj, Grids-1);				dlib::matrix<double> End_Contact_Force_Traj = dlib::colm(Contact_Force_Traj, Grids-1);

				dlib::matrix<double> D_q_End, B_q_End, C_q_qdot_End, Jac_Full_End, Ones_End_Dlib;
				Dynamics_Matrices(End_StateNDot, D_q_End, B_q_End, C_q_qdot_End, Jac_Full_End);

				dlib::matrix<double> Trivial_Acc_Matrix, Jac_Full_Trans_End;						Jac_Full_Trans_End = dlib::trans(Jac_Full_End);
				Trivial_Acc_Matrix = Jac_Full_Trans_End * End_Contact_Force_Traj + B_q_End * End_Ctrl_Traj - C_q_qdot_End;
				ObjNConstraint_ValNType_Update(Trivial_Acc_Matrix, ObjNConstraint_Val, ObjNConstraint_Type, 0);

				ObjNConstraint_Val.push_back(KE_ref - KE_j);
				ObjNConstraint_Type.push_back(1);

		// At the last frame

		// ObjNConstraint_Val.push_back((Robot_StateNDot_i.theta -PI/6.0) * (Robot_StateNDot_i.theta -PI/6.0));
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back(-Robot_StateNDot_i.q1 + Robot_StateNDot_i.q4 - mini);
		// ObjNConstraint_Type.push_back(1);
		//
		// ObjNConstraint_Val.push_back((Robot_StateNDot_i.q2 - PI/2.0) * (Robot_StateNDot_i.q2 - PI/2.0));
		// ObjNConstraint_Type.push_back(0);s
		//
		// ObjNConstraint_Val.push_back((Robot_StateNDot_i.q5) * (Robot_StateNDot_i.q5));
		// ObjNConstraint_Type.push_back(0);
			}
			else
			{
		// ObjNConstraint_Val[0] = -T;
		// ObjNConstraint_Val[0] = Objective_Function_Cal(StateNDot_Traj, Opt_Type_Flag, sigma_i_child);
				double KE_End = Kinetic_Energy_End_Frame(StateNDot_Traj);
		// ObjNConstraint_Val.push_back(60 - KE_End);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(KE_End - 90.66);
		// ObjNConstraint_Type.push_back(1);
		//
		// ObjNConstraint_Val.push_back(35.09 - KE_End);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(KE_End - 30.08);
		// ObjNConstraint_Type.push_back(1);


		// ObjNConstraint_Val[0] = KE_End;

		// ObjNConstraint_Val[0] = (KE_End - 90) * (KE_End - 90);

		// ObjNConstraint_Val[0] = (KE_End - 55) * (KE_End - 55);

				dlib::matrix<double> Robostate_Dlib_End = dlib::colm(StateNDot_Traj, Grids-1);
				Robot_StateNDot Robot_StateNDot_End = DlibRobotstate2StateNDot(Robostate_Dlib_End);

		// ObjNConstraint_Val[0] = (Robot_StateNDot_End.theta) * (Robot_StateNDot_End.theta);

				std::vector<double> rA_opt = Ang_Pos_fn(Robot_StateNDot_End, "rA");
				std::vector<double> rB_opt = Ang_Pos_fn(Robot_StateNDot_End, "rB");
				std::vector<double> rC_opt = Ang_Pos_fn(Robot_StateNDot_End, "rC");
				std::vector<double> rD_opt = Ang_Pos_fn(Robot_StateNDot_End, "rD");
		// std::vector<double> rE_opt = Ang_Pos_fn(Robot_StateNDot_End, "rE");
		// std::vector<double> rF_opt = Ang_Pos_fn(Robot_StateNDot_End, "rF");
		// std::vector<double> rT_opt = Ang_Pos_fn(Robot_StateNDot_End, "rT");
		// std::vector<double> vCOM_opt = Ang_Vel_fn(Robot_StateNDot_End, "vCOM");

		// ObjNConstraint_Val[0] = Robot_StateNDot_End.thetadot * Robot_StateNDot_End.thetadot;
		// ObjNConstraint_Val.push_back(-Robot_StateNDot_End.q9 - mini);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(rD_opt[0] - rA_opt[0] - 0.10);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(0.2 - (rD_opt[0] - rA_opt[0]));
		// ObjNConstraint_Type.push_back(1);

		// ObjNConstraint_Val.push_back((Robot_StateNDot_End.q10 + PI/2.0) * (Robot_StateNDot_End.q10 + PI/2.0));
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back((Robot_StateNDot_End.q2-Robot_StateNDot_End.q5) * (Robot_StateNDot_End.q2-Robot_StateNDot_End.q5));
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back(rB_opt[0] - rC_opt[0] - 2 * mini);
		// ObjNConstraint_Type.push_back(1);

		// ObjNConstraint_Val[0] = 0*Traj_Variation(StateNDot_Traj) + 0.0 * Variation_Cal(Contact_Force_Traj);
				ObjNConstraint_Val[0] = Traj_Variation(StateNDot_Traj) + 0 * Quadratic_Sum(Contact_Force_Traj);
		// ObjNConstraint_Val[0] = KE_Variation_fn(StateNDot_Traj);
		// ObjNConstraint_Val[0] = Variation_Cal(Contact_Force_Traj);
		// ObjNConstraint_Val[0] = Objective_Function_Cal(StateNDot_Traj, Opt_Type_Flag, sigma_i_child) + 0 * Variation_Cal(Contact_Force_Traj);
			}
		}

		double Objective_Function_Cal(dlib::matrix<double> &StateNDot_Traj, int Opt_Type_Flag, std::vector<double> &sigma_i_child)
		{
	// 7. Objective function value: the first value is fixed
			double KE_i;					double KE_tot = 0.0;
			double Impulse_Mag = 0.0;		double Obj_val = 0.0;
	// double Obj_Val = 0.0;
	// // for (int i = 0; i < 13; i++) {
	// // 	Obj_Val = Obj_Val + End_Robotstate(i+13) * End_Robotstate(i+13);
	// // }
	//
	// for (int i = 0; i < Grids; i++) {
		// 7. Objective function value: the first value is fixed
	// 	dlib::matrix<double> Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);
	// 	Robot_StateNDot Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
	// 	KE_i = Kinetic_Energy_fn(Robot_StateNDot_i);
	// 	KE_tot = KE_tot + KE_i;
	// }
			if(Opt_Type_Flag == 1)
			{
		// In this case, there is impact mapping involved
				Robot_StateNDot Robot_StateNDot_End_Grid;
				Robot_StateNDot_End_Grid = Impact_Mapping_fn(StateNDot_Traj, Impulse_Mag, sigma_i_child);
			}
			Obj_val = Kinetic_Energy_End_Frame(StateNDot_Traj) + 10 * Impulse_Mag;
			return Obj_val;
		}

		double Kinetic_Energy_End_Frame(dlib::matrix<double> &StateNDot_Traj)
		{
			dlib::matrix<double> End_Robotstate;
			End_Robotstate = dlib::colm(StateNDot_Traj, Grids-1);
			Robot_StateNDot Robot_StateNDot_i = DlibRobotstate2StateNDot(End_Robotstate);
			double KE_i = Kinetic_Energy_fn(Robot_StateNDot_i);
			return KE_i;
		}

		Robot_StateNDot Impact_Mapping_fn(dlib::matrix<double> &StateNDot_Traj, double &Impulse_Mag, std::vector<double> &sigma_i_child)
		{
			dlib::matrix<double> End_State_Traj = dlib::colm(StateNDot_Traj,Grids-1);
			Robot_StateNDot End_StateNDot = DlibRobotstate2StateNDot(End_State_Traj);
			std::vector<double> End_Statevec = StateNDot2StateVec(End_StateNDot);
			dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full;
			Dynamics_Matrices(End_StateNDot, D_q, B_q, C_q_qdot, Jac_Full);

	// Here sigma_i_child is used to select the full row rank Jacobian matrix
			std::vector<double> Jac_Act_Index = Full_Row_Rank_Index(sigma_i_child);
			const int Jac_Row_Number = Jac_Act_Index.size();
			dlib::matrix<double> Jac_Act;					Jac_Act = dlib::zeros_matrix<double>(Jac_Row_Number,13);
			for (int j = 0; j < Jac_Row_Number; j++) {
				for (int k = 0; k < 13; k++) {
					Jac_Act(j,k) = Jac_Full(Jac_Act_Index[j],k);
				}
			}
			dlib::matrix<double> Jac_Act_Trans = dlib::trans(Jac_Act);
			dlib::matrix<double> D_q_Inv = dlib::inv(D_q);
			dlib::matrix<double> Pre_Impact_Vel, Post_Impact_Vel;
			Pre_Impact_Vel = dlib::zeros_matrix<double>(13,1);
			Post_Impact_Vel = Pre_Impact_Vel;
			for (int i = 0; i < 13; i++)
			{
				Pre_Impact_Vel(i) = End_State_Traj(i+13);
			}
			dlib::matrix<double> Impulse_Lamda;
	// cout<<D_q<<endl;			cout<<B_q<<endl;			cout<<C_q_qdot<<endl;			cout<<Jac_Act<<endl;
			Impulse_Lamda = -dlib::pinv(Jac_Act * D_q_Inv * Jac_Act_Trans) * Jac_Act * Pre_Impact_Vel;
	// cout<<Impulse_Lamda<<endl;
			Post_Impact_Vel = D_q_Inv * Jac_Act_Trans * Impulse_Lamda + Pre_Impact_Vel;
	// cout<<Post_Impact_Vel<<endl;
			std::vector<double> Robotstate_End(26);
			for (int i = 0; i < 13; i++) {
				Robotstate_End[i] = End_State_Traj(i);
				Robotstate_End[i+13] = Post_Impact_Vel(i);
			}
			Robot_StateNDot Robot_StateNDot_End;
			Robot_StateNDot_End = StateVec2StateNDot(Robotstate_End);


			Impulse_Mag = 0.0;
			for (int i = 0; i < Impulse_Lamda.nr(); i++) {
				Impulse_Mag = Impulse_Mag + Impulse_Lamda(i) * Impulse_Lamda(i);
			}
			return Robot_StateNDot_End;
		}
		std::vector<double> Full_Row_Rank_Index(std::vector<double> &sigma_i_child)
		{
			std::vector<double> Row_Rank_Index_Vec;
			if(sigma_i_child[0]==1)
			{
				Row_Rank_Index_Vec.push_back(0);
				Row_Rank_Index_Vec.push_back(1);
				Row_Rank_Index_Vec.push_back(3);
			}
			if(sigma_i_child[1]==1)
			{
				Row_Rank_Index_Vec.push_back(4);
				Row_Rank_Index_Vec.push_back(5);
				Row_Rank_Index_Vec.push_back(7);
			}
			if(sigma_i_child[2]==1)
			{
				Row_Rank_Index_Vec.push_back(8);
				Row_Rank_Index_Vec.push_back(9);
			}
			if(sigma_i_child[3]==1)
			{
				Row_Rank_Index_Vec.push_back(10);
				Row_Rank_Index_Vec.push_back(11);
			}
			return Row_Rank_Index_Vec;
		}
		void Contact_Force_Proj(dlib::matrix<double> &StateNDot_Traj, dlib::matrix<double> &Contact_Force_Traj, std::vector<double> &Normal_Force_vec, std::vector<double> &Tange_Force_vec, int Grid_Index, std::vector<double> &Full_Normal)
		{
	// This function is used to project the contact force in to the normal and tangential direction
			Robot_StateNDot Robot_StateNDot_i;
			dlib::matrix<double> Robostate_Dlib_i = dlib::colm(StateNDot_Traj, Grid_Index);
			Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);

			dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;							dlib::matrix<double,6,1> End_Effector_Dist;
			std::vector<int> End_Effector_Obs(6);

			double Contact_Force_i_x, Contact_Force_i_y;											std:vector<double> Normal_Force, Tange_Force;
			double Normal_Force_1, Normal_Force_2, Normal_Force_3, Normal_Force_4;					double Tange_Force_1, Tange_Force_2, Tange_Force_3, Tange_Force_4;

			End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);			End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);
			dlib::matrix<double> Contact_Force_i = dlib::colm(Contact_Force_Traj,Grid_Index);

			for (int j = 0; j < Contact_Force_i.nr()/2; j++)
			{
				Contact_Force_i_x = Contact_Force_i(2*j);											Contact_Force_i_y = Contact_Force_i(2*j+1);
				Normal_Force.push_back(Contact_Force_i_x * Envi_Map_Normal(End_Effector_Obs[j],0) + Contact_Force_i_y * Envi_Map_Normal(End_Effector_Obs[j],1));
				Tange_Force.push_back(Contact_Force_i_x * Envi_Map_Tange(End_Effector_Obs[j],0) + Contact_Force_i_y * Envi_Map_Tange(End_Effector_Obs[j],1));
			}

			Normal_Force_1 = Normal_Force[0] + Normal_Force[1];		Normal_Force_2 = Normal_Force[2] + Normal_Force[3];		Normal_Force_3 = Normal_Force[4];	Normal_Force_4 = Normal_Force[5];
			Tange_Force_1 = Tange_Force[0] + Tange_Force[1];		Tange_Force_2 = Tange_Force[2] + Tange_Force[3];		Tange_Force_3 = Tange_Force[4];		Tange_Force_4 = Tange_Force[5];

			Normal_Force_vec[0] = Normal_Force_1;				Normal_Force_vec[1] = Normal_Force_2;				Normal_Force_vec[2] = Normal_Force_3;				Normal_Force_vec[3] = Normal_Force_4;
			Tange_Force_vec[0] = Tange_Force_1;					Tange_Force_vec[1] = Tange_Force_2;					Tange_Force_vec[2] = Tange_Force_3;					Tange_Force_vec[3] = Tange_Force_4;

	// Full normal
			Full_Normal[0] =  Normal_Force[0];				Full_Normal[1] =  Normal_Force[1];				Full_Normal[2] =  Normal_Force[2];
			Full_Normal[3] =  Normal_Force[3];				Full_Normal[4] =  Normal_Force[4];				Full_Normal[5] =  Normal_Force[5];

		}
		void End_Effector_Upper_Vel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double> &End_Effector_Upper_Normal_Speed)
		{
	// This function is used to calculate the velocity of the upper joint of the end effector
	// Along the normal direction of its link
	// For the foot link, we only consider its contact between the floor
			End_Effector_Upper_Normal_Speed = dlib::zeros_matrix<double>(4,1);

			std::vector<double> vG = Ang_Vel_fn(StateNDot_Init_i, "vG");
			std::vector<double> vJ = Ang_Vel_fn(StateNDot_Init_i, "vJ");
			std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");
			std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");
			std::vector<double> vI = Ang_Vel_fn(StateNDot_Init_i, "vI");
			std::vector<double> vL = Ang_Vel_fn(StateNDot_Init_i, "vL");

			std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
			std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
			std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
			std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
			std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
			std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
			std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
			std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");

			End_Effector_Upper_Normal_Speed(0) = Velocity_Projection(rG, rI, vI);
			End_Effector_Upper_Normal_Speed(1) = Velocity_Projection(rJ, rI, vI);
			End_Effector_Upper_Normal_Speed(2) = Velocity_Projection(rE, rL, vL);
			End_Effector_Upper_Normal_Speed(3) = Velocity_Projection(rF, rL, vL);

		}
		double Velocity_Projection(std::vector<double> &Pos_A, std::vector<double> &Pos_B, std::vector<double> &Vel_B)
		{
	// This function is used to project the upper velocity along the normal direction
			double End_Effector_x, End_Effector_y, End_Effector_Upper_x, End_Effector_Upper_y, Slope_Angle;
			End_Effector_x = Pos_A[0];					End_Effector_y = Pos_A[1];
			End_Effector_Upper_x = Pos_B[0];			End_Effector_Upper_y = Pos_B[1];
			Slope_Angle = atan2(End_Effector_Upper_y - End_Effector_y, End_Effector_Upper_x - End_Effector_x);
			double Lamda_x = cos(Slope_Angle);			double Lamda_y = sin(Slope_Angle);
			double Projected_Vel = Vel_B[0] * Lamda_x + Vel_B[1] * Lamda_y;
			return Projected_Vel;
		}

		double Variation_Cal(dlib::matrix<double> &Dlib_Matrix)
		{
	// This function is used to calcualte the variation of the state trajectory
			dlib::matrix<double> StateNDot_Traj_Front, StateNDot_Traj_Back, Matrix_result;
			double Traj_Variation_Val = 0.0;
			for (int i = 0; i < Dlib_Matrix.nc()-1; i++) {
				StateNDot_Traj_Front = dlib::colm(Dlib_Matrix, i);
				StateNDot_Traj_Back = dlib::colm(Dlib_Matrix, i+1);
				Matrix_result = StateNDot_Traj_Front - StateNDot_Traj_Back;
				for (int j = 0; j < Matrix_result.nr(); j++) {
					Traj_Variation_Val = Traj_Variation_Val +  Matrix_result(j) * Matrix_result(j);
				}
			}
	// Traj_Variation_Val = Traj_Variation_Val * Traj_Variation_Val;
			return Traj_Variation_Val;
		}

		double Quadratic_Sum(dlib::matrix<double> &Dlib_Matrix)
		{
	// This function is used to calcualte the variation of the state trajectory
			double Quadratic_Sum_Val = 0.0;
			for (int i = 0; i < Dlib_Matrix.nc(); i++)
			{
				for (int j = 0; j < Dlib_Matrix.nr(); j++)
				{
					Quadratic_Sum_Val = Quadratic_Sum_Val + Dlib_Matrix(i,j) * Dlib_Matrix(i,j);
				}
			}
			return Quadratic_Sum_Val;
		}

		double Traj_Variation(dlib::matrix<double> &StateNDot_Traj)
		{
	// This function is used to calcualte the variation of the state trajectory
			dlib::matrix<double> StateNDot_Traj_Front, StateNDot_Traj_Back, Matrix_result;
			double Traj_Variation_Val = 0.0;
			for (int i = 0; i < StateNDot_Traj.nc()-1; i++) {
				StateNDot_Traj_Front = dlib::colm(StateNDot_Traj, i);
				StateNDot_Traj_Back = dlib::colm(StateNDot_Traj, i+1);
				Matrix_result = StateNDot_Traj_Front - StateNDot_Traj_Back;
				for (int j = 0; j < Matrix_result.nr()/2; j++) {
					Traj_Variation_Val = Traj_Variation_Val +  Matrix_result(j) * Matrix_result(j);
				}
			}
	// Traj_Variation_Val = Traj_Variation_Val * Traj_Variation_Val;
			return Traj_Variation_Val;
		}
		Robot_StateNDot DlibRobotstate2StateNDot(dlib::matrix<double> &DlibRobotstate)
		{
	// This function is used to convert the dlib matrix robot state to Robot_StateNDot type
			std::vector<double> Robot_StateNDot_vec;
			for (int i = 0; i < DlibRobotstate.nr(); i++) {
				Robot_StateNDot_vec.push_back(DlibRobotstate(i));}
				Robot_StateNDot Robot_StateNDot_i(Robot_StateNDot_vec);
				return Robot_StateNDot_i;
			}
			void Robot_StateNDot_MidNAcc(double T, const Robot_StateNDot &Robot_StateNDot_Front, const Robot_StateNDot &Robot_StateNDot_Back, const dlib::matrix<double> &Ctrl_Front, const dlib::matrix<double> &Ctrl_Back, const dlib::matrix<double> &Contact_Force_Front, const dlib::matrix<double> &Contact_Force_Back, Robot_StateNDot &Robot_StateNDot_Mid, dlib::matrix<double> &Robotstate_Mid_Acc,std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
			{
				std::vector<double> Robotstate_Vec_Front, Robotstate_Vec_Back;
				dlib::matrix<double> D_q_Front, B_q_Front, 	C_q_qdot_Front, Jac_Full_Front, Jac_Full_Trans_Front, Acc_Front;
				dlib::matrix<double> D_q_Back, 	B_q_Back, 	C_q_qdot_Back, 	Jac_Full_Back, 	Jac_Full_Trans_Back, Acc_Back;

				Dynamics_Matrices(Robot_StateNDot_Front, D_q_Front, B_q_Front, C_q_qdot_Front, Jac_Full_Front);
				Jac_Full_Trans_Front = dlib::trans(Jac_Full_Front);
				Dynamics_Matrices(Robot_StateNDot_Back, D_q_Back, B_q_Back, C_q_qdot_Back, Jac_Full_Back);
				Jac_Full_Trans_Back = dlib::trans(Jac_Full_Back);

	// cout<<D_q_Front<<endl;				cout<<B_q_Front<<endl;				cout<<C_q_qdot_Front<<endl;			cout<<Jac_Full_Front<<endl;

				Acc_Front = dlib::inv(D_q_Front) * (Jac_Full_Trans_Front * Contact_Force_Front + B_q_Front * Ctrl_Front - C_q_qdot_Front);
				Acc_Back =  dlib::inv(D_q_Back) *  (Jac_Full_Trans_Back *  Contact_Force_Back +  B_q_Back *  Ctrl_Back -  C_q_qdot_Back);

	// For each variable in the state, we will calculate the Cubic spline coefficients and the middle state N Acc
				std::vector<double> Robotstate_Vec_Mid(26);										Robotstate_Mid_Acc = dlib::zeros_matrix<double>(13,1);
				Robotstate_Vec_Front = StateNDot2StateVec(Robot_StateNDot_Front);
				Robotstate_Vec_Back = StateNDot2StateVec(Robot_StateNDot_Back);

				double x_init, x_end, xdot_init, xdot_end, xddot_init, xddot_end;				std::vector<double> CubicSpline_Coeff;

	// 1. Calculate the Robotstate_Mid_Acc Pos
				for (int i = 0; i < 13; i++) {
					x_init = Robotstate_Vec_Front[i];				x_end = Robotstate_Vec_Back[i];
					xdot_init = Robotstate_Vec_Front[i+13];			xdot_end = Robotstate_Vec_Back[i+13];
					CubicSpline_Coeff = CubicSpline_Coeff_fn(T, x_init, x_end, xdot_init, xdot_end);
					Robotstate_Vec_Mid[i] = CubicSpline_Evaluation_fn(CubicSpline_Coeff, 0.5);
				}
	// 2. Calculate the Robotstate_Mid_Acc Vel and Acc
				for (int i = 0; i < 13; i++) {
					xdot_init = Robotstate_Vec_Front[i+13];			xdot_end = Robotstate_Vec_Back[i+13];
					xddot_init = Acc_Front(i);						xddot_end = Acc_Back(i);
					CubicSpline_Coeff = CubicSpline_Coeff_fn(T, xdot_init, xdot_end, xddot_init, xddot_end);

					double xdot_max = max(xdot_init, xdot_end);
					double xdot_min = min(xdot_init, xdot_end);
					double Pos_Real_Change = Robotstate_Vec_Back[i] - Robotstate_Vec_Front[i];
					ObjNConstraint_Val.push_back(xdot_max*T - Pos_Real_Change);
					ObjNConstraint_Type.push_back(1);
					ObjNConstraint_Val.push_back(Pos_Real_Change - xdot_min*T);
					ObjNConstraint_Type.push_back(1);

					Robotstate_Vec_Mid[i+13] = CubicSpline_Evaluation_fn(CubicSpline_Coeff, 0.5);
					Robotstate_Mid_Acc(i) = CubicSpline_1stOrder_Evaluation_fn(CubicSpline_Coeff, 0.5, T);
				}
	// cout<<Acc_Front<<endl;			cout<<Acc_Back<<endl;			cout<<Robotstate_Mid_Acc<<endl;
	// The acceleration constraint is not of first priority
				Robot_StateNDot_Mid = StateVec2StateNDot(Robotstate_Vec_Mid);
				for (int i = -1; i < 10; i++) {
					ObjNConstraint_Val.push_back(Acc_Front(i+3));
					ObjNConstraint_Type.push_back(2);
					ObjNConstraint_Val.push_back(Robotstate_Mid_Acc(i+3));
					ObjNConstraint_Type.push_back(2);
					ObjNConstraint_Val.push_back(Acc_Back(i+3));
					ObjNConstraint_Type.push_back(2);
				}
			}
			void Quadratic_Angular_Sum_Cal(std::vector<double> &Robot_Vel,double &Quadratic_Angular_Sum)
			{
				for (int i = 0; i < Robot_Vel.size(); i++) {
					Quadratic_Angular_Sum = Quadratic_Angular_Sum + Robot_Vel[i] * Robot_Vel[i];
				}
			}
			dlib::matrix<double> Quadratic_Minus(dlib::matrix<double> &Mat_A, dlib::matrix<double> &Mat_B)
			{
	// This function is used to take care of the situation where the SNOPT cannot be easily used to optimize the simple minus constraint
				dlib::matrix<double> Matrix_Minus, Matrix_result; Matrix_result = Mat_A;
				Matrix_Minus = Mat_A - Mat_B;
				for (int i = 0; i < Matrix_Minus.nr(); i++) {
					Matrix_result(i) = Matrix_Minus(i) * Matrix_Minus(i);
				}
				return Matrix_result;
			}
			std::vector<double> Opt_Soln_Load()
			{
	// This function is used to load the computed optimal solution for data analysis
				std::vector<double> Opt_Seed;
	ifstream Opt_Soln_File;              // This is to read the initial angle and angular velocities
	Opt_Soln_File.open("4.txt");
	if(Opt_Soln_File.is_open())
	{
		double data_each_line = 0.0;

		while(Opt_Soln_File>>data_each_line)
		{
			Opt_Seed.push_back(data_each_line);
		}
		Opt_Soln_File.close();
	}
	return Opt_Seed;
}
double KE_Variation_fn(dlib::matrix<double> &StateNDot_Traj)
{	double KE_Variation = 0.0;
	std::vector<double> KE_tot;
	dlib::matrix<double> Robostate_Dlib_i;
	Robot_StateNDot Robot_StateNDot_i;
	for (int i = 0; i < StateNDot_Traj.nc(); i++) {
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);
		Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		KE_tot.push_back(Kinetic_Energy_fn(Robot_StateNDot_i));
	}
	for (int i = 0; i < KE_tot.size(); i++)
	{
		KE_Variation = KE_Variation + KE_tot[i];
	}
	return KE_Variation;
}
void Sigma_TransNGoal(std::vector<double> & sigma_i, std::vector<double> & sigma_i_child, int &Opt_Type_Flag, int &Crit_Grid)
{
	// The role of this function is actually very important because making contact is easier than retracting contact since retracting contact requires an accumulation of
	// the necessary kinetic energy of certain body parts before a retract can be made
	double sigma_result = 0.0;
	Opt_Type_Flag = 0;
	for (int i = 0; i < sigma_i.size(); i++) {
		sigma_result = sigma_result + sigma_i_child[i] - sigma_i[i];
	}
	if(sigma_result>0)
	{
		// In this case, it is making the contact so the Critical frame is the ending frame
		Opt_Type_Flag = 1;
		Crit_Grid = Grids-1;
	}
	else
	{
		if (sigma_result ==0)
		{
			// In this case, it is the self-opt process
			Opt_Type_Flag = 0;
			Crit_Grid = Grids;
		}
		else
		{
			// In this case, it is the retracting contact process: the whole process is divided into two subprocesses, energy accumulation, energy release
			Opt_Type_Flag = -1;
			Crit_Grid = Grids-1;
		}
	}
}

dlib::matrix<double> StateVec2DlibMatrix_fn(const std::vector<double> &StateVec)
{	const int dim = StateVec.size();
	dlib::matrix<double> StateVec2DlibMatrix;
	StateVec2DlibMatrix = dlib::zeros_matrix<double>(dim,1);
	for (int i = 0; i < dim; i++) {
		StateVec2DlibMatrix(i) = StateVec[i];}
		return StateVec2DlibMatrix;
	}
	void CtrlNContact_ForcefromCtrlNContact_Force_Coeff(dlib::matrix<double> &Ctrl_Coeff,dlib::matrix<double> &Contact_Force_Coeff, int Grid_Ind, double s, dlib::matrix<double> &Ctrl_i,  dlib::matrix<double> &Contact_Force_i)
	{	double x_a, x_b;
		Ctrl_i = dlib::zeros_matrix<double>(10,1);				Contact_Force_i = dlib::zeros_matrix<double>(12,1);
		for (int i = 0; i < 10; i++) {
			x_a = Ctrl_Coeff(2*i, Grid_Ind);
			x_b = Ctrl_Coeff(2*i+1, Grid_Ind);
			Ctrl_i(i) = x_a * s + x_b;}
			for (int i = 0; i < 12; i++) {
				x_a = Contact_Force_Coeff(2*i, Grid_Ind);
				x_b = Contact_Force_Coeff(2*i+1, Grid_Ind);
				Contact_Force_i(i) = x_a * s + x_b;}
			}
			int Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double> &Opt_Soln_Output)
{	// This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
	// However, the constraint will be set to use the direct collocation method
	int Opt_Flag = 0;
	Structure_P.Node_i = Node_i;		Structure_P.Node_i_child = Node_i_child;
	double Time_Interval = 0.02;		int Total_Num = 41;
	std::vector<double> Time_Queue = Time_Seed_Queue_fn(Time_Interval, Total_Num);
	std::vector<double> Opt_Soln_Tot;
	int Feasible_Num = 0;

	for (int j = 0; j < Time_Queue.size(); j++) {
		std::vector<double> Opt_Soln, ObjNConstraint_Val, ObjNConstraint_Type;
		Time_Seed = Time_Queue[j];								// I forget what is the meaning of my Time_Seed but I would say that is the duration between two knots
		cout<<"--------------------------------- Time Guess Iteration: "<<j<<" --------------------------------- "<<endl;
		cout<<"The Guess for the time duaration is :" <<Time_Seed<<" s"<<endl;
		Opt_Soln = Nodes_Optimization_Inner_Opt(Node_i, Node_i_child);
		Nodes_Optimization_ObjNConstraint(Opt_Soln, ObjNConstraint_Val, ObjNConstraint_Type);
		double ObjNConstraint_Violation_Val = ObjNConstraint_Violation(ObjNConstraint_Val, ObjNConstraint_Type);
		cout<<endl;
		cout<<"-----------------------------ObjNConstraint_Violation_Val: "<<ObjNConstraint_Violation_Val<<"-------------------------------"<<endl;
		if(ObjNConstraint_Violation_Val<0.1)
		{
			Opt_Flag = 1;
			Feasible_Num = Feasible_Num + 1;
			// However, due to the small constraint violation, the robot may not 100% satisfy the holonomic constraint so here we would like to check the robot final state
			// Opt_Soln =  Final_State_Opt(Opt_Soln, Node_i, Node_i_child);
			for (int i = 0; i < Opt_Soln.size(); i++) {
				Opt_Soln_Tot.push_back(Opt_Soln[i]);}
				time_t now = time(0);
				// convert now to string form
				char* dt = ctime(&now);

				ofstream output_file;
				std::string pre_filename = "From_Node_";
				std::string Node_i_name = to_string(Node_i.Node_Index);
				std::string mid_filename = "_Expansion_At_Feasible_";
				std::string Node_i_Child_name = to_string(Feasible_Num);
				std::string post_filename = ".txt";

				std::string filename = pre_filename + Node_i_name + mid_filename + Node_i_Child_name + dt + post_filename;
				output_file.open(filename, std::ofstream::out);
				for (int i = 0; i < Opt_Soln.size(); i++)
				{
					output_file<<Opt_Soln[i]<<endl;
				}
				output_file.close();
			}
			else
			{	// The differentiation between the pure seed guess and the failure guess needs to be conducted
				// Then this is a failure attempt
				if(ObjNConstraint_Violation_Val<100)
				{
					time_t now = time(0);
					// convert now to string form
					char* dt = ctime(&now);
					ofstream output_file;
					std::string pre_filename = "Node_";
					std::string Node_i_name = to_string(Node_i.Node_Index);
					std::string mid_filename = "_Expansion_Quasi_At_";
					std::string post_filename = ".txt";

					std::string filename = pre_filename + Node_i_name + mid_filename + dt + post_filename;
					output_file.open(filename, std::ofstream::out);
					for (int i = 0; i < Opt_Soln.size(); i++)
						{	output_file<<Opt_Soln[i]<<endl;}
					output_file.close();
				}
			}
		}
		if(Opt_Flag == 1)
		{
			// Now it is time to select the best one from all feasible solutions
			time_t now = time(0);
			// convert now to string form
			char* dt = ctime(&now);

			ofstream output_file;
			std::string pre_filename = "From_Node_";
			std::string Node_i_name = to_string(Node_i.Node_Index);
			std::string mid_filename = "_Expansion_At_";
			std::string Node_i_Child_name = dt;
			std::string post_filename = "_Opt_Soln.txt";

			std::string filename = pre_filename + Node_i_name + mid_filename + Node_i_Child_name + post_filename;
			output_file.open(filename, std::ofstream::out);
			for (int i = 0; i < Opt_Soln_Tot.size(); i++)
			{
				output_file<<Opt_Soln_Tot[i]<<endl;
			}
			output_file.close();
			int Start_Ind;
			std::vector<double> Kinetic_Energy_CMP;
			for (int i = 0; i < Feasible_Num; i++)
			{
				std::vector<double> Opt_Soln_Temp;
				Start_Ind = i * Variable_Num;
				for (int k = Start_Ind; k < Start_Ind + Variable_Num; k++)
				{
					Opt_Soln_Temp.push_back(Opt_Soln_Tot[k]);
				}
				dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj; double T_tot;
				Opt_Seed_Unzip(Opt_Soln_Temp, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);
				double KE_i = Kinetic_Energy_End_Frame(StateNDot_Traj);
				Kinetic_Energy_CMP.push_back(KE_i);
			}
			int Min_Ind = Minimum_Index(Kinetic_Energy_CMP);
			Start_Ind = Min_Ind *Variable_Num;
			int Trivial_Index = 0;
			for (int h = Start_Ind; h < Start_Ind + Variable_Num; h++)
			{
				Opt_Soln_Output[Trivial_Index] = Opt_Soln_Tot[h];
				Trivial_Index = Trivial_Index + 1;
			}
		}

	// for (int j = 0; j < Time_Queue.size(); j++) {
	// 	std::vector<double> Opt_Soln_Seed, Opt_Soln, ObjNConstraint_Val, ObjNConstraint_Type;
	// 	Time_Seed = Time_Queue[j];
	// 	cout<<" --------------------------------- Time Guess Iteration: "<<j<<" --------------------------------- "<<endl;
	// 	cout<<"The Guess for the time duaration is :" <<Time_Seed<<" s"<<endl;
	// 	Opt_Soln = Nodes_Optimization_Inner_Opt(Node_i, Node_i_child, Opt_Soln_Seed);
	// 	Nodes_Optimization_ObjNConstraint(Opt_Soln, ObjNConstraint_Val, ObjNConstraint_Type);
	// 	// For the plot of the unsuccessful case
	// 	// Opt_Soln_Write2Txt(Node_i, Node_i, Opt_Soln);
	// 	double ObjNConstraint_Violation_Val = ObjNConstraint_Violation(ObjNConstraint_Val, ObjNConstraint_Type);
	// 	cout<<endl;
	// 	cout<<" ----------------------------- ObjNConstraint_Violation_Val: "<<ObjNConstraint_Violation_Val<<" ------------------------------- "<<endl;
	// 	if(ObjNConstraint_Violation_Val<0.1)
	// 	{
	// 		Opt_Flag = 1;
	// 		time_t now = time(0);
	// 		// convert now to string form
	// 		char* dt = ctime(&now);
	// 		ofstream output_file;
	// 		std::string pre_filename = "Node_";
	// 		std::string Node_i_name = to_string(Node_i.Node_Index);
	// 		std::string mid_filename = "_Expansion_At_";
	// 		std::string post_filename = ".txt";
	//
	// 		std::string filename = pre_filename + Node_i_name + mid_filename + dt + post_filename;
	// 		output_file.open(filename, std::ofstream::out);
	// 		for (int i = 0; i < Opt_Soln.size(); i++)
	// 		{	output_file<<Opt_Soln[i]<<endl;}
	// 		output_file.close();
	// 	}
	// 	else
	// 	{	// The differentiation between the pure seed guess and the failure guess needs to be conducted
	// 		double Seed_offset = 0.0;
	// 		for (int i = 0; i < Opt_Soln.size(); i++) {
	// 			double state_minus = (Opt_Soln[i] - Opt_Soln_Seed[i]) * (Opt_Soln[i] - Opt_Soln_Seed[i]);
	// 			// cout<<"Opt_Soln[i]"<<Opt_Soln[i]<<endl;
	// 			Seed_offset = max(Seed_offset, state_minus);}
	// 		// if(Seed_offset>0.1)
	// 		// {	// Then this is a failure attempt
	// 		// 	time_t now = time(0);
	// 		// 	// convert now to string form
	// 		// 	char* dt = ctime(&now);
	// 		// 	ofstream output_file;
	// 		// 	std::string pre_filename = "Node_";
	// 		// 	std::string Node_i_name = to_string(Node_i.Node_Index);
	// 		// 	std::string mid_filename = "_Expansion_Failure_At_";
	// 		// 	std::string post_filename = ".txt";
	// 		//
	// 		// 	std::string filename = pre_filename + Node_i_name + mid_filename + dt + post_filename;
	// 		// 	output_file.open(filename, std::ofstream::out);
	// 		// 	for (int i = 0; i < Opt_Soln.size(); i++)
	// 		// 	{	output_file<<Opt_Soln[i]<<endl;}
	// 		// 	output_file.close();
	// 		// }
	// 	}
	// }
		return Opt_Flag;
	}

	std::vector<double> Final_State_Opt(std::vector<double> &Opt_Soln, Tree_Node &Node_i, Tree_Node &Node_i_child)
	{
		dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj; double T_tot;

		Opt_Seed_Unzip(Opt_Soln, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);

		dlib::matrix<double> Final_State = dlib::colm(StateNDot_Traj, StateNDot_Traj.nc()-1);

		std::vector<double> Robot_State_Final;			Robot_StateNDot StateNDot_Final_Opt;

		for (int i = 0; i < Final_State.nr(); i++)
		{
			Robot_State_Final.push_back(Final_State(i));
		}
		Robot_StateNDot Node_i_Ori = Node_i.Node_StateNDot;
		Structure_P.Node_i = Node_i;
		Structure_P.Node_i_child = Node_i_child;
		Node_i.Node_StateNDot = StateVec2StateNDot(Robot_State_Final);

		Robot_State_Final = Seed_Guess_Gene_Robotstate(Node_i, Node_i_child);

		for (int i = 0; i < Final_State.nr(); i++) {
			StateNDot_Traj(i, StateNDot_Traj.nc()-1) = Robot_State_Final[i];
		}
	// Now it is time to rewrite this vector back into the zip
		std::vector<double> Opt_Seed;
		Opt_Seed.push_back(T_tot);
		Opt_Seed_Zip(Opt_Seed, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);

		Node_i.Node_StateNDot = Node_i_Ori;
		return Opt_Seed;
	}
	void Opt_Soln_Write2Txt(Tree_Node &Node_i,Tree_Node &Node_i_child, std::vector<double> &Opt_Soln)
	{
		ofstream output_file;
		std::string pre_filename = "From_Node_";
		std::string Node_i_name = to_string(Node_i.Node_Index);
		std::string mid_filename = "_To_Node_";
		std::string Node_i_Child_name = to_string(Node_i_child.Node_Index);
		std::string post_filename = "_Opt_Soln.txt";
		std::string filename = pre_filename + Node_i_name + mid_filename + Node_i_Child_name + post_filename;
		output_file.open(filename, std::ofstream::out);
		for (int i = 0; i < Opt_Soln.size(); i++)
		{
			output_file<<Opt_Soln[i]<<endl;
		}
		output_file.close();
	}
	std::vector<double> Time_Seed_Queue_fn(double Time_Interval, int Total_Num)
	{
	// This function is used to generate the Time_Seed_Queue
		double Time_Center = 0.5;
		int One_Side_Num = (Total_Num - 1)/2;
		std::vector<double> Time_Seed_Queue;
		Time_Seed_Queue.push_back(Time_Center);
		for (int i = 1; i < One_Side_Num; i++) {
			Time_Seed_Queue.push_back(Time_Center + Time_Interval * i);
			Time_Seed_Queue.push_back(Time_Center - Time_Interval * i);
		}
	// std::vector<double> Time_Seed_Queue;
	// Time_Seed_Queue.push_back(Time_Center);
	// for (int i = 0; i < Total_Num; i++) {
	// 	Time_Seed_Queue.push_back(Time_Center + (i+1) * Time_Interval);
	// }
		return Time_Seed_Queue;
	}
	double ObjNConstraint_Violation(const std::vector<double> &ObjNConstraint_Val, const std::vector<double> &ObjNConstraint_Type)
	{
		double ObjNConstraint_Violation_Val = 0.0;
		for (int i = 0; i < ObjNConstraint_Val.size(); i++)
		{
			if(ObjNConstraint_Type[i]==0)
			{
				if(abs(ObjNConstraint_Val[i])>ObjNConstraint_Violation_Val)
				{
					ObjNConstraint_Violation_Val = abs(ObjNConstraint_Val[i]);
				}
			}
		}
		return ObjNConstraint_Violation_Val;
	}
	std::vector<double> Nodes_Optimization_Inner_Opt(Tree_Node &Node_i, Tree_Node &Node_i_child)
	{
		std::vector<double> Opt_Seed = Seed_Guess_Gene(Node_i, Node_i_child);
		std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
		Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	snoptProblem Nodes_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = Opt_Seed.size();
	integer neF = ObjNConstraint_Val.size();
	integer lenA  =  n * neF;

	Structure_P.Opt_Val_No = Opt_Seed.size();			Structure_P.ObjNConst_No = ObjNConstraint_Val.size();

	integer *iAfun = new integer[lenA];              	integer *jAvar = new integer[lenA];					doublereal *A  = new doublereal[lenA];
	integer lenG   = lenA;								integer *iGfun = new integer[lenG];					integer *jGvar = new integer[lenG];
	doublereal *x      = new doublereal[n];				doublereal *xlow   = new doublereal[n];				doublereal *xupp   = new doublereal[n];
	doublereal *xmul   = new doublereal[n];				integer    *xstate = new    integer[n];

	doublereal *F      = new doublereal[neF];			doublereal *Flow   = new doublereal[neF];			doublereal *Fupp   = new doublereal[neF];
	doublereal *Fmul   = new doublereal[neF];			integer    *Fstate = new integer[neF];

	integer nxnames = 1;								integer nFnames = 1;					char *xnames = new char[nxnames*8];					char *Fnames = new char[nFnames*8];

	integer    ObjRow = 0;								doublereal ObjAdd = 0;

	for (int i = 0; i < n; i++) {
		xlow[i] =-Inf;
		xupp[i] = Inf;
	}
	xlow[0] = 0.25;		xupp[0] = 3.5;
	int Index_Count = 1;
	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < 26; j++) {
			xlow[Index_Count] = xlow_vec(j);
			xupp[Index_Count] = xupp_vec(j);
			Index_Count = Index_Count + 1;}}

			for (int i = 0; i < Grids; i++) {
				for (int j = 0; j < 10; j++) {
					xlow[Index_Count] = ctrl_low_vec(j);
					xupp[Index_Count] = ctrl_upp_vec(j);
					Index_Count = Index_Count + 1;}}

					for (int i = 0; i < n; i++) {
						xstate[i] = 0.0;
		x[i] = Opt_Seed[i];  	// Initial guess
	}

	for(int i = 0; i<neF; i++)
	{
		if (ObjNConstraint_Type[i]>1.0)
		{
			Flow[i] = -Acc_max;
			Fupp[i] = Acc_max;
		}
		else
		{
			Flow[i] = 0.0;
			if(ObjNConstraint_Type[i]>0.0)
			{
				Fupp[i] = Inf;
			}
			else
			{
				Fupp[i] = 0.0;
			}
		}
	}

	// Here the linear initial condition matching is conducted as bounds
	std::vector<double> Init_Config = StateNDot2StateVec(Node_i.Node_StateNDot);

	for (int i = 0; i < 26; i++)
	{
		Flow[i+1] = Init_Config[i];
		Fupp[i+1] = Init_Config[i];
	}

	// Load the data for ToyProb ...
	Nodes_Optimization_Pr.setPrintFile  ( "Nodes_Optimization_Pr.out" );
	Nodes_Optimization_Pr.setProblemSize( n, neF );
	Nodes_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
	Nodes_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
	Nodes_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
	Nodes_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
	Nodes_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
	Nodes_Optimization_Pr.setXNames     ( xnames, nxnames );
	Nodes_Optimization_Pr.setFNames     ( Fnames, nFnames );
	Nodes_Optimization_Pr.setProbName   ( "Nodes_Optimization_Pr" );
	Nodes_Optimization_Pr.setUserFun    ( Nodes_Optimization_Pr_fn);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	Nodes_Optimization_Pr.computeJac    ();
	Nodes_Optimization_Pr.setIntParameter( "Derivative option", 0 );
	Nodes_Optimization_Pr.setIntParameter("Major iterations limit", 250);
	Nodes_Optimization_Pr.setIntParameter("Minor iterations limit", 50000);
	Nodes_Optimization_Pr.setIntParameter("Iterations limit", 2000000);
	Nodes_Optimization_Pr.setIntParameter("Minor print level", 1);
	// Nodes_Optimization_Pr.setSpecsFile   ( ".\Nodes_Optimization_Pr.spc" );

	// string QPsolver_opt = "QN";
	// Nodes_Optimization_Pr.setRealParameter("QPsolver", QPsolver_opt);


	integer Cold = 0, Basis = 1, Warm = 2;
	Nodes_Optimization_Pr.solve( Cold );

	std::vector<double> Opt_Soln;

	for (int i = 0; i < n; i++)
	{
		// cout<<"x[i]"<<x[i]<<endl;
		Opt_Soln.push_back(x[i]);
	}
	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	return Opt_Soln;
}

void Contact_Force_Bounds(std::vector<double> &sigma, std::vector<double> &Contact_Force_Status_i)
{
	Contact_Force_Status_i[0] = !sigma[0];				Contact_Force_Status_i[1] = !sigma[0];
	Contact_Force_Status_i[2] = !sigma[0];				Contact_Force_Status_i[3] = !sigma[0];

	Contact_Force_Status_i[4] = !sigma[1];				Contact_Force_Status_i[5] = !sigma[1];
	Contact_Force_Status_i[6] = !sigma[1];				Contact_Force_Status_i[7] = !sigma[1];

	Contact_Force_Status_i[8] = !sigma[2];				Contact_Force_Status_i[9] = !sigma[2];
	Contact_Force_Status_i[10] = !sigma[3];				Contact_Force_Status_i[11] = !sigma[3];
}
int Nodes_Optimization_Pr_fn(integer    *Status, integer *n,    doublereal x[],
	integer    *needF,  integer *neF,  doublereal F[],
	integer    *needG,  integer *neG,  doublereal G[],
	char       *cu,     integer *lencu,
	integer    iu[],    integer *leniu,
	doublereal ru[],    integer *lenru )
{	 std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
	const int Opt_Val_No = Structure_P.Opt_Val_No ;				const int ObjNConst_No = Structure_P.ObjNConst_No;
	for (int i = 0; i < Opt_Val_No; i++)
	{
		Opt_Seed.push_back(x[i]);
	}

	Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	for (int i = 0; i < ObjNConst_No; i++)
	{
		F[i] = ObjNConstraint_Val[i];
	}

	return 0;
}

std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
	// This function will generate the spline coefficients needed for the further optimization
	// Here Time_Seed is a global variable whose value is changed in the outer loop

	double T = Time_Seed;

	// The first step is to generate a feasible configuration that can satisfy the contact mode in the node i child
	std::vector<double> Init_Config = StateNDot2StateVec(Node_i.Node_StateNDot);
	std::vector<double> Seed_Config = Seed_Guess_Gene_Robotstate(Node_i, Node_i_child);
	const int StateNDot_len = Init_Config.size();			const int Control_len = 10;			const int Contact_Force_len = 12;

	dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;
	StateNDot_Traj = dlib::zeros_matrix<double>(StateNDot_len, Grids);				// each variable traj will be in a row fashion
	Ctrl_Traj = dlib::zeros_matrix<double>(Control_len, Grids);
	Contact_Force_Traj = dlib::zeros_matrix<double>(Contact_Force_len, Grids);

	dlib::matrix<double> StateNDot_Coeff, Ctrl_Coeff, Contact_Force_Coeff;
	StateNDot_Coeff = dlib::zeros_matrix<double>(4*StateNDot_len, Grids-1);
	Ctrl_Coeff = dlib::zeros_matrix<double>(2*Control_len, Grids-1);
	Contact_Force_Coeff = dlib::zeros_matrix<double>(2*Contact_Force_len, Grids-1);

	// The only step that needs to be changed with the new idea is the following code
	// Now we have the Init_Config and the Seed_Config so how to interpolate the intermediate trajectories
	//
	// double T_tot = T * (Grids - 1);
	// std::vector<double> State_Traj_Coeff_i;
	// double Pos_i_Init, Pos_i_Goal, Vel_i_Init, Vel_i_Goal;
	// for (int i = 0; i < 13; i++)
	// {
	// 	Pos_i_Init = Init_Config[i];				Pos_i_Goal = Seed_Config[i];
	// 	Vel_i_Init = Init_Config[i+13];			Vel_i_Goal = Seed_Config[i+13];
	// 	State_Traj_Coeff_i = CubicSpline_Coeff_fn(T_tot, Pos_i_Init, Pos_i_Goal, Vel_i_Init, Vel_i_Goal); // 4 by 1 vector: a, b, c, d
	// 	// After this step, we hacve the cubic spline for the whole trajectories. Then it is time to discretize them
	// 	double ds = 1.0/(Grids * 1.0 - 1.0), s_j;
	// 	for (int j = 0; j < Grids; j++)
	// 	{
	// 		s_j = ds * (j);
	// 		// cout<<s_j<<endl;
	// 		StateNDot_Traj(i,j) = CubicSpline_Evaluation_fn(State_Traj_Coeff_i, s_j);
	// 		StateNDot_Traj( i+ 13,j) = CubicSpline_1stOrder_Evaluation_fn(State_Traj_Coeff_i, s_j, T_tot);
	// 	}
	// }
	// cout<<StateNDot_Traj<<endl;

	dlib::matrix<double> Robot_State_Interpol_i;
	for (int i = 0; i < StateNDot_len; i++){
		Robot_State_Interpol_i = dlib::linspace(Init_Config[i], Seed_Config[i], Grids);
		for (int j = 0; j < Grids; j++){
			StateNDot_Traj(i,j) = Robot_State_Interpol_i(j);}}

	// Here is the initialization of the spline coefficients for the stateNdot, control and contact force
	// First is to initialize: StateNDot
			double x_init, x_end, xdot_init, xdot_end;
			std::vector<double> CubicSpline_Coeff_Pos_vec, CubicSpline_Coeff_Vel_vec, PVA_init, PVA_end;
			for (int i = 0; i < Grids-1; i++)
			{
				for (int j = 0; j < StateNDot_len/2; j++)
				{
		  // Position first
					x_init = StateNDot_Traj(j,i);
					x_end = StateNDot_Traj(j,i+1);
					xdot_init = StateNDot_Traj(j+StateNDot_len/2,i);
					xdot_end = StateNDot_Traj(j+StateNDot_len/2,i+1);
			CubicSpline_Coeff_Pos_vec = CubicSpline_Coeff_fn(T, x_init, x_end, xdot_init, xdot_end); // 4 by 1 vector: a, b, c, d

			// Velocity second
			PVA_init = CubicSpline_PosVelAcc4(T, CubicSpline_Coeff_Pos_vec[0], CubicSpline_Coeff_Pos_vec[1], CubicSpline_Coeff_Pos_vec[2], CubicSpline_Coeff_Pos_vec[3], 0);
			PVA_end  = CubicSpline_PosVelAcc4(T, CubicSpline_Coeff_Pos_vec[0], CubicSpline_Coeff_Pos_vec[1], CubicSpline_Coeff_Pos_vec[2], CubicSpline_Coeff_Pos_vec[3], 1);
			CubicSpline_Coeff_Vel_vec = CubicSpline_Coeff_fn(T, xdot_init,xdot_end, PVA_init[2], PVA_end[2]); // 4 by 1 vector: a, b, c, d

			StateNDot_Coeff(8*j,i)   = CubicSpline_Coeff_Pos_vec[0];
			StateNDot_Coeff(8*j+1,i) = CubicSpline_Coeff_Pos_vec[1];
			StateNDot_Coeff(8*j+2,i) = CubicSpline_Coeff_Pos_vec[2];
			StateNDot_Coeff(8*j+3,i) = CubicSpline_Coeff_Pos_vec[3];
			StateNDot_Coeff(8*j+4,i) = CubicSpline_Coeff_Vel_vec[0];
			StateNDot_Coeff(8*j+5,i) = CubicSpline_Coeff_Vel_vec[1];
			StateNDot_Coeff(8*j+6,i) = CubicSpline_Coeff_Vel_vec[2];
			StateNDot_Coeff(8*j+7,i) = CubicSpline_Coeff_Vel_vec[3];
		}
	}
	// cout<<StateNDot_Coeff<<endl;

	// Second is to initialize: Control and Contact Force
	std::vector<double> Robot_Pos(13), Robot_Vel(13), Robot_VelfromPos(13), Robotstate_Vec_i;		dlib::matrix<double> Robot_Acc; Robot_Acc = dlib::zeros_matrix<double>(13,1);
	dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full, Jac_Full_Trans, Dynamics_LHS, Dynamics_RHS, Dynamics_RHS_Matrix;		Robot_StateNDot Robot_StateNDot_i;
	for (int i = 0; i < Grids-1; i++){
		Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(T, StateNDot_Coeff, i, 0, Robot_Pos, Robot_Vel, Robot_Acc, Robot_VelfromPos);
		Robotstate_Vec_i = PosNVel2StateVec(Robot_Pos, Robot_Vel);
		Robot_StateNDot_i = StateVec2StateNDot(Robotstate_Vec_i);
		Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
		Dynamics_LHS = D_q * Robot_Acc + C_q_qdot;
		Dynamics_RHS_Matrix = Dynamics_RHS_Matrix_fn(Jac_Full, B_q);
		Dynamics_RHS = dlib::pinv(Dynamics_RHS_Matrix) * Dynamics_LHS;
		for (int j = 0; j < Dynamics_RHS.nr(); j++) {
			if(j<Contact_Force_Traj.nr())
				{	Contact_Force_Traj(j,i) = Dynamics_RHS(j);}
			else
			{
				Ctrl_Traj(j - Contact_Force_Traj.nr(),i) = Dynamics_RHS(j);
			}
		}
	}
	// cout<<StateNDot_Traj<<endl;			cout<<Contact_Force_Traj<<endl;			cout<<Ctrl_Traj<<endl;			cout<<StateNDot_Coeff<<endl;
	// THen the last column
	Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(T, StateNDot_Coeff, Grids-2, 1, Robot_Pos, Robot_Vel, Robot_Acc, Robot_VelfromPos);
	Robotstate_Vec_i = PosNVel2StateVec(Robot_Pos, Robot_Vel);
	Robot_StateNDot_i = StateVec2StateNDot(Robotstate_Vec_i);
	Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
	Dynamics_LHS = D_q * Robot_Acc + C_q_qdot;
	Dynamics_RHS_Matrix = Dynamics_RHS_Matrix_fn(Jac_Full, B_q);
	Dynamics_RHS = dlib::pinv(Dynamics_RHS_Matrix) * Dynamics_LHS;
	for (int j = 0; j < Dynamics_RHS.nr(); j++) {
		if(j<Contact_Force_Traj.nr())
			{	Contact_Force_Traj(j,Grids-1) = Dynamics_RHS(j);}
		else
		{
			Ctrl_Traj(j - Contact_Force_Traj.nr(),Grids-1) = Dynamics_RHS(j);
		}
	}
	// cout<<StateNDot_Traj<<endl;			cout<<Ctrl_Traj<<endl;				cout<<Contact_Force_Traj<<endl;
	// The calculation of the coefficients of the control and contact force is easier compared to the robot state due to the assumption of the linear equation
	// Ctrl_Contact_Force_Coeff_fn(Ctrl_Traj, Contact_Force_Traj, Ctrl_Coeff, Contact_Force_Coeff);

	// The final task is to pile them into a single vector
	std::vector<double> Opt_Seed;
	Opt_Seed.push_back(T * (Grids - 1));
	Opt_Seed_Zip(Opt_Seed, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);
	// cout<<StateNDot_Coeff<<endl;				cout<<Ctrl_Coeff<<endl;			cout<<Contact_Force_Coeff<<endl;
	return Opt_Seed;
}

void Opt_Seed_Unzip(std::vector<double> &Opt_Seed, double &T_tot, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj)
{
	const int NumOfStateNDot_Traj = 26;
	const int NumOfCtrl_Traj = 10;
	const int NumOfContactForce_Traj = 12;

	T_tot = Opt_Seed[0];			int Opt_Seed_Index = 1;

	StateNDot_Traj = dlib::zeros_matrix<double>(NumOfStateNDot_Traj, Grids);
	Ctrl_Traj = dlib::zeros_matrix<double>(NumOfCtrl_Traj, Grids);
	Contact_Force_Traj = dlib::zeros_matrix<double>(NumOfContactForce_Traj, Grids);

	// 1. Retrieve the StateNDot_Traj matrix
	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < NumOfStateNDot_Traj; j++) {
			StateNDot_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<StateNDot_Traj<<endl;
	// 2. Retrieve the control matrix
			for (int i = 0; i < Grids; i++) {
				for (int j = 0; j < NumOfCtrl_Traj; j++) {
					Ctrl_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
					Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<Ctrl_Traj<<endl;
	// 3. Retrieve the contact force matrix
					for (int i = 0; i < Grids; i++) {
						for (int j = 0; j < NumOfContactForce_Traj; j++) {
							Contact_Force_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
							Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<Contact_Force_Traj<<endl;
							return ;
						}
						void Opt_Seed_Zip(std::vector<double> &Opt_Seed, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj)
{	// This function is used to stack the coefficient matrices into a column vector
	for (int i = 0; i < StateNDot_Traj.nc(); i++) {
		for (int j = 0; j < StateNDot_Traj.nr(); j++) {
			Opt_Seed.push_back(StateNDot_Traj(j,i));
		}
	}
	for (int i = 0; i < Ctrl_Traj.nc(); i++) {
		for (int j = 0; j < Ctrl_Traj.nr(); j++) {
			Opt_Seed.push_back(Ctrl_Traj(j,i));
		}
	}
	for (int i = 0; i < Contact_Force_Traj.nc(); i++) {
		for (int j = 0; j < Contact_Force_Traj.nr(); j++) {
			Opt_Seed.push_back(Contact_Force_Traj(j,i));
		}
	}
}
void Ctrl_Contact_Force_Coeff_fn(dlib::matrix<double> &Ctrl_Traj, dlib::matrix<double> &Contact_Force_Traj, dlib::matrix<double> &Ctrl_Coeff, dlib::matrix<double> &Contact_Force_Coeff)
{
	// y = a * s + b
	double Ctrl_init, Ctrl_end, Contact_Force_init, Contact_Force_end;
	for (int i = 0; i < Grids-1; i++){
		// Computation of the control coeff

		for (int j = 0; j < 10; j++) {
			Ctrl_init = Ctrl_Traj(j,i);
			Ctrl_end = Ctrl_Traj(j,i+1);
			Ctrl_Coeff(2*j,i) = Ctrl_end - Ctrl_init;
			Ctrl_Coeff(2*j+1,i) = Ctrl_init;
		}
		// Computation of the contact force coeff
		for (int j = 0; j < 12; j++) {
			Contact_Force_init = Contact_Force_Traj(j,i);
			Contact_Force_end = Contact_Force_Traj(j,i+1);
			Contact_Force_Coeff(2*j,i) = Contact_Force_end - Contact_Force_init;
			Contact_Force_Coeff(2*j+1,i) = Contact_Force_init;
		}
	}
	return;
}
dlib::matrix<double> Dynamics_RHS_Matrix_fn(dlib::matrix<double> &Jac_Full, dlib::matrix<double> &B_q)
{
	dlib::matrix<double> Jac_Full_Trans, Dynamics_RHS_Matrix;
	Jac_Full_Trans = dlib::trans(Jac_Full);

	const int Dynamics_RHS_Matrix_Row = Jac_Full_Trans.nr();
	const int Dynamics_RHS_Matrix_Col = Jac_Full_Trans.nc() + B_q.nc();
	Dynamics_RHS_Matrix = dlib::zeros_matrix<double>(Dynamics_RHS_Matrix_Row, Dynamics_RHS_Matrix_Col);
	for (int i = 0; i < Dynamics_RHS_Matrix_Col; i++) {
		if (i<Jac_Full_Trans.nc())
		{

			dlib::set_colm(Dynamics_RHS_Matrix, i) = dlib::colm(Jac_Full_Trans,i);
		}
		else{
			dlib::set_colm(Dynamics_RHS_Matrix, i) = dlib::colm(B_q,i-Jac_Full_Trans.nc());
		}
	}
	return Dynamics_RHS_Matrix;
}
std::vector<double> PosNVel2StateVec(std::vector<double> & Pos, std::vector<double> & Vel)
{
	std::vector<double> StateVec_i;
	for (int i = 0; i < Pos.size(); i++) {
		StateVec_i.push_back(Pos[i]);
	}
	for (int i = 0; i < Vel.size(); i++) {
		StateVec_i.push_back(Vel[i]);
	}
	return StateVec_i;
}
void Dynamics_Matrices(const Robot_StateNDot &Node_StateNDot, dlib::matrix<double> &D_q, dlib::matrix<double> &B_q, dlib::matrix<double> &C_q_qdot, dlib::matrix<double> &Jac_Full)
{
	D_q = D_q_fn(Node_StateNDot);
	B_q = B_q_fn();
	C_q_qdot = C_q_qdot_fn(Node_StateNDot);
	Jac_Full = Jac_Full_fn(Node_StateNDot);
}
void Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(double T, dlib::matrix<double> &StateNDot_Coeff, int Grid_Ind, double s, std::vector<double> &Robot_Config,  std::vector<double> &Robot_Vel, dlib::matrix<double> &Robot_Acc, std::vector<double> &Robot_VelfromPos)
{
	//Here Grid_Ind denotes which Grid we are talking about and s is a value between 0 and 1
	Robot_Acc = dlib::zeros_matrix<double>(13,1);
	std::vector<double> PVAVP_i;
	double x_a, x_b, x_c, x_d, xdot_a, xdot_b, xdot_c, xdot_d;
	for (int i = 0; i < 13; i++) {
		x_a = StateNDot_Coeff(8*i, Grid_Ind);
		x_b = StateNDot_Coeff(8*i+1, Grid_Ind);
		x_c = StateNDot_Coeff(8*i+2, Grid_Ind);
		x_d = StateNDot_Coeff(8*i+3, Grid_Ind);

		xdot_a = StateNDot_Coeff(8*i+4, Grid_Ind);
		xdot_b = StateNDot_Coeff(8*i+5, Grid_Ind);
		xdot_c = StateNDot_Coeff(8*i+6, Grid_Ind);
		xdot_d = StateNDot_Coeff(8*i+7, Grid_Ind);

		PVAVP_i = CubicSpline_PosVelAcc8(T, x_a, x_b, x_c, x_d, xdot_a, xdot_b, xdot_c, xdot_d, s);
		Robot_Config[i] = PVAVP_i[0];
		Robot_Vel[i] = PVAVP_i[1];
		Robot_Acc(i) = PVAVP_i[2];
		Robot_VelfromPos[i] = PVAVP_i[3];
	}
	return;
}
double CubicSpline_Evaluation_fn(const std::vector<double> &CubicSpline_Coeff, double s)
{
	double a, b, c, d, y;
	a = CubicSpline_Coeff[0];			b = CubicSpline_Coeff[1];
	c = CubicSpline_Coeff[2];			d = CubicSpline_Coeff[3];
	y = a * s * s * s + b * s * s + c * s + d;
	return y;
}
double CubicSpline_1stOrder_Evaluation_fn(const std::vector<double> &CubicSpline_Coeff, double s, double T)
{ 	// This function is used to output the first order derivative value of a given spline at s position
	double a, b, c, d, y;
	a = CubicSpline_Coeff[0];			b = CubicSpline_Coeff[1];
	c = CubicSpline_Coeff[2];			d = CubicSpline_Coeff[3];
	y = (3 * a * s * s + 2 * b * s + c)/T;
	return y;
}

std::vector<double> CubicSpline_PosVelAcc4(double T, double a, double b, double c, double d, double s)
{
	//# This function is used to calcualte the position, velocity and acceleration given the spline coefficients
    //# Here T is the duration constant, s is the path variable
	std::vector<double> PVA(3);
	double Pos = a*s*s*s + b*s*s + c*s + d;
	double Vel = (3*a*s*s + 2*b*s + c)/T;
	double Acc = (6*a*s+2*b)/(T*T);
	PVA[0] = Pos; 			PVA[1] = Vel; 			PVA[2] = Acc;
	return PVA;
}
std::vector<double> CubicSpline_PosVelAcc8(double T, double x_a, double x_b, double x_c, double x_d, double xdot_a, double xdot_b, double xdot_c, double xdot_d, double s)
{
	std::vector<double> PVAVP(4);
	double Pos = x_a*s*s*s + x_b*s*s + x_c*s + x_d;
	double Vel =  xdot_a*s*s*s + xdot_b*s*s + xdot_c*s + xdot_d;
	double Acc = (3*xdot_a*s*s + 2*xdot_b*s + xdot_c)/T;
	double VelfromPos = (3*x_a*s*s + 2*x_b*s + x_c)/T;
	PVAVP[0] = Pos;				PVAVP[1] = Vel;				PVAVP[2] = Acc;				PVAVP[3] = VelfromPos;
	return PVAVP;
}
std::vector<double> CubicSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init, double xdot_end)
{	// This function is used to calcualte the coefficients for the cubic spline
	// The cubic spline is expressed to be : y(s) = a*s^3 + b*s^2 + c*s + d
	std::vector<double> CubicSpline_Coeff_vec(4);
	double a, b, c, d;
	a = 2*x_init - 2*x_end + T*xdot_end + T*xdot_init;
	b = 3*x_end - 3*x_init - T*xdot_end - 2*T*xdot_init;
	c = T*xdot_init;
	d = x_init;
	CubicSpline_Coeff_vec[0] = a;
	CubicSpline_Coeff_vec[1] = b;
	CubicSpline_Coeff_vec[2] = c;
	CubicSpline_Coeff_vec[3] = d;
	// cout<<a<<endl;
	// cout<<b<<endl;
	// cout<<c<<endl;
	// cout<<d<<endl;

	return CubicSpline_Coeff_vec;
}
std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
	// This function is used to generate a configuration to initialize the optimization
	std::vector<double> Robot_State_Seed = StateNDot2StateVec(Node_i.Node_StateNDot);
	std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
	Seed_Conf_Optimization_ObjNConstraint(Robot_State_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	snoptProblem Seed_Conf_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = Robot_State_Seed.size();
	integer neF = ObjNConstraint_Val.size();
	integer lenA  =  n * neF;

	integer *iAfun = new integer[lenA];              integer *jAvar = new integer[lenA];					doublereal *A  = new doublereal[lenA];

	integer lenG   = lenA;							integer *iGfun = new integer[lenG];						integer *jGvar = new integer[lenG];

	doublereal *x      = new doublereal[n];			doublereal *xlow   = new doublereal[n];					doublereal *xupp   = new doublereal[n];
	doublereal *xmul   = new doublereal[n];			integer    *xstate = new    integer[n];

	doublereal *F      = new doublereal[neF];		doublereal *Flow   = new doublereal[neF];				doublereal *Fupp   = new doublereal[neF];
	doublereal *Fmul   = new doublereal[neF];		integer    *Fstate = new integer[neF];

	integer nxnames = 1;							integer nFnames = 1;				char *xnames = new char[nxnames*8];					char *Fnames = new char[nFnames*8];

	integer    ObjRow = 0;							doublereal ObjAdd = 0;

	for (int i = 0; i < n; i++) {
		xlow[i] = xlow_vec(i);						xupp[i] = xupp_vec(i);				xstate[i] = 0.0;							x[i] = Robot_State_Seed[i];  	// Initial guess
	}

	for(int i = 0; i<neF; i++){
		// The lower bound is the same
		Flow[i] = 0.0;
		if(ObjNConstraint_Type[i]>0)	// Inequality constraint
			{	Fupp[i] = Inf;}
		else{
			Fupp[i] = 0.0;}
		}

	// Load the data for ToyProb ...
		Seed_Conf_Optimization_Pr.setPrintFile  ( "Seed_Conf_Optimization_Pr.out" );
		Seed_Conf_Optimization_Pr.setProblemSize( n, neF );
		Seed_Conf_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
		Seed_Conf_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
		Seed_Conf_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
		Seed_Conf_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
		Seed_Conf_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
		Seed_Conf_Optimization_Pr.setXNames     ( xnames, nxnames );
		Seed_Conf_Optimization_Pr.setFNames     ( Fnames, nFnames );
		Seed_Conf_Optimization_Pr.setProbName   ( "Seed_Conf_Optimization_Pr" );
		Seed_Conf_Optimization_Pr.setUserFun    ( Seed_Conf_Optimization_Pr_fn_);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
		Seed_Conf_Optimization_Pr.computeJac    ();
		Seed_Conf_Optimization_Pr.setIntParameter( "Derivative option", 0 );
		Seed_Conf_Optimization_Pr.setIntParameter( "Print level", 0 );
		Seed_Conf_Optimization_Pr.setIntParameter( "Major print level", 0 );
		Seed_Conf_Optimization_Pr.setIntParameter( "Minor print level", 0 );

		integer Cold = 0, Basis = 1, Warm = 2;
		Seed_Conf_Optimization_Pr.solve( Cold );
		for (int i = 0; i < 26; i++){Robot_State_Seed[i] = x[i];}
			Robot_StateNDot Robot_StateNDot_Seed(Robot_State_Seed);
		std::string input_name = "Seed Configuration in Opt";
		Robot_Plot_fn(Robot_StateNDot_Seed, input_name);
		delete []iAfun;  delete []jAvar;  delete []A;		delete []iGfun;  delete []jGvar;
		delete []x;      delete []xlow;   delete []xupp;	delete []xmul;   delete []xstate;
		delete []F;      delete []Flow;   delete []Fupp;	delete []Fmul;   delete []Fstate;
		delete []xnames; delete []Fnames;
		return Robot_State_Seed;
	}
	int Seed_Conf_Optimization_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
		integer    *needF,  integer *neF,  doublereal F[],
		integer    *needG,  integer *neG,  doublereal G[],
		char       *cu,     integer *lencu,
		integer    iu[],    integer *leniu,
		doublereal ru[],    integer *lenru )
	{	 std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
		for (int i = 0; i < 26; i++){
			Opt_Seed.push_back(x[i]);}
			Seed_Conf_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
			for (int i = 0; i < ObjNConstraint_Val.size(); i++){
				F[i] = ObjNConstraint_Val[i];}
				return 0;
			}
			void Seed_Conf_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
			{
				Robot_StateNDot StateNDot_Init_i(Opt_Seed);		dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel, End_Effector_Pos_ref, End_Effector_Vel_ref;
				End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);

				End_Effector_Pos_ref = Structure_P.Node_i.End_Effector_Pos;
				End_Effector_Vel_ref = Structure_P.Node_i.End_Effector_Vel;

				Robot_StateNDot StateNDot_Init_ref = Structure_P.Node_i.Node_StateNDot;
				std::vector<double> StateVec_Init_ref = StateNDot2StateVec(StateNDot_Init_ref);

				std::vector<double> sigma_i = Structure_P.Node_i.sigma;
				std::vector<double> sigma_i_child = Structure_P.Node_i_child.sigma;

				std::vector<double> rCOM_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rCOM");
				std::vector<double> rA_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rA");
				std::vector<double> rB_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rB");
				std::vector<double> rC_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rC");
				std::vector<double> rD_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rD");
				std::vector<double> rE_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rE");
				std::vector<double> rF_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rF");
				std::vector<double> rT_ref = Ang_Pos_fn(Structure_P.Node_i.Node_StateNDot, "rT");

	// cout<<rA_ref[0]<<endl;
	// cout<<rB_ref[0]<<endl;
	// cout<<rC_ref[0]<<endl;
	// cout<<rD_ref[0]<<endl;
	// cout<<rE_ref[0]<<endl;
	// cout<<rF_ref[0]<<endl;

				std::vector<double> rCOM_opt = Ang_Pos_fn(StateNDot_Init_i, "rCOM");
				std::vector<double> rA_opt = Ang_Pos_fn(StateNDot_Init_i, "rA");
				std::vector<double> rB_opt = Ang_Pos_fn(StateNDot_Init_i, "rB");
				std::vector<double> rC_opt = Ang_Pos_fn(StateNDot_Init_i, "rC");
				std::vector<double> rD_opt = Ang_Pos_fn(StateNDot_Init_i, "rD");
				std::vector<double> rE_opt = Ang_Pos_fn(StateNDot_Init_i, "rE");
				std::vector<double> rF_opt = Ang_Pos_fn(StateNDot_Init_i, "rF");
				std::vector<double> rT_opt = Ang_Pos_fn(StateNDot_Init_i, "rT");

				std::vector<double> vCOM_ref = Ang_Vel_fn(StateNDot_Init_i,"vCOM");


				ObjNConstraint_Val.push_back(0);
				ObjNConstraint_Type.push_back(1);

				int Opt_Type_Flag, Crit_Grid;

				Sigma_TransNGoal(sigma_i, sigma_i_child, Opt_Type_Flag, Crit_Grid);

				double sigma_diff = 0.0; double Obj_val = 0.0;
				for (int i = 0; i < 4; i++) {
					sigma_diff = sigma_diff + (sigma_i_child[i] - sigma_i[i]);
				}
				if(abs(sigma_diff)>0)
				{
					Obj_val = 0.0;

		// In this case, there must a contact modification during this whole process
					for (int i = 0; i < StateVec_Init_ref.size(); i++)
					{
						Obj_val = Obj_val + (StateVec_Init_ref[i] - Opt_Seed[i]) * (StateVec_Init_ref[i] - Opt_Seed[i]);
					}

					if(sigma_diff<0)
					{
			// Contact break/retract this part is a little tough to handle. However, in our problem, we only consider the foot contact retract
						int Contact_Index = Sigma_Change(sigma_i, sigma_i_child);
			// double vCOM_sign = vCOM_ref[0]/abs(vCOM_ref[0]);
						if(Contact_Index<2)
						{
				// In this case, the robot will lift up left/right foot
							if(Contact_Index == 0)
							{
					// In this case, the robot will lift AB foot so we would like to move the COM to CD
								ObjNConstraint_Val.push_back(rCOM_opt[0] - rD_ref[0] - mini);
								ObjNConstraint_Type.push_back(1);
								ObjNConstraint_Val.push_back(rC_ref[0] - rCOM_opt[0] - mini);
								ObjNConstraint_Type.push_back(1);
							}
							else
							{
					// In this case, the robot will lift CD foot so we would like to move the COM to AB
								ObjNConstraint_Val.push_back(rCOM_opt[0] - rB_ref[0] - mini);
								ObjNConstraint_Type.push_back(1);
								ObjNConstraint_Val.push_back(rA_ref[0] - rCOM_opt[0] - mini);
								ObjNConstraint_Type.push_back(1);
							}
						}
					}
					else
					{
			// ObjNConstraint_Val.push_back((rCOM_opt[0] - rCOM_ref[0]) * (rCOM_opt[0] - rCOM_ref[0]));
			// ObjNConstraint_Type.push_back(0);

			// Making contact
						int Contact_Index = Sigma_Change(sigma_i, sigma_i_child);
			// cout<<Contact_Index<<endl;
			// The next step is to figure out which foot to step and which direction to go
						double vCOM_sign = vCOM_ref[0]/abs(vCOM_ref[0]);
						if(Contact_Index<2)
						{
							if(Contact_Index > 0)
							{
								if(vCOM_sign>0)
								{
						// ObjNConstraint_Val.push_back(rD_opt[0] - rA_opt[0] - mini);
						// ObjNConstraint_Type.push_back(1);
									ObjNConstraint_Val.push_back(rD_opt[0] - rD_ref[0] - mini);
									ObjNConstraint_Type.push_back(0);

								}
								else
								{
									ObjNConstraint_Val.push_back(rB_opt[0] - rC_opt[0] - mini);
									ObjNConstraint_Type.push_back(1);
								}
							}
							else
							{
								if(vCOM_sign>0)
								{
									ObjNConstraint_Val.push_back(rB_opt[0] - rC_opt[0] - mini);
						// cout<<rB_opt[0] - rC_opt[0]<<endl;
									ObjNConstraint_Type.push_back(1);
								}
								else
								{
									ObjNConstraint_Val.push_back(rD_opt[0] - rA_opt[0] - mini);
									ObjNConstraint_Type.push_back(1);
								}

								ObjNConstraint_Val.push_back((rT_opt[0] - rT_ref[0]) * (rT_opt[0] - rT_ref[0]));
								ObjNConstraint_Type.push_back(0);
							}
						}
					}

		// Obj_val = Kinetic_Energy_fn(StateNDot_Init_i);
		// Obj_val = (rCOM_ref[0] - rCOM_opt[0]) * (rCOM_ref[0] - rCOM_opt[0]) + (rCOM_ref[1] - rCOM_opt[1]) * (rCOM_ref[1] - rCOM_opt[1]);

		// if(Opt_Type_Flag==-1)
		// {
		// 	Obj_val = (Opt_Seed[2] - PI/2.0) * (Opt_Seed[2] - PI/2.0) ;
		// }
				}
				else
				{
		// This is used for the self-stasbilization process so kinetic energy is the only cost
		// for (int i = 0; i < StateVec_Init_ref.size(); i++)
		// {
		// 	Obj_val = Obj_val + (StateVec_Init_ref[i] - Opt_Seed[i]) * (StateVec_Init_ref[i] - Opt_Seed[i]);
		// }

					Obj_val = Kinetic_Energy_fn(StateNDot_Init_i);
		// double KE_End = Kinetic_Energy_fn(StateNDot_Init_i);
		// ObjNConstraint_Val.push_back(0.001 -  KE_End);
		// ObjNConstraint_Type.push_back(1);

		// ObjNConstraint_Val.push_back(rB_opt[0] - rC_opt[0] - mini);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(rB_opt[0] - rC_opt[0] - mini);
		// ObjNConstraint_Type.push_back(1);
		// ObjNConstraint_Val.push_back(Opt_Seed[2] * Opt_Seed[2]);
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.rIx - 0.2351)*(StateNDot_Init_i.rIx - 0.2351)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.rIy - 0.5363)*(StateNDot_Init_i.rIy - 0.5363)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.theta - 0.2521) * (StateNDot_Init_i.theta - 0.2521)) ;
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q1 + 1.8634)*(StateNDot_Init_i.q1 + 1.8634)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q2 - 2.4521)*(StateNDot_Init_i.q2 - 2.4521)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q3 + 1.1836)*(StateNDot_Init_i.q3 + 1.1836)) ;
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q4 + 1.0312)*(StateNDot_Init_i.q4 + 1.0312)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q5 - 1.6671)*(StateNDot_Init_i.q5 - 1.6671)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q6 + 0.8849)*(StateNDot_Init_i.q6 + 0.8849)) ;
		// ObjNConstraint_Type.push_back(0);
		//
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q7 - 1.0470)*(StateNDot_Init_i.q7 - 1.0470)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q8 + 1.0)*(StateNDot_Init_i.q8 + 1.0)) ;
		// ObjNConstraint_Type.push_back(0);
		// ObjNConstraint_Val.push_back((StateNDot_Init_i.q9 + 0.6965)*(StateNDot_Init_i.q9 + 0.6965)) ;
		// ObjNConstraint_Type.push_back(0);

				}
	// ObjNConstraint_Val.push_back(-StateNDot_Init_i.q1 + StateNDot_Init_i.q4 - mini);
	// ObjNConstraint_Type.push_back(1);
	//
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q2 -PI/4)*(StateNDot_Init_i.q2 -PI/6)) ;
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q4 + PI/2)*(StateNDot_Init_i.q4 + PI/2)) ;
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.theta - PI/2.0)*(StateNDot_Init_i.theta - PI/2.0)) ;
	// ObjNConstraint_Type.push_back(0);

				double Hand_max = max(rE_ref[0], rF_ref[0]);
				Hand_max = max(Hand_max, rD_ref[0]);
				Hand_max = max(Hand_max, rC_ref[0]);
				Hand_max = max(Hand_max, rA_ref[0]);
				Hand_max = max(Hand_max, rT_opt[0]);

				ObjNConstraint_Val[0] = Obj_val;
	//
	// ObjNConstraint_Val.push_back((Opt_Seed[2]) * (Opt_Seed[2]));
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q7- PI/2.0) * (StateNDot_Init_i.q7 - PI/2.0));
	// ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q8- 0*PI/2.0) * (StateNDot_Init_i.q8 - 0*PI/2.0));
	// ObjNConstraint_Type.push_back(0);
	// //
	// ObjNConstraint_Val.push_back((Opt_Seed[2] - PI/2) * (Opt_Seed[2] - PI/2));
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q5);
	// ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back(StateNDot_Init_i.q2);
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back((rF_ref[0] - rF_opt[0])*(rF_ref[0] - rF_opt[0]));
	// ObjNConstraint_Type.push_back(0);
	// ObjNConstraint_Val.push_back((StateNDot_Init_i.q8- 0*PI/2.0) * (StateNDot_Init_i.q8 - 0*PI/2.0));
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back((!sigma_i[2]) * (sigma_i_child[2]) * (rE_opt[0] - Hand_max) * (rE_opt[0] - Hand_max));
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back((!sigma_i[3]) * (sigma_i_child[3]) * (rF_opt[0] - Hand_max) * (rF_opt[0] - Hand_max));
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back((!sigma_i[2]) * (sigma_i_child[2]) * (rE_opt[1] - rT_opt[1]) * (rE_opt[1] - rT_opt[1]));
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back((!sigma_i[3]) * (sigma_i_child[3]) * (rF_opt[1] - rT_opt[1]) * (rF_opt[1] - rT_opt[1]));
	// ObjNConstraint_Type.push_back(0);

	// ObjNConstraint_Val.push_back((rE_opt[0] - Hand_max) * (rE_opt[0] - Hand_max));
	// ObjNConstraint_Type.push_back(0);
	//
	// ObjNConstraint_Val.push_back((rF_opt[0] - Hand_max) * (rF_opt[0] - Hand_max));
	// ObjNConstraint_Type.push_back(0);


				dlib::matrix<double,6,1> End_Effector_Dist;
				std::vector<int> End_Effector_Obs(6);

				End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

				dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Eqn_Maint_Matrix, Matrix_result;
				std::vector<double> sigma_temp, sigma_real;
				if(Opt_Type_Flag == -1)
				{
		// In this case, sigma_i has to be maintained.
					sigma_real = sigma_i;
				}
				else
				{
		// In other cases, sigma_i_child has to be maintained.
					sigma_real = sigma_i_child;
				}
				sigma_temp = Sigma2Pos(sigma_real, 0);			Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
				sigma_temp = Sigma2Pos(sigma_real, 1);			Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
				sigma_temp = Sigma2Vel(sigma_real);				Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);

	// 1. Active constraints have to be satisfied: Position and Velocity
				Matrix_result = Eqn_Pos_Matrix * End_Effector_Dist;
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

				Matrix_result = Eqn_Vel_Matrix * End_Effector_Vel;
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	// 2. Inactive constraints have to be strictly away from the obstacle
				dlib::matrix<double> ones_vector, temp_matrix;
				ones_vector = ONES_VECTOR_fn(6);
				Matrix_result = Ineqn_Pos_Matrix * (End_Effector_Dist - ones_vector * mini);
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

	// 3. Middle joints have to be strictly away from the obs
				temp_matrix = Middle_Joint_Obs_Dist_Fn(StateNDot_Init_i);
				Matrix_result = temp_matrix - ones_vector * mini;
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

	// 4. One more constraint to be added is to maintain the active unchanged constraint
				if(Opt_Type_Flag == -1)
				{
					Eqn_Maint_Matrix = Eqn_Maint_Matrix_fn(sigma_i, sigma_i);
				}
				else
				{
					Eqn_Maint_Matrix = Eqn_Maint_Matrix_fn(sigma_i, sigma_i_child);
				}
				Matrix_result = Eqn_Maint_Matrix * (End_Effector_Pos_ref - End_Effector_Pos);
				ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
				return;
			}

			int Sigma_Change(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child)
			{
	// This function is only used when a contact change is made for sure
				int Sigma_Index = 0;
				double Sigma_Result = 0;
				for (int i = 0; i < 4; i++)
				{
					Sigma_Result = sigma_i[i] - sigma_i_child[i];
					if(abs(Sigma_Result)==1)
					{
						Sigma_Index = i;
						break;
					}
				}
				return Sigma_Index;
			}
			dlib::matrix<double> Eqn_Maint_Matrix_fn(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child)
			{
	// This function is used to generate the contact maintenance matrix
				std::vector<double> sigma_maint;
				sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
				sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
				sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
				sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
				sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
				sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
				sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
				sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
				sigma_maint.push_back(sigma_i[2] * sigma_i_child[2]);
				sigma_maint.push_back(sigma_i[2] * sigma_i_child[2]);
				sigma_maint.push_back(sigma_i[3] * sigma_i_child[3]);
				sigma_maint.push_back(sigma_i[3] * sigma_i_child[3]);

				dlib::matrix<double> Eqn_Maint_Matrix;
				Eqn_Maint_Matrix = Diag_Matrix_fn(sigma_maint);
				return Eqn_Maint_Matrix;
			}
			std::vector<double> Sigma2Pos(std::vector<double> &sigma, int EqOrIneq)
			{
				std::vector<double> sigma_pos;
				if(EqOrIneq==0)
				{
					sigma_pos.push_back(sigma[0]);
					sigma_pos.push_back(sigma[0]);
					sigma_pos.push_back(sigma[1]);
					sigma_pos.push_back(sigma[1]);
					sigma_pos.push_back(sigma[2]);
					sigma_pos.push_back(sigma[3]);
				}
				else
				{
					sigma_pos.push_back(!sigma[0]);
					sigma_pos.push_back(!sigma[0]);
					sigma_pos.push_back(!sigma[1]);
					sigma_pos.push_back(!sigma[1]);
					sigma_pos.push_back(!sigma[2]);
					sigma_pos.push_back(!sigma[3]);
				}
				return sigma_pos;
			}
			std::vector<double> Sigma2Vel(std::vector<double> &sigma)
			{	std::vector<double> sigma_pos;
				sigma_pos.push_back(sigma[0]);
				sigma_pos.push_back(sigma[0]);
				sigma_pos.push_back(sigma[0]);
				sigma_pos.push_back(sigma[0]);
				sigma_pos.push_back(sigma[1]);
				sigma_pos.push_back(sigma[1]);
				sigma_pos.push_back(sigma[1]);
				sigma_pos.push_back(sigma[1]);
				sigma_pos.push_back(sigma[2]);
				sigma_pos.push_back(sigma[2]);
				sigma_pos.push_back(sigma[3]);
				sigma_pos.push_back(sigma[3]);
				return sigma_pos;
			}
			dlib::matrix<double> Diag_Matrix_fn(std::vector<double> &diag_vec)
			{
	// This function is used to generate a diagonal matrix within diag_vec to its diagonal elements
				int dim = diag_vec.size();
				dlib::matrix<double> Diag_Matrix;
				Diag_Matrix = dlib::zeros_matrix<double>(dim,dim);
				for (int i = 0; i < dim; i++)
				{
					Diag_Matrix(i,i) = diag_vec[i];
				}
				return Diag_Matrix;
			}
			dlib::matrix<double> Node_Expansion_fn(const Tree_Node &Node_i, int &Adjacent_Number)
			{
	// Checked! Oct.4th.2018 9:46PM
	// This function is used to conduct the node expansion for a given parent node
	// The basic consideration is not to have the flying-in-air phase
				std::vector<double> sigma_i = Node_i.sigma;
				dlib::matrix<double> Nodes_Sigma_Matrix;
				Nodes_Sigma_Matrix = dlib::zeros_matrix<double>(4,4);
				double sigma_sum;
				Adjacent_Number = 0;
				for (int i = 0; i < 4; i++)
				{
					std::vector<double> sigma_t = sigma_i;
					sigma_t[i] = !sigma_t[i];
					sigma_sum = sigma_t[0] + sigma_t[1] + sigma_t[2] + sigma_t[3];
					if(sigma_sum==0)
					{
						continue;
					}
					else
					{
						Nodes_Sigma_Matrix(Adjacent_Number, 0) = sigma_t[0];
						Nodes_Sigma_Matrix(Adjacent_Number, 1) = sigma_t[1];
						Nodes_Sigma_Matrix(Adjacent_Number, 2) = sigma_t[2];
						Nodes_Sigma_Matrix(Adjacent_Number, 3) = sigma_t[3];
						Adjacent_Number = Adjacent_Number + 1;
					}
				}
	// cout<<Nodes_Sigma_Matrix<<endl;
				return Nodes_Sigma_Matrix;
			}

			std::vector<double> End_RobotNDot_Extract(std::vector<double> &Opt_Soln, std::vector<double> &sigma, std::vector<double>&sigma_i_child)
			{
	// This function is used to extract the end robot state out from the Opt_Soln while conducting Impact Mapping if there exists
				int Opt_Type_Flag, Critical_Grid_Index;										double T_tot;
				dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;
				Opt_Seed_Unzip(Opt_Soln, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);
				Sigma_TransNGoal(sigma, sigma_i_child, Opt_Type_Flag, Critical_Grid_Index);

				std::vector<double> Robotstate_vec;
				if(Opt_Type_Flag ==1)
				{
		// In this case, there is impact mapping needed to be conducted
					Robot_StateNDot Robot_End_State;
					double Impulse_Mag;
					Robot_End_State = Impact_Mapping_fn(StateNDot_Traj, Impulse_Mag, sigma_i_child);
					Robotstate_vec = StateNDot2StateVec(Robot_End_State);
				}
				else
				{
					dlib::matrix<double> End_RobotstateDlib;
					End_RobotstateDlib = dlib::colm(StateNDot_Traj, Grids-1);
					for (int i = 0; i < End_RobotstateDlib.nr(); i++) {
						Robotstate_vec.push_back(End_RobotstateDlib(i));}
					}
					return Robotstate_vec;
				}
