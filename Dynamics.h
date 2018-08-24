#pragma once

#include "Eigen\Dense"
//#include "MathTool.h"
#include "DHTable.h"
//#include "Jacobian.h"
//#include "Log.h"
#define DEG2RAD 0.01745329251994
#define RAD2DEG 57.2957795130823
#define PI 3.1415926535897932384

using namespace Eigen;

class Dynamics
{
public:
	Dynamics(int dof);
	~Dynamics(void);

	//Inverse Dynamics:: get joint torque based on Newton Euler recursive Method.
	void     inv(double* out_Torq,
				 const double* q,	    const double* qd,    const double* qdd,       //Note: assume 6 axis.
				 const double F_ex[3], const double N_ext[3]);

	VectorXd inv(const VectorXd q,	   const VectorXd qd,    const VectorXd qdd,
				 const Vector3d F_ex, const Vector3d N_ext);
	VectorXd g(const VectorXd q);
	VectorXd auxilaryTorque(const VectorXd qd);
	/*void invKin();*/

	Matrix4d get_A( double theta, const int Axis_No);
	Matrix4d get_A( int Axis_No);
	void cal_allA(const VectorXd q);

	int sign(double value);

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	DYNTable DYN;

	int dof;
	Matrix3d Identity_3x3;
	Vector3d Zero_3;
	Vector3d z0;		// Z not. = [0,0,1]T

	Matrix3d *Roi;		// Rotation Matrix from root to flange i.
	Matrix3d *Rt;		// Rotation Matrix from i to i-1.
	Matrix4d *T;		// Trans    Matrix from i-1 to i.
	Vector3d *bi;
	Vector3d *P;		// Position Vector from i-1 to i.
	Vector3d *Omg;		// Omega	Vector i.
	Vector3d *Omg_d;	// Omega dot Vector i.
	Vector3d *Ae, *Ac;	// Acc		Vector i at end  /center of link i.
	Vector3d *Fs, *Fc;	// Force	Vector i at start/center of link i.
	Vector3d *Ns, *Nc;	// Torque	Vector i at start/center of link i.

	// temperary variables to convert user-input 'double[]' to 'Eigen style'.  And convert output back to double[].
	VectorXd q_Vec, qd_Vec, qdd_Vec, torq_Vec;
	Vector3d F_ex_Vec, N_ex_Vec;

	Matrix4d *matA; // =getA(q[i-1], i);

	//Matrix3d Ro;
	//
	//VectorXd  *f/*joint force*/,*n/*joint troque*/,
	//			/, F_ex,N_ex;
	//Vector3d z0, gravity;

	//VectorXd _q,_qp,_qpp;
	//double _theta[7];



	//int _linkNum,_dof;
	//int _jointNum;


	//int joint_type;
	//double _a,_d,_alpha,joint_offset,_mass;
	//double _Im,_Gr,_B,_Cf; //joint & motor values
	//VectorXd rc, p;
	//Matrix3d I;
};

