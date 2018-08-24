#pragma once

#include "Eigen\Dense"
//#include "MathTool.h"
#include "Dynamics.h"
#include "Jacobian.h"
//#include "Kinematics.h"

#define epis 1e-5
#define SAMPLING_T_e3 1e-3

using namespace Eigen;

//define a Pose, and its derivatives.
struct PoseDeriv {

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	// --- position ---
	Vector3d p;  
	Vector3d pd;
	Vector3d pdd;

	VectorXd pr;

	// --- orientation ---
	Matrix3d R;
	Quaterniond qt; // quaternion of R
	Quaterniond qtd;
	Vector3d w;
	Vector3d wd;

	PoseDeriv() {
		p  = pd = pdd = w = wd = Vector3d::Zero();
		pr = VectorXd::Zero(6);
		
		qt	= Quaterniond::Identity();
		qtd = Quaterniond::Identity();
	}
};

//class Jacobian;
class Impedance
{
public:
	Impedance(void);
	~Impedance(void);


	void calculateForceTorq( double* out_Torq,	const PoseDeriv & PoseD,	const double* F_d,
							const double* N_d,  const double* q,	const double* qd,
							const double* F_ex,	const double* N_ex );
	// ==== main function ====

	VectorXd calculateForceTorq(const PoseDeriv & PoseD,const Vector3d F_d, const Vector3d N_d,	
								const VectorXd & q,		const VectorXd & qd, const Vector3d & F_ex,
								const Vector3d & N_ex); //direct force control

	void gravityCompensation(double* out_Torque,const double* q);

	
	void Impedance::initialzero();//fu
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	int dof;

	//Kin_6Axis* kin; ///Kinematics*	kin

	Jacobian* jcb;
	Dynamics* dyn;
	Quaterniond quate;//endeffector orientation

	// outer loop gain.  use(Recommended Starting Value)
	Matrix3d Mp; // translational	Inertia Matrix		
	Matrix3d Mo; // rotational		Inertia matrix	

	Matrix3d Dp; // translational	damping matrix		
	Matrix3d Do; // rotational		damping matrix

	Matrix3d Kp; // translational	stifness matrix	
	Matrix3d Ko; // rotational		stifness matrix		
	
//~~~~~~~~~~iner motion loop gains. starting value(adjust based or particular manipulator environment interaction)
	Matrix3d K_Pp;
	Matrix3d K_Po;
	
	Matrix3d K_Vp; 
	Matrix3d K_Vo;  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	////stiffeness gains
	//Matrix3d kpp,kpO;
	//MatrixXd kv;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//direct force control gains
	Matrix3d K_Af,K_Vf;
	Matrix3d K_Am,K_Vm;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	MatrixXd M;
	VectorXd C;
	VectorXd G;

	Matrix3d Identity_3x3;
	MatrixXd matJ;  // = jcb->getJ(q);
	MatrixXd _matJInv;
	MatrixXd matJd; // = jcb->getJd(q, qd);

	//double half_sampling_T;

	// temperary variables to convert user-input 'double[]' to 'Eigen style'.  And convert output back to double[].
	VectorXd p_Vec,q_Vec, qd_Vec, torq_Vec,torq_Vec2,q_Vec2;
	Vector3d F_ext_Vec, N_ext_Vec,fd_Vec,nd_Vec;

	PoseDeriv PoseC, PoseCD;

	Vector3d wcd,wc,wcd_prev;

	//---------- old static parameters in calculateForceTorq() --------
	PoseDeriv PoseC_prev;
	PoseDeriv PoseE; //End Effector Pose (from real-world manupulator)
	VectorXd f_accs;
	VectorXd f_qdd; 
	VectorXd f_torque;

	// ==== sub function ====
	void		forcePose		(PoseDeriv &PoseC_prev, const PoseDeriv &PoseD, const Vector3d &F_d,  const Vector3d &N_d, const Vector3d &F_ex, const Vector3d &N_ex); //force control
	PoseDeriv	FKI				(const VectorXd &q, const VectorXd &qd); //temp
	VectorXd	forcePoseControl	(const PoseDeriv &PoseD,const PoseDeriv &PoseE);

	// ===== tool functions ======
	Quaterniond quatDerivative	(const Vector3d &w, Quaterniond &qt);
	void R2Quaternion(const Matrix3d R);
	Matrix3d Impedance::skew(const Vector3d &v);
	Vector3d	integrate		(const Vector3d &past, const Vector3d &present);
	short		integrate		(Quaterniond &qt, const Quaterniond &past, const Quaterniond &present);
};



