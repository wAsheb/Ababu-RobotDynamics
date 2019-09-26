#include "Dynamics.h"
#include <iostream>
using namespace std;



Dynamics::Dynamics(int dof)
{
	this->dof	= dof;
	Identity_3x3= Matrix3d::Identity();
	Zero_3		= Vector3d::Zero();
	z0		= Vector3d(0,0,1);

	Roi		= new Matrix3d[dof+1];
	Rt		= new Matrix3d[dof+1];
	T		= new Matrix4d[dof+1];
	bi		= new Vector3d[dof+1];
	P		= new Vector3d[dof+1];
	Omg		= new Vector3d[dof+1];
	Omg_d	        = new Vector3d[dof+1];
	Ae		= new Vector3d[dof+1];
	Ac		= new Vector3d[dof+1];
	Fs		= new Vector3d[dof+1];
	Fc		= new Vector3d[dof+1];
	Ns		= new Vector3d[dof+1];
	Nc		= new Vector3d[dof+1];
	matA	        = new Matrix4d[dof+1];

	q_Vec	= VectorXd(dof);
	qd_Vec	= VectorXd(dof);
	qdd_Vec	= VectorXd(dof);
	torq_Vec= VectorXd(dof);

	/*Roi[0]	= Identity_3x3;
	T[0]	= Matrix4d::Identity();
	P[0]	= Zero_3;
	Omg[0]	= Zero_3;
	Omg_d[0]= Zero_3;
	
	Ac[0]	= Zero_3;*/
	
	
	for(int i = 0; i< dof+1; i++) {

		Roi	[i]	= Matrix3d::Identity();
		Rt	[i]	= Matrix3d::Identity();
		T	[i]	= Matrix4d::Identity();
		bi	[i] = Zero_3;
		P	[i]	= Zero_3;
		Omg	[i]	= Zero_3;
		Omg_d[i] =Zero_3;
		Ae	[i]	= Zero_3;
		Ac	[i]	= Zero_3;
		Fs	[i]	= Zero_3;
		Fc	[i]	= Zero_3;
		Ns	[i]	= Zero_3;
		Nc	[i]	= Zero_3;
		matA[i]	= Matrix4d::Identity();
	}

	Ae[0]	= Vector3d(0,0,9.81);
	
}
Dynamics::~Dynamics(void)
{
}

void Dynamics::inv(
	double* out_Torq,
	const double* q,	   const double* qd,    const double* qdd,       
	const double F_ex[3], const double N_ex[3])
{
	// convert double[] to Eigen.
	for(int i = 0; i< dof; i++)
	{
		q_Vec(i)	= q[i];
		qd_Vec(i)	= qd[i];
		qdd_Vec(i)	= qdd[i];
	}
	for(int i = 0; i< 3; i++)
	{
		F_ex_Vec(i) = F_ex[i];
		N_ex_Vec(i) = N_ex[i];
	}

	// call algorithm.
	torq_Vec = inv(q_Vec, qd_Vec, qdd_Vec, F_ex_Vec, N_ex_Vec);

	// convert Eigen to double[].
	for(int i = 0; i< dof; i++)
	{
		out_Torq[i]	= torq_Vec(i);
	}
}

VectorXd Dynamics::inv(
	const VectorXd q,     const VectorXd qd,   const VectorXd qdd,
	const Vector3d F_ex, const Vector3d N_ex)
{
	
	
	/*Matrix4d T_flange;
	T_flange<<1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1;
	Matrix3d R_flange;
	R_flange<<1,0,0,0,-1,0,0,0,-1;*/
	Matrix3d tempRt;

	for(int i=1; i <= dof; i++){
		
		tempRt	= T[i].block<3,3>(0,0);

		//if(i==dof){
		//	T[i]	= get_A(q(i-1),i)/**T_flange*/;
		//	Rt[i]	= (tempRt/**R_flange*/).transpose();
		//	Roi[i]	*=tempRt/**R_flange*/;
		//}
		T[i]	= get_A(i);//get_A(q(i-1),i);
		Rt[i]	= tempRt.transpose();
		Roi[i]  *= tempRt;

		bi[i] = Roi[i].transpose()*Roi[i-1]*z0;

		//Roi[i]	= Roi[i-1] * T[i].block<3,3>(0,0);
		P[i]	= T[i].block<3,1>(0,3);		
	}

	
	/*static int count = 0;
	if(count++ >100) count = 0;
	if(count == 0)
		Log::write(3, "bi[i],%lf, %lf,%lf,\n ", bi[2](0), bi[2](1),bi[2](2));*/

	// joint angular velocity & acceleration, and acceleration of center of mass
	Vector3d v1, v2, v3, v4;

	for(int i = 1; i <= dof; i++){
		
		//-----------------------------------------
		v1 = Rt[i]*Omg[i-1];
		v2 = Rt[i]*Omg_d[i-1];

		v3 = qd(i-1) * z0;//qd(i-1) * bi[i];
		v4 = qdd(i-1)* z0;//qdd(i-1)* bi[i];

		Omg[i]   = v1 + v3;
		Omg_d[i] = v2 + v1.cross(v3) + v4;

		//Omg[i]   = Rt[i]*Omg[i-1]	 + qd(i-1)*z0;
		//Omg_d[i] = Rt[i]*Omg_d[i-1] + (Rt[i]*Omg[i-1]).cross(z0*qd(i-1)) + qdd(i-1)*z0;
		//--------------------------------------------
		
		Vector3d tempV = /*P[i] -*/ DYN.LinkMassCenter(i);

		Ae[i] = Rt[i]*(Ae[i-1] +  Omg_d[i-1].cross(P[i])    + Omg[i-1].cross(Omg[i-1].cross(P[i])));							
		Ac[i] = Ae[i]		   +  Omg_d[i].cross(tempV)		+ Omg[i].cross(Omg[i].cross(tempV));

		

	}
		
	// backwards iteration starting from last link(TCP) to base joint
	VectorXd out_torque = VectorXd::Zero(6);
	MatrixXd temp;
	Matrix3d t1;

	VectorXd auxilary(6);
	auxilary = auxilaryTorque(qd);

	for(int i = dof; i >= 1; i--) {

		Fc[i] = Ac[i]*DYN.LinkMass(i);
		Nc[i] = DYN.LinkI(i)*Omg_d[i] + Omg[i].cross( DYN.LinkI(i)*Omg[i]);
		

		if(i == dof) {
			Fs[i] = Fc[i] + F_ex;			
			Ns[i] = DYN.LinkMassCenter(i).cross(Fs[i]) + Nc[i] + N_ex;		
		} 
		else {			
			t1 = T[i+1].block<3,3>(0,0);
			Fs[i] = t1*Fs[i+1] + Fc[i];			
			Ns[i] = t1*Ns[i+1] + (P[i+1]).cross(t1*Fs[i+1]) + DYN.LinkMassCenter(i).cross(Fc[i]) + Nc[i]; 
			//Fs[i] = T[i+1].block<3,3>(0,0)*Fs[i+1] + Fc[i];			
			//Ns[i] = T[i+1].block<3,3>(0,0)*Ns[i+1] + (P[i+1]).cross(T[i+1].block<3,3>(0,0)*Fs[i+1]) + DYN.LinkMassCenter(i).cross(Fc[i]) + Nc[i]; 
		}
		

				temp = (z0.transpose()*Ns[i]);
				
				double jgr = DYN.JointGearRatio(i);
				
				/*int SIGN;
				
				if(qd(i-1) > 0){
					
					SIGN = 1;
				
				}if (qd(i-1) < 0){
					
					SIGN = -1;
				}else
					SIGN = 0;*/
				
		out_torque(i-1) = temp(0,0)
						+ DYN.MotorInertia(i)*jgr*jgr*qdd(i-1)
						+ jgr*(jgr*DYN.JointDampingB(i)*qd(i-1) + DYN.JointFrictionCf(i)*sign(qd(i-1))); /*+ auxilary(i-1) + DYN.StartFriction(i)*sign(qd(i-1))*/
	}
	
	return out_torque;
}
//this emperical torque to account unmodeled effect such as friction at joints
VectorXd Dynamics::auxilaryTorque(const VectorXd qd){

	VectorXd auxilaryAid = VectorXd::Zero(6);
	
		double a_Factor[6]				= { 0.1,0.1,0.1,0.1,0.1,0.1}; //a
		double K_Factor[6]				= {1000,100,100,1000,100,1000}; //k
		double aidTorqueAxis[6]			= {300,300,300,300,300,300}; //J
		double frictionPositive[6]		= {50,50,50,50,50,50};
		double frictionNegative[6]      = {-50,-50,-50,-50,-50,-50};//

		//T_aid = T0*(w/(aw^2+K)) w -angular speed

		for(int i = 0; i < 6; i++)
		{
		auxilaryAid(i) = aidTorqueAxis[i] * qd(i) / (a_Factor[i] * qd(i) * qd(i) + K_Factor[i]);
		}
		return auxilaryAid;

}
//this is to considering only the gravity term of total torque at each joint
//call this in hand guide mode
VectorXd Dynamics::g(const VectorXd q){

	/*Matrix4d T_flange;
	T_flange<<1,0,0,0,
			  0,-1,0,0,
			  0,0,-1,0,
			  0,0,0,1;
	Matrix3d R_flange;
	R_flange<<1,0,0,	
			  0,-1,0,
			  0,0,-1;*/
	Matrix3d tempRt;

	for(int i=1; i <= dof; i++){
		
		tempRt	= T[i].block<3,3>(0,0);

		//if(i==dof){
		//	T[i]	= get_A(q(i-1),i)/**T_flange*/;
		//	Rt[i]	= (tempRt/**R_flange*/).transpose();
		//	Roi[i]	*=tempRt/**R_flange*/;
		//}
		T[i]	= get_A(i);//get_A(q(i-1),i);
		Rt[i]	= tempRt.transpose();
		Roi[i]  *= tempRt;

		bi[i] = Roi[i].transpose()*Roi[i-1]*z0;

		//Roi[i]	= Roi[i-1] * T[i].block<3,3>(0,0);
		P[i]	= T[i].block<3,1>(0,3);		
	}

	/*static int count = 0;
	if(count++ >100) count = 0;
	if(count == 0)
		Log::write(3, "bi[i],%lf, %lf,%lf,\n ", bi[2](0), bi[2](1),bi[2](2));*/

	// joint angular velocity & acceleration, and acceleration of center of mass
	//Vector3d v1, v2, v3, v4;

	for(int i = 1; i <= dof; i++){
		
	//-----------------------------------------
		/*v1 = Rt[i]*Omg[i-1];
		v2 = Rt[i]*Omg_d[i-1];

		v3 = qd(i-1) *bi[i];//z0;
		v4 = qdd(i-1)*bi[i];//z0;

		Omg[i]   = v1 + v3;
		Omg_d[i] = v2 + v1.cross(v3) + v4;

		Omg[i]   = Rt[i]*Omg[i-1]	 + qd(i-1)*z0;
		Omg_d[i] = Rt[i]*Omg_d[i-1] + (Rt[i]*Omg[i-1]).cross(z0*qd(i-1)) + qdd(i-1)*z0;
		*/
		
		//--------------------------------------------
		
		Vector3d tempV = /*P[i] -*/ DYN.LinkMassCenter(i);

		Ae[i] = Rt[i]*(Ae[i-1]) ;							
		Ac[i] = Ae[i]		   ;

		

	}
	
	
	// backwards iteration
	VectorXd out_torque = VectorXd::Zero(6);
	MatrixXd temp;
	Matrix3d t1;

	
	for(int i = dof; i >= 1; i--) {

		Fc[i] = Ac[i]*DYN.LinkMass(i);
		//Nc[i] = DYN.LinkI(i)*Omg_d[i] + Omg[i].cross( DYN.LinkI(i)*Omg[i]);
		

		if(i == dof) {
			Fs[i] = Fc[i] /*+ F_ex*/;			
			Ns[i] = DYN.LinkMassCenter(i).cross(Fs[i]) /*+ Nc[i] + N_ex*/;		
		} 
		else {			
			t1 = T[i+1].block<3,3>(0,0);
			Fs[i] = t1*Fs[i+1] + Fc[i];			
			Ns[i] = t1*Ns[i+1] + (P[i+1]).cross(t1*Fs[i+1]) + DYN.LinkMassCenter(i).cross(Fc[i]) /*+ Nc[i]*/; 
			//Fs[i] = T[i+1].block<3,3>(0,0)*Fs[i+1] + Fc[i];			
			//Ns[i] = T[i+1].block<3,3>(0,0)*Ns[i+1] + (P[i+1]).cross(T[i+1].block<3,3>(0,0)*Fs[i+1]) + DYN.LinkMassCenter(i).cross(Fc[i]) + Nc[i]; 
		}
		
		temp = (z0.transpose()*Ns[i]);
		
		out_torque(i-1) = temp(0,0);
	
	}
	
	return out_torque;

}
//Transformation matrix from joint (i-1) to i
Matrix4d Dynamics::get_A( int Axis_No)
{
	return matA[Axis_No];
}
void Dynamics::cal_allA(const VectorXd q)
{
	for(int i = 1; i<=dof; i++)
		matA[i] = get_A(q[i-1], i);
}
//helper
int Dynamics::sign( double value)
{
	if (value > 0)
	
		return 1;
	
	else if (value < 0)
	
		return -1;
	
	else 
		return 0;
}
//static define 6-axis robot arm used for testing. change it to your robot physical property 
double d[7] = {-0.423, 0.0,0.0, -0.284, 0.0,0.2305,0};// index 6 not used d6 =0.0905
double a[7] = {0,0.105,0.2800,0.065,0.0,0.0,0.0};// index 6 not used
double alpha[7] = {PI, PI/2.0, 0.0, PI/2.0, -PI/2.0, -PI/2.0, 0.0}; //is alpha[5] is - / +

Matrix4d Dynamics::get_A( double theta, const int Axis_No) {
		
	double cosA = cos( alpha[Axis_No-1] );
	double sinA = sin( alpha[Axis_No-1] );
	double cosT = cos( theta );
	double sinT = sin( theta );
	double D = d[Axis_No-1];
	double A = a[Axis_No-1];

	

	Matrix4d result;
	result << cosT,			-sinT,		0,	 	A,
			sinT*cosA,	cosT*cosA,	-sinA,	-sinA*D,
			sinT*sinA,	cosT*sinA,	cosA,	cosA*D,
			0,			0,			0,		1;

	return result;	
}

   
