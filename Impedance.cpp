#include "Impedance.h"
#include <iostream>
using namespace std;

//Matrix4d get_A2( double theta, const int Axis_No);

Impedance::Impedance(void)
{
	dof = 6;

	//kin				= new Kin_6Axis(0);
	jcb				= new Jacobian(dof);
	dyn				= new Dynamics(dof);
	jcb->dyn = dyn; //only to let jcb get "getA(i)" result.

	//half_sampling_T = 0.5*SAMPLING_T_e3;
	Identity_3x3	= Matrix3d::Identity();

	f_accs = VectorXd(6);
	f_qdd  = VectorXd(6);
	f_torque= VectorXd(6);

//~~~~~~~~~~~~~~~~~~~~~~~~~~ Outer Loop(Impedance) gain~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Mp		= Identity_3x3 * 6;
	Mo		= Identity_3x3 * 0.4;

	Dp<<200,0,0,   0,200,0,  0,0,2000; //???

	Do<<3.5,0,0, 0,3.5,0, 0,0,3.5; //???	


	Kp<<100,0,0,  0,100,0,   0,0,50;

	Ko<<2.5,0,0,  0,2.5,0,   0,0,2.5;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Inner Loop gain (motion) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//usually characteristics of the robot i.e. don't change once tuned 
	K_Pp<< 2500,0,0	 ,0,2275,0	,0,0,2375;//4950//2025;
	K_Po<< 2600,0,0,	0,2400,0,	0,0,2500;//4500


	K_Vp<< 25,0,0,	0,25,0,	0,0,25;
	K_Vo<< 5,0,0,	0,5,0,	0,0,5;

	p_Vec	= VectorXd(dof);
	q_Vec	= VectorXd(dof);
	q_Vec2	= VectorXd(dof);  // added for direct force control (maps egine vector to double vector)
	qd_Vec	= VectorXd(dof);
	torq_Vec= VectorXd::Zero(dof);

	wcd		 = Vector3d::Zero();
	wc		 = Vector3d::Zero();
	wcd_prev = Vector3d::Zero();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~stiffness gains~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//kpp = kpO = Matrix3d::Identity();
	//kpp(0,0) = 250;kpp(1,1) = 190; kpp(2,2) = 750;		// <<350,0,0,		0,190,0,	0,0,550;
	//kpO(0,0) = 2.1;kpO(1,1) = 2.1;kpO(2,2) = 2.1;	// <<2.1,0,0,		0,2.1,0,	0,0,2.1;

	//kv = MatrixXd::Identity(6,6);
	//kv(0,0) = 100;			kv(1,1) = 76;		kv(2,2) = 300;
	//kv(3,3) = 25;			kv(4,4) = 15;		kv(5,5) = 35;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Direct force Control gain s~~~~~~~~~~~~~~~~~~~~
	K_Af << 8,0,0,	0,8,0,	0,0, 8;
	K_Vf << 200,0,0,	0,200,0,	0,0, 100;

	K_Am <<1,0,0,	0,1,0,	0,0,1;
	K_Vm <<	85,0,0,	0,85,0,	0,0,85;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
}
Impedance::~Impedance(void)
{
	//delete kin;
	delete jcb;
	delete dyn;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Direct Force control ~~~~~~~~~~~~~~~~~~~~~~~~~
void Impedance::calculateForceTorq( double* out_Torq,
							 const PoseDeriv & PoseD,	const double* F_d, const double* N_d,
							 const double* q,		const double* qd,     const double* F_ex,
							 const double* N_ex )
{
	// convert double[] to Eigen.
	for(int i = 0; i< 6; i++)
	{
		q_Vec(i)	= DEG2RAD*q[i];
		qd_Vec(i)	= DEG2RAD*qd[i];
	}
	for(int i = 0; i< 3; i++)
	{
		F_ext_Vec(i) = F_ex[i];
		N_ext_Vec(i) = N_ex[i];

		fd_Vec(i) = F_d[i];
		nd_Vec(i) = N_d[i];
	}
	/*static int count = 0;
	if(count++ >100) count = 0;
	if(count == 0)
		Log::write(1, "q,qd,%lf, %lf,\n", q_Vec(2),  qd_Vec(2));*/
	// call algorithm.

	dyn->cal_allA(q_Vec);
	jcb->calJ(q_Vec);
	torq_Vec = calculateForceTorq(PoseD, fd_Vec, nd_Vec, q_Vec, qd_Vec, F_ext_Vec, N_ext_Vec);

	// convert Eigen to double[].
	for(int i = 0; i< dof; i++)
		out_Torq[i]	= torq_Vec(i);
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End Direct Force~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//gravity compensation 
void Impedance::gravityCompensation(double* out_Torque,const double* q){

	for(int i=0;i<6;i++)
		q_Vec2(i) = DEG2RAD*q[i];
	
	//torq_Vec2 = getG(q_Vec2);
	dyn->cal_allA(q_Vec2);

	torq_Vec2 = dyn->g(q_Vec2); //TODO choose between these two

	for(int i = 0; i<6; i++)
		out_Torque[i]	= torq_Vec2(i);
	
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	START Force Control Block~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START Direct Force Control Block ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VectorXd Impedance::calculateForceTorq(
	const PoseDeriv & PoseD,	const Vector3d F_d,		const Vector3d N_d,
	const VectorXd & q,			const VectorXd & qd,	const Vector3d & F_ex,	
	const Vector3d & N_ex)
{

	// 量時間
	/*LARGE_INTEGER timeStart;	
	RtGetClockTime(CLOCK_2,&timeStart);*/

	_matJInv  = jcb->DLSJacobInverse(q);   // use this for robust result near singular configuration 
	matJd = jcb->getJd(q,qd);

	
	// 量時間
	/*LARGE_INTEGER timeEnd1;
	double calTime1;
	RtGetClockTime(CLOCK_2,&timeEnd1);
	calTime1 = double(timeEnd1.QuadPart - timeStart.QuadPart)/(10);*/


	//compliencePose(PoseC_prev, PoseD, F_ex, N_ex);
	PoseE = FKI(q, qd);


	// 量時間
	/*LARGE_INTEGER timeEnd2;
	double calTime2;
	RtGetClockTime(CLOCK_2,&timeEnd2);
	calTime2 = double(timeEnd2.QuadPart - timeEnd1.QuadPart)/(10);*/
	

	forcePose(PoseC_prev,PoseD,F_d,N_d,F_ex, N_ex);

	// 量時間
	/*LARGE_INTEGER timeEnd3;
	double calTime3;
	RtGetClockTime(CLOCK_2,&timeEnd3);
	calTime3 = double(timeEnd3.QuadPart - timeEnd2.QuadPart)/(10);*/

	f_accs  = forcePoseControl(PoseD,PoseE); 
	
	// 量時間
	/*LARGE_INTEGER timeEnd4;
	double calTime4;
	RtGetClockTime(CLOCK_2,&timeEnd4);
	calTime4 = double(timeEnd4.QuadPart - timeEnd3.QuadPart)/(10);*/


	// *** qdd = matJ.inverse()*(accs - matJd * qd); //***
	//MatrixXd matJInv = matJ.inverse();      // unstable near singularity use DLS inverse of the Jacobian
	VectorXd vec1 = matJd * qd; 
	vec1 = f_accs - vec1; 
	f_qdd = _matJInv*vec1;
	
	
//==================== control law =====================================================================================================
// torque = B(q)*jacobian.inverse()*(a - jacobian_dot*dq/dt) + c(q,dq/dt)*dq/dt + d(q,dq/dt) + g(q) + jacobian.transpose()*h  --------(1)
// in place of equation (1) use torque calculator in dynamics class
//======================================================================================================================================
	
	
	// 量時間
	/*LARGE_INTEGER timeEnd5;
	double calTime5;
	RtGetClockTime(CLOCK_2,&timeEnd5);
	calTime5 = double(timeEnd5.QuadPart - timeEnd4.QuadPart)/(10);*/

	f_torque = dyn->inv(q, qd, f_qdd, F_ex, N_ex);

	// 量時間
	/*LARGE_INTEGER timeEnd6;
	double calTime6;
	RtGetClockTime(CLOCK_2,&timeEnd6);
	calTime6 = double(timeEnd6.QuadPart - timeEnd5.QuadPart)/(10);	*/

	//static int count = 0;
	//count++;
	////if(count++ >500)count = 0;
	//if(count< 50){
	//Log::write(1, "Time:calculateTorq: ,%lf, %lf, %lf, %lf, %lf, %lf", calTime1, calTime2, calTime3, calTime4, calTime5, calTime6);
	//Log::writeln(1, ", TotalTime, %lf",calTime1 + calTime2 + calTime3 + calTime4 + calTime5 + calTime6);
	//}
	
	
	return f_torque;

}
	
void Impedance::initialzero()
{//force control
		PoseC.p  = Vector3d::Zero();
		PoseC.pd = Vector3d::Zero();
		PoseC.pdd = Vector3d::Zero();

		PoseC_prev.p = Vector3d::Zero();
		PoseC_prev.pd = Vector3d::Zero();
		PoseC_prev.pdd = Vector3d::Zero();
		//Log::writeln(2, "zero");

} 
void Impedance::forcePose(PoseDeriv &PoseC_prev,const PoseDeriv &PoseD,const Vector3d &F_d,const Vector3d &N_d,const Vector3d &F_ex, const Vector3d &N_ex) //force control
{
	//K_Vf(2,2) = shm->bounceK;
	//static int countt=0;
	//if(countt++ >3000) countt = 0;
	//if(countt == 0)
	//	Log::writeln(0, "kv = %lf\n", K_Vf(2,2));

	//Quaterniond qt;
	//static bool first = true;
	//if(first)
	//{
	//	PoseC.p			= PoseD.p;
	//	PoseC.qt		= PoseD.qt;
	//	PoseC.qtd		= quatDerivative(PoseC.w, PoseC.qt);
	//	PoseC_prev.qtd	= PoseC.qtd;

	//	PoseC.wd		= PoseD.wd;
	//	PoseC.w			= PoseD.w;
	//	first = false; 
	//}

	// set new value //20170829
	/*K_Af(3,3) = shm->afValue;
	K_Vf(3,3) = shm->vfValue; */

	//if(PoseC.qt.dot(PoseD.qt) < 0){ // if <0, multiply PoseD.qt
	//	
	//		qt.w()	 = PoseD.qt.w()*-1;
	//		qt.vec() = PoseD.qt.vec()*-1;

	//}
	//else {
	//	
	//	qt = PoseD.qt;
	//}
	//Qdc = Qc.inverse()*Qd  quaternion based endeffector orientation
	// use PoseCD for consistance with Impedance case, otherwise in actual PoseDC(PoseD - PoseC) is computed

	PoseC_prev = PoseC;

	PoseCD.qt.w()	= PoseD.qt.w()*PoseC.qt.w() + (PoseC.qt.vec().transpose())*PoseD.qt.vec();
	PoseCD.qt.vec() = PoseC.qt.w()*PoseD.qt.vec() - PoseD.qt.w()*PoseC.qt.vec() - skew(PoseC.qt.vec())*PoseD.qt.vec();
	PoseCD.w	= PoseD.w	- PoseC.w;
	PoseCD.wd	= PoseD.wd	- PoseC.wd;
	/*
	//force control/ trajectory

	PoseF.pdd = K_Af.inverse()*((F_d - F_ex) - K_Vf*PoseF.pd);
	PoseF.pd = PoseF_prev.pd + integrate(PoseF_prev.pdd, PoseF.pdd);  
	PoseF.p += integrate(PoseF_prev.pd, PoseF.pd); 
	*/

	PoseC.pdd = K_Af.inverse()*((F_d - F_ex) - K_Vf*PoseC.pd);//PoseC.pdd = K_Af.inverse()*((F_d - F_ex) - K_Vf*PoseC.pd);
	PoseC.pd = PoseC_prev.pd + integrate(PoseC_prev.pdd,PoseC.pdd);
	PoseC.p += integrate(PoseC_prev.pd,PoseC.pd);


	//moment control/trajectory
	/*
	PoseF.wd = K_Am.inverse()*((N_d - N_ex) - (K_Vm*PoseF.w)); //PoseF = PoseC -- compliant frame
	PoseF.w  += integrate(PoseF_prev.wd, PoseF.wd); //check consistency on refrence frame used
	PoseF.qtd = quatDerivative(PoseF.w, PoseF.qt);  //note: qtd is not normalized.			
	integrate(PoseF.qt, PoseF_prev.qtd, PoseF.qtd);
	*/
	 
	Matrix3d RCTranspose = PoseC.R.transpose();
	Vector3d wc		= RCTranspose*PoseC.w;

	Vector3d wc_d = K_Am*((N_d - N_ex) - (K_Vm*wc));//Vector3d wc_d = K_Am.inverse()*((N_d - N_ex) - (K_Vm*wc));
	PoseC.wd = PoseC.R*wc_d;

	PoseC.w  += integrate(PoseC_prev.wd, PoseC.wd);

	PoseC.qtd = quatDerivative(PoseC.w, PoseC.qt);
	integrate(PoseC.qt, PoseC_prev.qtd, PoseC.qtd);

	//to prevent the jittering at start and end of new motion


	//PoseF_prev = PoseF;
	//PoseC_prev = PoseC;

	/*static int countt = 0;
	if(countt++ >3000) countt = 0;
	if(countt == 0)*/
	//Log::writeln(2, "last:,%lf, %lf,%lf",PoseC.p(2),PoseC.pd(2),PoseC.pdd(2));

}
VectorXd Impedance::forcePoseControl(const PoseDeriv &PoseD,const PoseDeriv &PoseE)
{
	PoseDeriv Pose_r; //force refernce pose
	Quaterniond qt_temp;
	Vector3d quat_v_error;
	R2Quaternion(PoseE.R);
	//PoseE.qt = Quaterniond(PoseE.R);

	/*Quaterniond qt;*/
	//static bool first = true;
	//if(first)
	//{
	//	Pose_r.p		= PoseE.p;
	//	Pose_r.qt		= PoseE.qt;
	//	Pose_r.qtd		= quatDerivative(Pose_r.w, Pose_r.qt);
	//	/*Pose_r_prev.qtd	= Pose_r.qtd;*/

	//	Pose_r.wd		= PoseE.wd;
	//	Pose_r.w		= PoseE.w;
	//	first = false; 
	//}
	
	if(/*quatDotProduct(PoseE.qt,PoseC.qt)*/PoseE.qt.dot(Pose_r.qt) < 0){
		qt_temp.w()	= Pose_r.qt.w()*-1;
		qt_temp.vec() = Pose_r.qt.vec()*-1;
	}
   else
	  qt_temp = Pose_r.qt;
	
	//qt_temp.vec() = Vector3d::Zero();
		
	
	

	//parallel composition
	Pose_r.p	= PoseC.p	+ PoseD.p;
	Pose_r.pd	= PoseC.pd	+ PoseD.pd;
	Pose_r.pdd	= PoseC.pdd	+ PoseD.pdd;

	//Log::writeln(2, "first:,%lf, %lf,%lf",PoseC.p(2),PoseC.pd(2),PoseC.pdd(2));


	Pose_r.qt = PoseC.qt*PoseCD.qt ;//(PoseD.qt.inverse()*PoseC.qt);
	Pose_r.w  = PoseC.w  + PoseCD.w;//(PoseD.w - PoseC.w);
	Pose_r.wd = PoseC.wd + PoseCD.wd;//(PoseD.wd - PoseC.wd);

	quat_v_error = PoseE.R*((PoseE.qt.inverse()*Pose_r.qt).vec()); //(PoseE.qt.w()*qt_temp.vec() - qt_temp.w()*PoseE.qt.vec() - skew(PoseE.qt.vec())*qt_temp.vec());//quaternion diffrence 

	Vector3d ap, ao;

	VectorXd a(6);

	ap = Pose_r.pdd + K_Vp*(Pose_r.pd - PoseE.pd) + K_Pp*(Pose_r.p - PoseE.p); 
	ao = Pose_r.wd  + K_Vo*(Pose_r.w - PoseE.w)   + K_Po*(quat_v_error);//PoseE.R*qt_CE.vec();
	
	
	a.block<3,1>(0,0) = ap;
	a.block<3,1>(3,0) = ao;	

	static int countt = 0;
	if(countt++ >3000) countt = 0;
	if(countt == 0)
		Log::writeln(0, "apao:,%lf, %lf,%lf,%lf,%lf,%lf",ap(0),ap(1),ap(2),ao(0),ao(1),ao(2));
	
	//Log::writeln(3, "apao:,%lf, %lf,%lf,%lf,%lf,%lf",ap(0),ap(1),ap(2),ao(0),ao(1),ao(2));
	return a;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  END Force Control Block ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



PoseDeriv Impedance::FKI(const VectorXd& q, const VectorXd& qd)
{
	PoseDeriv PoseE;

	
	Matrix4d T = Matrix4d::Identity();

	for(int i = 1; i<= 6; i++)
		T *= dyn->get_A(i); //get_A2(q(i-1), i);


	PoseE.R  = T.block<3,3>(0,0);
	PoseE.p  = T.block<3,1>(0,3);
	
	VectorXd Ve = VectorXd::Zero(6);
	Ve			= jcb->Jmat*qd;
	PoseE.pd	= Ve.topRows(3);//(jcb->getJ() * qd).topRows(3);//jcb->pd[6];
	PoseE.w		= Ve.bottomRows(3);//(jcb->getJ() * qd).bottomRows(3);//jcb->w[6];
	PoseE.qt	= Quaterniond(PoseE.R);
	

	return PoseE;
}



Quaterniond Impedance::quatDerivative(const Vector3d &w, Quaterniond &qt)
{
	Quaterniond qtd;
	MatrixXd temp;

	qt.normalize();
	Vector3d v = qt.vec();

	temp = -0.5 * v.transpose()*w;

	qtd.w() = temp(0,0);
	qtd.vec() = 0.5*(qt.w()*Identity_3x3 - skew(v)) * w;/* E_fact(quat)*w;*/

	//qtd.normalize();
	return qtd;
}

Vector3d Impedance::integrate(const Vector3d &past, const Vector3d &present)
{ 
	return 0.5 * SAMPLING_T_e3 * (past + present);
}
short Impedance::integrate(Quaterniond &qt, const Quaterniond &past, const Quaterniond &present)
{ 
	qt.w()	 += /*half_sampling_T*/0.5*SAMPLING_T_e3 * (past.w()   + present.w());
	qt.vec() += /*half_sampling_T*/0.5*SAMPLING_T_e3 * (past.vec() + present.vec());

	qt.normalize();

	return 0;
}
void	Impedance::R2Quaternion(const Matrix3d R){

	
	Vector4d temp;

	temp(0) = 0.5*sqrt(R(0,0) + R(1,1) + R(2,2) + 1);
	temp(1) = 0.5*dyn.sign(R(2,1)-R(1,2))*sqrt(R(0,0)-R(1,1)-R(2,2) + 1);
	temp(2) = 0.5*dyn.sign(R(0,2)-R(2,0))*sqrt(R(1,1)-R(2,2)-R(0,0) + 1);
	temp(3) = 0.5*dyn.sign(R(1,0)-R(0,1))*sqrt(R(2,2)-R(0,0)-R(1,1) + 1);
	
	quate.w() = temp(0);

	quate.vec()<< temp(1), temp(2),temp(3);
						
}
Matrix3d Impedance::skew(const Vector3d &v) {
	Matrix3d skew;
	skew<<	      0,   -1*v(2),       v(1),
			   v(2),	     0,	   -1*v(0),
			-1*v(1),      v(0),	         0;
	return skew;  
}

//===========================================================
//Jacobian::Jacobian(int dof)
//{
//	this->dof = dof;
//	Jmat  = MatrixXd(6, dof);
//	Jmat_inv = MatrixXd(6, dof);
//	Jmatd = MatrixXd(6, dof);
//	Jmat_prev = MatrixXd::Zero(6,dof);  //MatrixXd(6, dof);
//
//	R  = new Matrix3d[dof + 1];
//	T  = new Matrix4d[dof + 1];
//	z  = new Vector3d[dof + 1];
//	w  = new Vector3d[dof + 1];
//	pd = new Vector3d[dof + 1];
//	Oi = new Vector3d[dof + 1];
//
//	R[0] = Matrix3d::Identity();
//	T[0] = Matrix4d::Identity();
//	w[0] = Vector3d::Zero();
//	z[0] = Vector3d::Zero();
//	pd[0]= Vector3d::Zero();
//}
//Jacobian::~Jacobian()
//{
//	delete[] R;
//	delete[] T;
//	delete[] z;
//	delete[] w;
//	delete[] pd;
//	delete[] Oi;
//}
//
//void Jacobian::calJ(const VectorXd &q){
//	
//	// jacobian is given as J = [Jvi Jwi].transpose
//	// Jvi = Zi-1 x (On-Oi-1), n = dof 
//	// Jwi = Zi-1 = Ri-1*k, k=[0,0,1]^T
//
//	Matrix4d tempMat;
//	
//	for(int i=1; i<= dof; i++){
//		
//		tempMat = dyn->get_A(i);
//
//		R[i]    = R[i-1]*tempMat.block<3,3>(0,0);
//		T[i]    = T[i-1]*tempMat;
//
//		z[i]  = R[i].rightCols(1);
//		Oi[i] = T[i].block<3,1>(0,3);
//	}	
//			
//	for(int i = 1;i <= dof; i++){	
//		Vector3d tempVec = z[i].cross( Oi[dof] - Oi[i] );
//		Jmat.block<3,1>(0,i-1) = tempVec;
//		Jmat.block<3,1>(3,i-1) = z[i];
//	}
//
//	for(int i = 0; i<6; i++){
//		for(int j = 0;j < dof;j++){
//			if(abs(Jmat(i,j)) < TINY_VALUE )
//				Jmat(i,j) = 0.0;
//		}
//	}
//}
//
//
//
//MatrixXd Jacobian::getJd( const VectorXd& q, const VectorXd& qd)
//{
//	Vector3d temp1,temp3;
//	Matrix4d temp;
//
//	Vector3d O_endeffectToJointi, pd_endeffectorJointiSpeed;
//	
//	/*Matrix4d T_flange;
//	T_flange<< 1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1;
//	Matrix3d R_flange;
//	R_flange<< 1,0,0,0,-1,0,0,0,-1;*/
//
//	for(int i= 1;i <=dof; i++){
//
//		temp = dyn->get_A(i);
//
//	/*	if(i==dof){
//			
//			T[i]   = T[i-1]*temp*T_flange;
//			R[i]   = R[i-1]*temp.block<3,3>(0,0)*R_flange;
//		}
//		
//		else{*/
//		T[i]   *= /*T[i-1]**/temp;
//		R[i]   *= /*R[i-1]**/temp.block<3,3>(0,0);
//		//}
//		z[i]   = R[i].rightCols(1);
//		w[i]   = w[i-1] + qd(i-1)*z[i];
//		
//		Oi[i] = T[i].block<3,1>(0,3);
//		pd[i] = pd[i-1] + R[i-1] * w[i-1].cross(Oi[i]-Oi[i-1]);
//	}	
//		
//	for(int i=1; i<=dof; i++){
//
//		//w[i] = w[i-1] + qd(i-1)*z[i];
//
//		O_endeffectToJointi = Oi[dof] - Oi[i];
//		pd_endeffectorJointiSpeed = pd[dof] - pd[i];
//				
//		temp1 = (R[i-1] * w[i-1]).cross(z[i]);
//		temp3 = temp1.cross(O_endeffectToJointi) + z[i].cross(pd_endeffectorJointiSpeed);
//		
//		Jmatd.block<3,1>(0,i-1) = temp3;
//		Jmatd.block<3,1>(3,i-1) = temp1;
//	}
//
//	for(int i=0; i<6; i++){
//		for(int j =0; j<dof; j++){
//
//			if(abs(Jmatd(i,j)) < TINY_VALUE) 
//				Jmatd(i,j) = 0.0;
//		}
//	}
//
//	return Jmatd;
//}
//
//MatrixXd Jacobian::DLSJacobInverse(VectorXd q){
//
//	//getJ(q);
//
//	double k = 0.2;
//
//	MatrixXd temp = MatrixXd::Identity(6,6);
//	temp = (Jmat*Jmat.transpose() + k*k*temp).inverse();
//
//	temp = Jmat.transpose()*temp;
//
//	return temp;
//}


