#pragma once
//#include "Def.h"
#include "Eigen\Dense"

#define MAX_MOTOR_PER_ROBOT = 7;

using namespace Eigen;

class DHTable
{
public:
	DHTable(void);
	~DHTable(void);
};

class DYNTable
{
public:
	DYNTable(void);
	~DYNTable(void);

	double		LinkMass(int ind)		{return linkMass[ind];}
	Vector3d	LinkMassCenter(int ind)	{return linkMassCenter[ind];}
	Matrix3d	LinkI(int ind)			{return linkI[ind];}
	double		MotorInertia(int ind)	{return motorInertia[ind];}
	double		JointGearRatio(int ind)	{return jointGearRatio[ind];}
	double		JointFrictionCf(int ind){return jointFrictionCf[ind];}
	double		StartFriction(int ind)  {return startFriction[ind];}
	double		JointDampingB(int ind)	{return jointDampingB[ind];}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:

	double linkMass			[7];

	Vector3d linkMassCenter	[7];
	/*double linkMassCenterX	[MAX_MOTOR_PER_ROBOT];
	double linkMassCenterY	[MAX_MOTOR_PER_ROBOT];
	double linkMassCenterZ	[MAX_MOTOR_PER_ROBOT];*/
	Matrix3d linkI			[7];
	/*double linkIxx			[MAX_MOTOR_PER_ROBOT];
	double linkIxy			[MAX_MOTOR_PER_ROBOT];
	double linkIxz			[MAX_MOTOR_PER_ROBOT];
	double linkIyy			[MAX_MOTOR_PER_ROBOT];
	double linkIyz			[MAX_MOTOR_PER_ROBOT];
	double linkIzz			[MAX_MOTOR_PER_ROBOT];*/

	double motorInertia		[7];
	double jointGearRatio	[7];
	double jointFrictionCf	[7];
	double startFriction	[7];
	double jointDampingB	[7];

	//ashe Add it for temp us ein Jacobian check

	//double m1,m2,m3,m4,m5,m6;
	//double q1,q2,q3,q4,q5,q6;
	//double a1,a2,a3,d1,d4,d6;
	//double I1x,I2x,I3x,I4x,I5x,I6x;
	//double I1y,I2y,I3y,I4y,I5y,I6y;
	//double I1z,I2z,I3z,I4z,I5z,I6z;
	//double rc1_x,rc2_x,rc3_x,rc4_x,rc5_x,rc6_x;
	//double rc1_y,rc2_y,rc3_y,rc4_y,rc5_y,rc6_y;
	//double rc1_z,rc2_z,rc3_z,rc4_z,rc5_z,rc6_z;

	//double I1_x,I2_x,I3_x,I4_x,I5_x,I6_x;
	//double I1_y,I2_y,I3_y,I4_y,I5_y,I6_y;
	//double I1_z,I2_z,I3_z,I4_z,I5_z,I6_z;


};

