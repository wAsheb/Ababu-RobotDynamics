#include "DHTable.h"


DHTable::DHTable(void)
{
}


DHTable::~DHTable(void)
{
}


//==================================

DYNTable::DYNTable(void)
{
	//temp
	linkMass[0] = 0.0;
	linkMass[1] = 12.11;//12.11
	linkMass[2] = 9.6;//9.04;
	linkMass[3] = 2.50;
	linkMass[4] = 3.65;
	linkMass[5] = 1.25;//1.12;
	linkMass[6] = 0.25;//0.0000001;

	linkI[0](0,0) = 0.0;
	linkI[1](0,0) = 0.36158821;                   
	linkI[2](0,0) = 0.03222410;
	linkI[3](0,0) = 0.01042924;
	linkI[4](0,0) = 0.06431164;
	linkI[5](0,0) = 0.00495915;
	linkI[6](0,0) = 0.0;

	linkI[0](0,1) = linkI[0](1,0) = 0.0;
	linkI[1](0,1) = linkI[1](1,0) = -0.02137294;
	linkI[2](0,1) = linkI[2](1,0) = 0.0;
	linkI[3](0,1) = linkI[3](1,0) = 0.0;
	linkI[4](0,1) = linkI[4](1,0) = 0.0;
	linkI[5](0,1) = linkI[5](1,0) = 0.0;
	linkI[6](0,1) = linkI[6](1,0) = 0.0;

	linkI[0](0,2) = linkI[0](2,0) = 0.0;
	linkI[1](0,2) = linkI[1](2,0) = 0.13885626;                   
	linkI[2](0,2) = linkI[2](2,0) = 0.02613594;
	linkI[3](0,2) = linkI[3](2,0) = 0.0;
	linkI[4](0,2) = linkI[4](2,0) = 0.0;
	linkI[5](0,2) = linkI[5](2,0) = 0.0;
	linkI[6](0,2) = linkI[6](2,0) = 0.0;
			  
	linkI[0](1,1) = 0.0;
	linkI[1](1,1) = 0.40824519;                   
	linkI[2](1,1) = 0.34944986;
	linkI[3](1,1) = 0.0973376;
	linkI[4](1,1) = 0.05863466;
	linkI[5](1,1) = 0.00172243;
	linkI[6](1,1) = 0.0;

	linkI[0](1,2) = linkI[0](2,1)  = 0.0;
	linkI[1](1,2) = linkI[1](2,1)  = 0.0;                   
	linkI[2](1,2) = linkI[2](2,1)  = 0.0;
	linkI[3](1,2) = linkI[3](2,1)  = 0.0;
	linkI[4](1,2) = linkI[4](2,1)  = 0.0;
	linkI[5](1,2) = linkI[5](2,1)  = 0.0;
	linkI[6](1,2) = linkI[6](2,1)  = 0.0;

	linkI[0](2,2) = 0.0;
	linkI[1](2,2) = 0.12200933;                   
	linkI[2](2,2) = 0.34108476;
	linkI[3](2,2) = 0.01368245;
	linkI[4](2,2) = 0.009079792;
	linkI[5](2,2) = 0.00374815;
	linkI[6](2,2) = 0.0;

	linkMassCenter[0].x() = 0.0;
	linkMassCenter[1].x() = 0.0571;
	linkMassCenter[2].x() = 0.20;//0.1530;
	linkMassCenter[3].x() = 0.0295;
	linkMassCenter[4].x() = 0.0;
	linkMassCenter[5].x() = 0.0;
	linkMassCenter[6].x() = 0.0;

	linkMassCenter[0].y() = 0.0;
	linkMassCenter[1].y() = 0.0;
	linkMassCenter[2].y() = 0.0;
	linkMassCenter[3].y() = 0.14;//0.031;
	linkMassCenter[4].y() = 0.0;
	linkMassCenter[5].y() = 0.07;//0.045;
	linkMassCenter[6].y() = 0.0;
					 
	linkMassCenter[0].z() =  0.0367;
	linkMassCenter[1].z() =  0.0367;//-0.1448;
	linkMassCenter[2].z() =  0.0;
	linkMassCenter[3].z() =  0.0;
	linkMassCenter[4].z() =  0.065;//-0.1081;
	linkMassCenter[5].z() =  0.0;//-0.0295;
	linkMassCenter[6].z() =  0.01;

	motorInertia[0] = 0.0000001;
	motorInertia[1] = 0.00055;
	motorInertia[2] = 0.00033;
	motorInertia[3] = 0.00011;
	motorInertia[4] = 0.000019;
	motorInertia[5] = 0.0000103;
	motorInertia[6] = 0.000008;

	jointGearRatio[0] = 1.0;//donot use
	jointGearRatio[1] = 100.0;
	jointGearRatio[2] = 120.0;
	jointGearRatio[3] = 120.0;
	jointGearRatio[4] = 120.0;
	jointGearRatio[5] = 50.0;
	jointGearRatio[6] = 50.0;

	jointFrictionCf[0] = 0.0002;
	jointFrictionCf[1] = 0.00002;
	jointFrictionCf[2] = 0.0013;
	jointFrictionCf[3] = 0.00097;
	jointFrictionCf[4] = 0.0001;
	jointFrictionCf[5] = 0.0006;
	jointFrictionCf[6] = 0.0002;

//starting friction
	startFriction[0]	=	0;
	startFriction[1]	=	0;
	startFriction[2]	=	20;
	startFriction[3]	=	20;
	startFriction[4]	=	0;
	startFriction[5]	=	20;
	startFriction[6]	=	0;

	
	jointDampingB[0] = 0.000001;
	jointDampingB[1] = 0.000001;
	jointDampingB[2] = 0.000005;
	jointDampingB[3] = 0.000005;
	jointDampingB[4] = 0.000001;
	jointDampingB[5] = 0.000005;
	jointDampingB[6] = 0.000003;

	
}
DYNTable::~DYNTable(void)
{
}




