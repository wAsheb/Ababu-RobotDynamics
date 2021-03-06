//Gravity Compensation 

g(const VectorXd q){

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

		
		T[i]	= get_A(i);//get_A(q(i-1),i);
		Rt[i]	= tempRt.transpose();
		Roi[i]  *= tempRt;

		bi[i] = Roi[i].transpose()*Roi[i-1]*z0;

		//Roi[i]	= Roi[i-1] * T[i].block<3,3>(0,0);
		P[i]	= T[i].block<3,1>(0,3);		
	}

	

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
			
		}
		
		temp = (z0.transpose()*Ns[i]);
		
		out_torque(i-1) = temp(0,0);
	
	}
	
	return out_torque;

}


get_A( double theta, const int Axis_No) {
		
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
