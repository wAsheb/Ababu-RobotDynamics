
#include "Dynamics.h"
#include "Eigen\Dense"
//#include "Log.h"

using namespace Eigen;

#define TINY_VALUE 1e-10;

//Mostly used from d P(q)/dq assuming TCP position is function of joint angle
class Jacobian
{
public:
	Jacobian( int dof);
	~Jacobian();

	void calJ ( const VectorXd& q);
	MatrixXd getJ (){return Jmat;}
	MatrixXd getJd( const VectorXd& q, const VectorXd& qd); // Derivative of geometric jacobian w.r.t time
	MatrixXd DLSJacobInverse(VectorXd q);  //Singular value decomposition for Jacobian Inverse
	Vector3d *z, *w, *pd, *Oi;

	MatrixXd Jmat,Jmat_inv;  // Jacobian Matrix.

	Dynamics* dyn; //only to let Jacobion get "get_A(i)" result.

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
	int dof;
	//MatrixXd Jmat,Jmat_inv;  // Jacobian Matrix.
	MatrixXd Jmatd; // derivative of Jmat
	MatrixXd Jmat_prev;

	Matrix3d *R;
	Matrix4d *T;
	//Vector3d *z, *w, *pd, *Oi;

	Matrix4d T_flange;
	Matrix3d R_flange;
	Matrix4d tempMat;
	Vector3d tempVec;
};
