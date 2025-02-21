#pragma once

#include <Eigen/Dense>
#include "DHTable.h"

constexpr double DEG2RAD = 0.01745329251994;
constexpr double RAD2DEG = 57.2957795130823;
constexpr double PI = 3.1415926535897932384;

using namespace Eigen;

class Dynamics {
public:
    Dynamics(int dof);
    ~Dynamics();

    // Inverse Dynamics: get joint torque based on Newton-Euler recursive method.
    void inv(double* out_Torq, const double* q, const double* qd, const double* qdd, 
             const double F_ex[3], const double N_ext[3]);

    VectorXd inv(const VectorXd& q, const VectorXd& qd, const VectorXd& qdd,
                 const Vector3d& F_ex, const Vector3d& N_ext);
    VectorXd g(const VectorXd& q);
    VectorXd auxiliaryTorque(const VectorXd& qd);

    Matrix4d get_A(double theta, int Axis_No);
    Matrix4d get_A(int Axis_No);
    void cal_allA(const VectorXd& q);

    int sign(double value);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
    DYNTable DYN;

    int dof;
    Matrix3d Identity_3x3;
    Vector3d Zero_3;
    Vector3d z0; // Z unit vector [0,0,1]^T

    Matrix3d* Roi; // Rotation Matrix from root to flange i.
    Matrix3d* Rt;  // Rotation Matrix from i to i-1.
    Matrix4d* T;   // Transformation Matrix from i-1 to i.
    Vector3d* bi;
    Vector3d* P;   // Position Vector from i-1 to i.
    Vector3d* Omg; // Angular velocity Vector i.
    Vector3d* Omg_d; // Angular acceleration Vector i.
    Vector3d* Ae, *Ac; // Acceleration Vector i at end/center of link i.
    Vector3d* Fs, *Fc; // Force Vector i at start/center of link i.
    Vector3d* Ns, *Nc; // Torque Vector i at start/center of link i.

    // Temporary variables for user-input conversion.
    VectorXd q_Vec, qd_Vec, qdd_Vec, torq_Vec;
    Vector3d F_ex_Vec, N_ex_Vec;

    Matrix4d* matA; // Stores transformation matrix A
};
