#pragma once

#include <Eigen/Dense>
#include <vector>
#include "DHTable.h"

constexpr double DEG2RAD = 0.01745329251994;
constexpr double RAD2DEG = 57.2957795130823;
constexpr double PI = 3.1415926535897932384;

using namespace Eigen;

class Dynamics {
public:
    explicit Dynamics(int dof);
    ~Dynamics();

    // Rule of Five: Copy/move constructors and assignment operators
    Dynamics(const Dynamics& other);
    Dynamics& operator=(const Dynamics& other);
    Dynamics(Dynamics&& other) noexcept;
    Dynamics& operator=(Dynamics&& other) noexcept;

    // Inverse Dynamics: get joint torque based on Newton-Euler recursive method.
    void inv(double* out_Torq, const double* q, const double* qd, const double* qdd, 
             const double F_ex[3], const double N_ext[3]);

    VectorXd inv(const VectorXd& q, const VectorXd& qd, const VectorXd& qdd,
                 const Vector3d& F_ex, const Vector3d& N_ext);
    VectorXd g(const VectorXd& q);
    VectorXd auxiliaryTorque(const VectorXd& qd);
    VectorXd frictionTorque(const VectorXd& qd);

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

    std::vector<Matrix3d> Roi; // Rotation Matrix from root to flange i.
    std::vector<Matrix3d> Rt;  // Rotation Matrix from i to i-1.
    std::vector<Matrix4d> T;   // Transformation Matrix from i-1 to i.
    std::vector<Vector3d> bi;
    std::vector<Vector3d> P;   // Position Vector from i-1 to i.
    std::vector<Vector3d> Omg; // Angular velocity Vector i.
    std::vector<Vector3d> Omg_d; // Angular acceleration Vector i.
    std::vector<Vector3d> Ae, Ac; // Acceleration Vector i at end/center of link i.
    std::vector<Vector3d> Fs, Fc; // Force Vector i at start/center of link i.
    std::vector<Vector3d> Ns, Nc; // Torque Vector i at start/center of link i.

    // Temporary variables for user-input conversion.
    VectorXd q_Vec, qd_Vec, qdd_Vec, torq_Vec;
    Vector3d F_ex_Vec, N_ex_Vec;

    std::vector<Matrix4d> matA; // Stores transformation matrix A

    // Additional parameters for gravity compensation and friction
    Vector3d gravity = {0, 0, -9.81};
    VectorXd coulomb_friction;
    VectorXd viscous_friction;
};
