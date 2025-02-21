#pragma once

#include <Eigen/Dense>

constexpr int MAX_MOTOR_PER_ROBOT = 7;

using namespace Eigen;

class DHTable {
public:
    DHTable() = default;
    ~DHTable() = default;
};

class DYNTable {
public:
    DYNTable() = default;
    ~DYNTable() = default;

    double LinkMass(int ind) const { return linkMass[ind]; }
    Vector3d LinkMassCenter(int ind) const { return linkMassCenter[ind]; }
    Matrix3d LinkI(int ind) const { return linkI[ind]; }
    double MotorInertia(int ind) const { return motorInertia[ind]; }
    double JointGearRatio(int ind) const { return jointGearRatio[ind]; }
    double JointFrictionCf(int ind) const { return jointFrictionCf[ind]; }
    double StartFriction(int ind) const { return startFriction[ind]; }
    double JointDampingB(int ind) const { return jointDampingB[ind]; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
    double linkMass[MAX_MOTOR_PER_ROBOT] = {0};
    Vector3d linkMassCenter[MAX_MOTOR_PER_ROBOT];
    Matrix3d linkI[MAX_MOTOR_PER_ROBOT];
    double motorInertia[MAX_MOTOR_PER_ROBOT] = {0};
    double jointGearRatio[MAX_MOTOR_PER_ROBOT] = {0};
    double jointFrictionCf[MAX_MOTOR_PER_ROBOT] = {0};
    double startFriction[MAX_MOTOR_PER_ROBOT] = {0};
    double jointDampingB[MAX_MOTOR_PER_ROBOT] = {0};
};
