#ifndef SHELL_H
#define SHELL_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class Shell
{    
public:
    Shell(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F);

    double energyStretch(int i, Eigen::Matrix<double, 9, 1>* grad=nullptr, Eigen::Matrix<double, 9, 9>* hess=nullptr) const;
    double energyBend(int i, Eigen::Matrix<double, 18, 1>*grad=nullptr, Eigen::Matrix<double, 18, 18>*hess=nullptr) const;

    // elasticity
    double energy() const;
    Eigen::VectorXd gradient() const;
    Eigen::SparseMatrix<double> hessian() const;

    // // homogenization
    // Eigen::SparseMatrix<double> massMatrix() const;
    // double elementVolume(int i) const;
    // Eigen::Vector3d moment() const;

    // // gravity
    // double volume() const;

    // analysis
    // std::array<Eigen::VectorXd, 2> principalCurvature() const;

    Eigen::MatrixX3d m_v;
    Eigen::MatrixX2d m_V;
    Eigen::MatrixX3i m_F;
    std::vector<std::array<int, 3>> m_flaps;
    std::vector<std::array<Eigen::Vector2d, 3>> m_ts;
    Eigen::VectorXd m_A;

    double m_h; // height
    Eigen::Matrix3d m_C; // StVK
    Eigen::Matrix3d m_B; // StVK

};

#endif