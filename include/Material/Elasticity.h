#ifndef MATERIAL_ELASTICITY_H
#define MATERIAL_ELASTICITY_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class Tet10;

class Elasticity
{
public:
    Elasticity() = delete;
    
    // Saint Venant-Kirchhoff
    static Eigen::Matrix<double, 6, 6> stiffnessTensor(double E, double nu);
    static double StVK(
        const Tet10& mesh,
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX3d& X, 
        const Eigen::Matrix<double, 6, 6>& C,
        Eigen::VectorXd* gradient=nullptr,
        Eigen::SparseMatrix<double>* hessian=nullptr);
    static Eigen::Matrix<double, 6, 1> secondPKStress(
        const Tet10& mesh, 
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX3d& X,
        const Eigen::Matrix<double, 6, 6>& C);
    static Eigen::Vector3d moment(
        const Tet10& mesh, 
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX3d& X,
        const Eigen::Matrix<double, 6, 6>& C);

    // Neo-Hookean
    static std::array<double, 2> lameParameters(double E, double nu);
    static double NeoHookean(
        const Tet10& mesh,
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX3d& X, 
        double lambda,
        double mu,
        Eigen::VectorXd* gradient=nullptr,
        Eigen::SparseMatrix<double>* hessian=nullptr);
    static Eigen::Matrix<double, 6, 1> secondPKStress(
        const Tet10& mesh, 
        const Eigen::MatrixX3d& x,
        const Eigen::MatrixX3d& X, 
        double lambda,
        double mu);
    static Eigen::Vector3d moment(
        const Tet10& mesh, 
        const Eigen::MatrixX3d& x,
        const Eigen::MatrixX3d& X, 
        double lambda,
        double mu);
};

#endif