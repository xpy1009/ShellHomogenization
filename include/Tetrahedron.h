#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class Tetrahedron
{    
public:
    Tetrahedron(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T);

    // basis functions
    Eigen::Matrix<double, 10, 1> basisFunction(const Eigen::Vector3d &q) const;
    Eigen::Matrix<double, 10, 3> basisFunctionGradient(const Eigen::Vector3d &q) const;

    // constitutive model
    enum class Model { StVK, NeoHookean };
    Model m_model = Model::StVK;
    // stiffness matrix
    Eigen::Matrix<double, 6, 6> C;
    // Lame parameters
    double m_mu, m_lambda;

    // elasticity
    double energy() const;
    Eigen::VectorXd gradient() const;
    Eigen::SparseMatrix<double> hessian() const;

    // homogenization
    Eigen::SparseMatrix<double> massMatrix() const;
    double elementVolume(int i) const;
    Eigen::Vector3d moment() const;

    // gravity
    double volume() const;

    Eigen::MatrixXd undeformed, deformed;
    Eigen::MatrixXi faces, tets;

    // weights and quadrature points
    std::vector<std::pair<double,Eigen::Vector3d>> qs;

};

#endif