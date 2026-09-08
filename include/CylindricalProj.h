#ifndef CYLINDRICAL_PROJ_H
#define CYLINDRICAL_PROJ_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class CylindricalProj
{
public:
    CylindricalProj() = delete;

    // project (x,y,z) to cylinder: ((k*z+1)/k*sin(kx), y, ((k*z+1)*cos(kx)-1)/k)
    static Eigen::MatrixX3d proj(const Eigen::MatrixX3d& w, double k, double alpha);
    static Eigen::SparseMatrix<double> gradient(const Eigen::MatrixX3d& w, double k, double alpha);
    static Eigen::SparseMatrix<double> hessian(const Eigen::MatrixX3d& w, double k, double alpha, const Eigen::VectorXd& gradient);

    static constexpr double m_sEps = 1e-10;
};

#endif