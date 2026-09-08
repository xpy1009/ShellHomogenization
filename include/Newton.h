#ifndef NEWTON_H
#define NEWTON_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class Newton
{
public:
    Newton() = delete;

    static Eigen::VectorXd descentDir(
        const Eigen::SparseMatrix<double> &hessian, 
        const Eigen::VectorXd& gradient,
        bool verbose,
        Eigen::SparseMatrix<double> *U=nullptr);

    static double lineSearch(
        std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        bool verbose);

    static void optSolve(
        std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
        Eigen::VectorXd &x0, 
        double maxIter, 
        double gradTol, 
        bool verbose=false,
        Eigen::SparseMatrix<double> *U=nullptr);

    // check gradient and hessian
    static void gradCheck(
        std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
        const Eigen::VectorXd &x0);
    static void gradConvergence(
        std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
        const Eigen::VectorXd &x0);
    static void hessCheck(
        std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
        const Eigen::VectorXd &x0);
    static void hessConvergence(
        std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
        const Eigen::VectorXd &x0);
    static constexpr double m_sEps = 1e-6;
};

#endif