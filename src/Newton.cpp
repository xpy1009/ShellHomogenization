#include <Newton.h>

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <iostream>

Eigen::VectorXd Newton::solve(Eigen::SparseMatrix<double> &A, const Eigen::VectorXd& b)
{
    double alpha = 1e-6;
    Eigen::VectorXd x;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    for (int i = 0; i < 50; ++i)
    {
        solver.factorize(A);
        if (solver.info() == Eigen::NumericalIssue) {
            A.diagonal().array() += alpha;
            alpha *= 10;
            std::cout << "factorize fails, add diag" << std::endl;
            continue;
        }

        x = solver.solve(b);

        if(((A*x - b).norm() / b.norm() < 1e-6) && (x.dot(b) > 0.0)) {
            break;
        }

        std::cout << "linear solve fails, add diag" << std::endl;
        A.diagonal().array() += alpha;
        alpha *= 10;
   }
    
    return x;
}

Eigen::VectorXd Newton::solve(Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &U, const Eigen::VectorXd& b)
{
    double alpha = 1e-6;
    Eigen::VectorXd x;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    for (int i = 0; i < 50; ++i)
    {
        solver.factorize(A);
        if (solver.info() == Eigen::NumericalIssue) {
            A.diagonal().array() += alpha;
            alpha *= 10;
            std::cout << "factorize fails, add diag" << std::endl;
            continue;
        }

        const Eigen::VectorXd Ainv_b = solver.solve(b);
        Eigen::MatrixXd Ainv_U(U.rows(), U.cols());
        for (int i = 0; i < U.cols(); ++i) {
            Ainv_U.col(i) = solver.solve(U.col(i));
        }
        const Eigen::Matrix3d Cinv = Eigen::Matrix3d::Identity();
        const Eigen::SparseMatrix<double> V = U.transpose();
        x = Ainv_b - Ainv_U * (Cinv + V * Ainv_U).inverse() * V * Ainv_b;

        if(((A*x+U*(V*x) - b).norm() / b.norm() < 1e-6)
            && ((A*Ainv_b - b).norm() / b.norm() < 1e-6)
            && ((A*Ainv_U - U).norm() / U.norm() < 1e-6)) {
            break;
        }

        std::cout << "linear solve fails, add diag" << std::endl;
        A.diagonal().array() += alpha;
        alpha *= 10;
   }
    
    return x;
}

void Newton::solveOpt(Objective& obj, double maxIter, double tol)
{
    for (int iter = 0; iter < maxIter; ++iter) {
        Eigen::VectorXd grad = obj.gradient();
        
        double grad_norm = grad.norm();
        std::cout << "[NEWTON] iter " << iter << "/" << maxIter << " gradient norm: " << grad_norm << " tol: " << tol << std::endl;

        if (grad_norm < tol || iter == maxIter) {
            return;
        }

        // linear solve
        Eigen::SparseMatrix<double> K = obj.hessian();
        const Eigen::VectorXd dir = Newton::solve(K, -grad);

        // line search
        constexpr int maxCnt = 15;
        const double E0 = obj.energy();
        double stepSize = 1;
        int cnt = 0;

        while (true) {
            Eigen::VectorXd dx = stepSize * dir;
            obj.update(dx);

            const double E1 = obj.energy();

            if (E1 - E0 < 0 || cnt > maxCnt) {
                if (cnt > maxCnt) {
                    std::cout << "line search max" << std::endl;
                }
                return;
            }

            dx *= -1;
            obj.update(dx);
            stepSize *= 0.5;
            cnt++;
        }        
    }
}