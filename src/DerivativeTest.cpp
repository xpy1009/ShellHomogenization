#include <DerivativeTest.h>
#include <Objective.h>

#include <iostream>

void DerivativeTest::checkGradient(Objective &obj)
{
    constexpr double eps = 1e-6;
    Eigen::VectorXd grad = obj.gradient();
    Eigen::VectorXd dx(grad.size()); dx.setZero();
    Eigen::VectorXd grad_fd(grad.size());

    for (int i = 0; i < grad.size(); ++i)
    {
        dx(i) = eps;
        obj.update(dx);
        double e0 = obj.energy();

        dx(i) = -2.0 * eps;
        obj.update(dx);
        double e1 = obj.energy();

        grad_fd(i) = (e0 - e1) / (2.0 * eps);
        dx(i) = eps;
        obj.update(dx);

        dx(i) = 0;

        if ((abs(grad_fd(i) - grad(i)) < 1e-4) || abs(grad_fd(i) - grad(i)) < abs(1e-3 * grad(i))) {
            continue;
        }
        std::cout << i << " fd: " << grad_fd(i) << "\t analytic: " << grad(i) << " diff: " << abs(grad_fd(i)-grad(i)) << std::endl;

    }
    std::cout << "checkGradient done, size: " << grad.size() << ", norm: " << grad.norm() << std::endl;
}


void DerivativeTest::checkHessian(Objective &obj)
{
    constexpr double eps = 1e-6;
    const Eigen::SparseMatrix<double> hess = obj.hessian();
    Eigen::VectorXd dx(hess.rows()); dx.setZero();

    for (int j = 0; j < hess.cols(); j++)
    {
        dx(j) = eps;
        obj.update(dx);
        Eigen::VectorXd g0 = obj.gradient();

        dx(j) = -2.0 * eps;
        obj.update(dx);
        Eigen::VectorXd g1 = obj.gradient();

        dx(j) = eps;
        obj.update(dx);

        dx(j) = 0;

        const Eigen::VectorXd dgdxj = (g0 - g1) / (2 * eps);

        for (int i = 0; i < hess.rows(); ++i)
        {
            if ((abs(hess.coeff(i, j)) < 1e-6 && abs(dgdxj(i)) < 1e-6) || abs(hess.coeff(i, j) - dgdxj(i)) < abs(1e-3 * dgdxj(i)))
                continue;
            std::cout << i << ", " << j  << " fd: " << dgdxj(i) << "\t analytic: " << hess.coeff(i, j) << " diff: " << abs(dgdxj(i)-hess.coeff(i, j)) << std::endl;            
        }
        
    }
    std::cout << "checkHessian done, size: " << hess.rows() << "x" << hess.cols() << ", norm: " << hess.norm() << std::endl;
}
