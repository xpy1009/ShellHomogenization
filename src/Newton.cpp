#include <Newton.h>

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <numeric>
#include <iostream>

Eigen::VectorXd Newton::descentDir(
    const Eigen::SparseMatrix<double> &hessian, 
    const Eigen::VectorXd& gradient,
    bool verbose,
    Eigen::SparseMatrix<double> *U)
{
    double alpha = 1e-6;
    Eigen::SparseMatrix<double> I(hessian.rows(), hessian.cols()); I.setIdentity();
    Eigen::SparseMatrix<double> H = hessian + alpha * I;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(H);
    solver.factorize(H);

    bool spd = true;
    while(solver.info() != Eigen::Success) {
        alpha *= 2;
        H = hessian + alpha * I;
        solver.factorize(H);
        if (verbose) {
            std::cout << "H is not SPD, add alpha * I, current alpha: " << alpha << std::endl;
        }
        if (alpha > 1e4) {
            spd = false;
            if (verbose) {
                std::cout << "alpha is too large, switch to steepest descent" << std::endl;
            }
        }
    }

    const Eigen::VectorXd negGrad = -gradient;
    Eigen::VectorXd dir = spd ? solver.solve(negGrad) : negGrad;

    // solve using Woodbury matrix identity
    if (U != nullptr) {
        Eigen::MatrixXd Ainv_U(U->rows(), U->cols());
        for (int i = 0; i < U->cols(); ++i) {
            Ainv_U.col(i) = solver.solve(U->col(i));
        }
        const Eigen::Matrix3d Cinv = Eigen::Matrix3d::Identity();
        const Eigen::SparseMatrix<double> V = U->transpose();
        dir = dir - Ainv_U * (Cinv + V * Ainv_U).inverse() * V * dir;
    }

    if (spd && solver.info() != Eigen::Success) {
        throw std::runtime_error("solve failed, something is wrong");
    }
    
    return dir;
}


double Newton::lineSearch(
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dir,
    bool verbose)
{
    constexpr double c = 0.2, rho = 0.5;
    constexpr int maxCnt = 50;
    double stepSize = 1;
    int cnt = 0;

    Eigen::VectorXd grad;
    const double f = objFunc(x, &grad, nullptr);
    const double cache = c * grad.dot(dir);
    Eigen::VectorXd xNew = x + stepSize * dir;
    double fNew = objFunc(xNew, nullptr, nullptr);

    while (fNew > f + cache * stepSize) {
        stepSize *= rho;
        xNew = x + stepSize * dir;
        fNew = objFunc(xNew, nullptr, nullptr);
        if (cnt > maxCnt) {
            if (verbose) {
                std::cout << "line search max" << std::endl;
            }
            break;
        }
        cnt++;
    }
    return stepSize;
}

void Newton::optSolve(
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
    Eigen::VectorXd &x0, 
    double maxIter, 
    double gradTol, 
    bool verbose,
    Eigen::SparseMatrix<double> *U)
{
    Eigen::VectorXd grad;
    Eigen::SparseMatrix<double> hessian;

    for (int iter = 0; iter < maxIter; ++iter) {
        const double f = objFunc(x0, &grad, &hessian);
        if (verbose) {
            std::cout << "[NEWTON] iter " << iter << "/" << maxIter << ", f: " << f <<  ", grad norm: " << grad.norm() << std::endl;
        }

        if (grad.norm() < gradTol) {
            return;
        }

        const Eigen::VectorXd dir = Newton::descentDir(hessian, grad, verbose, U);

        const double stepSize = lineSearch(objFunc, x0, dir, verbose);
        x0 += stepSize * dir;
    }
}


void Newton::gradCheck(
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
    const Eigen::VectorXd &x0)
{
    Eigen::VectorXd grad;
    objFunc(x0, &grad, nullptr);
    Eigen::VectorXd x = x0;

    for (int i = 0; i < grad.size(); ++i) {
        x(i) = x0(i) + m_sEps;
        const double e0 = objFunc(x, nullptr, nullptr);

        x(i) = x0(i) - m_sEps;
        const double e1 = objFunc(x, nullptr, nullptr);

        x(i) = x0(i);
        const double fdi = (e0 - e1) / (2.0 * m_sEps);

        if ((abs(fdi - grad(i)) < 1e-4) || abs(fdi - grad(i)) < abs(1e-3 * grad(i))) {
            continue;
        }
        std::cout << i << " fd: " << fdi << "\t analytic: " << grad(i) << " diff: " << abs(fdi-grad(i)) << std::endl;

    }
    std::cout << "gradient checked, size: " << grad.size() << ", norm: " << grad.norm() << std::endl;
}

void Newton::gradConvergence(
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
    const Eigen::VectorXd &x0)
{
    std::cout << "gradConvergence: " << std::endl;
    Eigen::VectorXd grad;
    objFunc(x0, &grad, nullptr);

    const Eigen::VectorXd dir = Eigen::VectorXd::Random(grad.size());
    const double analytic = grad.dot(dir);

    // 10^-5 ~ 10^-10
    std::vector<double> eps(6);
    std::iota(eps.begin(), eps.end(), 5);
    for (size_t i = 0; i < eps.size(); ++i) {
        const double ti = std::pow(10.0, -eps[i]);
        Eigen::VectorXd x = x0 + ti * dir;
        const double ep = objFunc(x, nullptr, nullptr);
        x = x0 - ti * dir;
        const double en = objFunc(x, nullptr, nullptr);
        const double fdi = (ep - en) / (2.0 * ti);
        std::cout << ti << " " << abs(analytic - fdi) / abs(analytic) << std::endl;
    }
}


void Newton::hessCheck(
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
    const Eigen::VectorXd &x0)
{
    Eigen::SparseMatrix<double> hess;
    Eigen::VectorXd grad;
    objFunc(x0, &grad, &hess);
    Eigen::VectorXd x = x0;

    for (int j = 0; j < hess.cols(); ++j) {
        x(j) = x0(j) + m_sEps;
        Eigen::VectorXd g0;
        objFunc(x, &g0, nullptr);

        x(j) = x0(j) - m_sEps;
        Eigen::VectorXd g1;
        objFunc(x, &g1, nullptr);

        x(j) = x0(j);
        const Eigen::VectorXd dgdxj = (g0 - g1) / (2 * m_sEps);

        for (int i = 0; i < hess.rows(); ++i) {
            if ((abs(hess.coeff(i, j)) < 1e-6 && abs(dgdxj(i)) < 1e-6) || abs(hess.coeff(i, j) - dgdxj(i)) < abs(1e-3 * dgdxj(i)))
                continue;
            std::cout << i << ", " << j  << " fd: " << dgdxj(i) << "\t analytic: " << hess.coeff(i, j) << " diff: " << abs(dgdxj(i)-hess.coeff(i, j)) << std::endl;            
        }
        
    }
    std::cout << "checkHessian done, size: " << hess.rows() << "x" << hess.cols() << ", norm: " << hess.norm() << std::endl;
}


void Newton::hessConvergence(
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> objFunc,
    const Eigen::VectorXd &x0)
{
    std::cout << "hessConvergence: " << std::endl;
    Eigen::SparseMatrix<double> hess;
    Eigen::VectorXd grad;
    objFunc(x0, &grad, &hess);

    const Eigen::VectorXd dir = Eigen::VectorXd::Random(hess.rows());
    const Eigen::VectorXd analytic = hess * dir;

    // 10^-3 ~ 10^-9
    std::vector<double> eps(7);
    std::iota(eps.begin(), eps.end(), 3);
    for (size_t i = 0; i < eps.size(); ++i) {
        const double ti = std::pow(10.0, -eps[i]);
        Eigen::VectorXd x = x0 + ti * dir;
        Eigen::VectorXd gp, gn;
        objFunc(x, &gp, nullptr);

        x = x0 - ti * dir;
        objFunc(x, &gn, nullptr);
        const Eigen::VectorXd fdi = (gp - gn) / (2.0 * ti);

        std::cout << ti << " " << (analytic - fdi).norm() / analytic.norm() << std::endl;
    }
}
