#include <Compression.h>
#include <Element/Tet10.h>
#include <Element/Shell.h>
#include <Material/Elasticity.h>
#include <Material/Neural.h>
#include <Newton.h>

#include <fstream>
#include <Eigen/Dense>


Eigen::MatrixX3d Compression::run(
    const Tet10& mesh, 
    const Eigen::MatrixX3d& restPos, 
    double E,
    double nu,
    std::ofstream &file,
    bool verbose)
{
    auto [lambda, mu] = Elasticity::lameParameters(E, nu);

    // move center to origin
    Eigen::MatrixX3d curPos = restPos;
    const Eigen::RowVector3d center = 0.5 * (curPos.colwise().maxCoeff() - curPos.colwise().minCoeff());
    curPos.rowwise() -= center;

    // boundary condition
    std::vector<int> fixedVar;
    const double xMin = curPos.col(0).minCoeff(), xMax = curPos.col(0).maxCoeff();
    for (int i = 0; i < curPos.rows(); ++i) {
        if ((abs(curPos(i,0)-xMin) < m_sEps || abs(curPos(i,0)-xMax) < m_sEps) && abs(curPos(i,2)) < m_sEps) {
            fixedVar.push_back(3*i);
            fixedVar.push_back(3*i+1);
            fixedVar.push_back(3*i+2);
        }
    }

    // mid-surface
    std::vector<int> midIdx;
    for (int i = 0; i < curPos.rows(); ++i) {
        if (abs(curPos(i,2)) < m_sEps) {
            midIdx.push_back(i);
        }
    }

    
    auto objFunc = [&](const Eigen::VectorXd &var, Eigen::VectorXd *grad, Eigen::SparseMatrix<double> *hessian) {
        const Eigen::MatrixX3d pos = var.reshaped<Eigen::RowMajor>(var.size()/3, 3);

        const double energy = Elasticity::NeoHookean(mesh, pos, restPos, lambda, mu, grad, hessian);

        if (grad != nullptr) {
            (*grad)(fixedVar).setZero();
        }

        if (hessian != nullptr) {
            for (const int& idx : fixedVar) {
                hessian->row(idx) *= 0.0;
                hessian->col(idx) *= 0.0;
                hessian->coeffRef(idx,idx) = 1.0;
            }
        }

        return energy;
    };

    Eigen::VectorXd var = curPos.reshaped<Eigen::RowMajor>();

    const double origL = restPos.col(0).maxCoeff() - restPos.col(0).minCoeff();
    double prevL = origL;
    for (size_t i = 0; i < 50; ++i) {
        const double L = prevL - 1e-3;
        // initial guess
        curPos.col(0) *= L / prevL;
        curPos.col(2) += 2e-3 * (curPos.col(0) / L * M_PI).array().cos().matrix();

        Eigen::VectorXd var = curPos.reshaped<Eigen::RowMajor>();
        Newton::optSolve(objFunc, var, m_sMaxIter, m_sTol, verbose);
        curPos = var.reshaped<Eigen::RowMajor>(var.size()/3, 3);
        prevL = L;


        // fit
        const Eigen::MatrixX3d midPos = curPos(midIdx, Eigen::placeholders::all);
        auto [k1, k2] = fit(midPos);

        const double dl = origL - L;
        file << dl << " " << k1 << " " << k2 << std::endl;
    }

    return curPos;
}


Eigen::MatrixX3d Compression::run(
    const Shell& mesh, 
    const Eigen::MatrixX2d& restPos, 
    const torch::jit::Module& stretchModel,
    const torch::jit::Module& bendModel,
    std::ofstream &file,
    bool verbose)
{
    // boundary condition
    std::vector<int> fixedVar;
    const double xMin = restPos.col(0).minCoeff(), xMax = restPos.col(0).maxCoeff();
    for (int i = 0; i < restPos.rows(); ++i) {
        if (abs(restPos(i,0)-xMin) < m_sEps || abs(restPos(i,0)-xMax) < m_sEps) {
            fixedVar.push_back(3*i);
            fixedVar.push_back(3*i+1);
            fixedVar.push_back(3*i+2);
        }
    }


    const double origL = restPos.col(0).maxCoeff() - restPos.col(0).minCoeff();
    double prevL = origL;
    Eigen::MatrixX3d curPos(restPos.rows(), 3);
    curPos.leftCols(2) = restPos;
    curPos.col(2).setZero();
    // move center to origin
    const Eigen::RowVector3d center = curPos.colwise().mean();
    curPos.rowwise() -= center;

    
    auto objFunc = [&](const Eigen::VectorXd &var, Eigen::VectorXd *grad, Eigen::SparseMatrix<double> *hessian) {
        const Eigen::MatrixX3d pos = var.reshaped<Eigen::RowMajor>(var.size()/3, 3);

        const double energy = Neural::elasticEnergy(mesh, pos, restPos, stretchModel, bendModel, grad, hessian);

        if (grad != nullptr) {
            (*grad)(fixedVar).setZero();
        }

        if (hessian != nullptr) {
            for (const int& idx : fixedVar) {
                hessian->row(idx) *= 0.0;
                hessian->col(idx) *= 0.0;
                hessian->coeffRef(idx,idx) = 1.0;
            }
        }

        return energy;
    };


    for (size_t i = 0; i < 50; ++i) {
        const double L = prevL - 1e-3;
        // initial guess
        curPos.col(0) *= L / prevL;
        curPos.col(2) += 2e-3 * (curPos.col(0) / L * M_PI).array().cos().matrix();


        Eigen::VectorXd var = curPos.reshaped<Eigen::RowMajor>();

        Newton::optSolve(objFunc, var, m_sMaxIter, m_sTol, verbose);
        curPos = var.reshaped<Eigen::RowMajor>(var.size()/3, 3);
        prevL = L;


        // fit
        // const Eigen::MatrixX3d midPos = curPos;
        // const double C = M_PI / L;
        // const Eigen::VectorXd sincy = (C * midPos.col(0)).array().cos();
        // const Eigen::VectorXd w = midPos.col(2);
        // Eigen::MatrixX3d M(midPos.rows(), midPos.cols());
        // M.col(0) = sincy;
        // M.col(1) = midPos.col(1).array().square() * sincy.array();
        // M.col(2).setOnes();
        // const Eigen::Vector3d x = (M.transpose()*M).inverse() * (M.transpose()*w);
        // const double A = x(0), B = x(1);
        // const double k1 = -A * C * C;
        // const double k2 = 2 * B;
        // std::cout << "k1: " << k1 << " / k2: " << k2 << " D: " << x(2) << std::endl;

        auto [k1, k2] = fit(curPos);

        const double dl = origL - L;
        file << dl << " " << k1 << " " << k2 << std::endl;
    }

    return curPos;
}

std::array<double, 2> Compression::fit(const Eigen::MatrixX3d& pos)
{
    const double L = pos.col(0).maxCoeff() - pos.col(0).minCoeff();
    const double C = M_PI / L;
    const Eigen::VectorXd sincy = (C * pos.col(0)).array().cos();
    const Eigen::VectorXd w = pos.col(2);

    Eigen::MatrixX3d M(pos.rows(), pos.cols());
    M.col(0) = sincy;
    M.col(1) = pos.col(1).array().square() * sincy.array();
    M.col(2).setOnes();
    const Eigen::Vector3d x = (M.transpose() * M).inverse() * (M.transpose() * w);

    const double A = x(0), B = x(1);
    const double k1 = -A * C * C;
    const double k2 = 2 * B;
    return {k1, k2};
}