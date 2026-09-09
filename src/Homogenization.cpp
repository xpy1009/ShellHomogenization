#include <Homogenization.h>
#include <Material/Elasticity.h>
#include <Newton.h>
#include <CylindricalProj.h>

#include <fstream>
#include <iostream>

Homogenization::Homogenization(
    const Tet10& mesh, 
    const Eigen::SparseMatrix<double>& proj, 
    const Eigen::VectorXd& var, 
    const std::array<double, 4>& trans,
    double E,
    double nu) : m_mesh(mesh), m_proj(proj), m_restVar(var)
{
    auto [lambda, mu] = Elasticity::lameParameters(E, nu);
    m_lambda = lambda, m_mu = mu;

    
    m_area = std::abs(trans[0] * trans[3] - trans[1] * trans[2]);
    m_restPos = (proj * var).reshaped<Eigen::RowMajor>(proj.rows()/3, 3);
    m_massMat = m_mesh.massMatrix(m_restPos);
    m_U = WoodburyMatrix(m_massMat);
    setDir(0.0);
    setWeight(1e4);


    m_objFunc = [&](const Eigen::VectorXd &var, Eigen::VectorXd *grad, Eigen::SparseMatrix<double> *hessian) {
        const Eigen::MatrixX3d w = (m_rProj * var).reshaped<Eigen::RowMajor>(m_rProj.rows()/3, 3);
        const Eigen::MatrixX3d pos = CylindricalProj::proj(w, m_curvature, m_angle);

        double energy = Elasticity::NeoHookean(m_mesh, pos, m_restPos, m_lambda, m_mu, grad, hessian);
        const Eigen::RowVector3d cs = (m_massMat * w).colwise().sum();
        energy += 0.5 * m_weight * cs.squaredNorm();

        if (grad) {
            // cylindrical projection
            const Eigen::SparseMatrix<double> dxdw = CylindricalProj::gradient(w, m_curvature, m_angle);
            *grad = dxdw.transpose() * (*grad);
            // constraints
            Eigen::MatrixX3d e(m_restPos.rows(), m_restPos.cols()); e.rowwise() = cs;
            *grad += m_weight * (m_massMat * e).reshaped<Eigen::RowMajor>();
            // project
            *grad = m_rProj.transpose() * (*grad);

            (*grad)(m_fixedVar).setZero();
        }

        if (hessian) {
            Eigen::VectorXd gradient;
            Elasticity::NeoHookean(m_mesh, pos, m_restPos, m_lambda, m_mu, &gradient);
            // cylindrical projection
            const Eigen::SparseMatrix<double> dxdw = CylindricalProj::gradient(w, m_curvature, m_angle);
            const Eigen::SparseMatrix<double> d2xdw2 = CylindricalProj::hessian(w, m_curvature, m_angle, gradient);
            *hessian = dxdw.transpose() * (*hessian) * dxdw + d2xdw2;
            // project
            *hessian = m_rProj.transpose() * (*hessian) * m_rProj;

            for (const int& idx : m_fixedVar) {
                hessian->row(idx) *= 0.0;
                hessian->col(idx) *= 0.0;
                hessian->coeffRef(idx,idx) = 1.0;
            }
        }

        return energy;
    };

    // move mesh center to origin, align midsurface to z=0
    Newton::optSolve(m_objFunc, m_restVar, m_sMaxIter, m_sTol, false, &m_Uproj);
    m_restPos = (m_rProj * m_restVar).reshaped<Eigen::RowMajor>(m_rProj.rows()/3, 3);
}

void Homogenization::setDir(double angle)
{
    m_rProj = m_proj;
    const int j = m_proj.cols() - 3;
    for (Eigen::SparseMatrix<double>::InnerIterator it(m_rProj, j); it; ++it) {
        const int i = it.row();
        const double dx = m_rProj.coeff(i, j), dy = m_rProj.coeff(i, j+1);

        m_rProj.coeffRef(i, j) = cos(angle) * cos(angle) * dx + cos(angle) * sin(angle) * dy;
        m_rProj.coeffRef(i, j+1) = (cos(angle) * cos(angle) - sin(angle) * sin(angle)) * dy - 2 * cos(angle) * sin(angle) * dx;
        m_rProj.coeffRef(i, j+2) = sin(angle) * sin(angle) * dx - cos(angle) * sin(angle) * dy;
        m_rProj.coeffRef(i+1, j) = cos(angle) * sin(angle) * dx + sin(angle) * sin(angle) * dy;
        m_rProj.coeffRef(i+1, j+1) = (cos(angle) * cos(angle) - sin(angle) * sin(angle)) * dx + 2 * cos(angle) * sin(angle) * dy;
        m_rProj.coeffRef(i+1, j+2) = cos(angle) * cos(angle) * dy - cos(angle) * sin(angle) * dx;
        ++it;
    }
    // Woodbury matrix needs to be reset
    setWeight(m_weight);
}

void Homogenization::setWeight(double weight)
{
    m_weight = weight;
    m_Uproj = m_rProj.transpose() * m_U;
    m_Uproj = sqrt(m_weight) * m_Uproj;
    for (const auto& i : m_fixedVar) {
        m_Uproj.row(i) *= 0.0;
    }
}

void Homogenization::solve(Eigen::VectorXd& var, bool verbose)
{
    Newton::optSolve(m_objFunc, var, m_sMaxIter, m_sTol, verbose, &m_Uproj);
    Eigen::MatrixX3d w = (m_rProj * var).reshaped<Eigen::RowMajor>(m_rProj.rows()/3, 3);
    Eigen::RowVector3d cs = (m_massMat * w).colwise().sum();
    while (cs.maxCoeff() > m_sConstraintTol && m_weight < m_sMaxWeight) {
        if (verbose) {
            std::cout << "Constraints: " << cs << " not statisfied, increase weight" << std::endl;
        }
        setWeight(m_weight * 1e2);
        Newton::optSolve(m_objFunc, var, m_sMaxIter, m_sTol, verbose, &m_Uproj);
        w = (m_rProj * var).reshaped<Eigen::RowMajor>(m_rProj.rows()/3, 3);
        cs = (m_massMat * w).colwise().sum();
    }
    if (cs.maxCoeff() > m_sConstraintTol) {
        throw std::runtime_error("constraints not satisfied");
    }
}

void Homogenization::stretch(std::ofstream &file, bool verbose)
{
    file << "# s1 s2 angle E00 E11 E01 S00 S11 S01 Psi" << std::endl;

    Eigen::VectorXd var = m_restVar;
    const int s1Idx = m_restVar.size() - 3;

    auto sim = [&](double s1, double s2, double angle){
        // initial guess
        var = m_restVar;
        Eigen::Matrix2d T;
        T << 1.0 + s1, 0.0, 0.0, 1.0 + s2;
        Eigen::Matrix2d R; R << cos(angle), -sin(angle), sin(angle), cos(angle);
        T = R * T * R.transpose();
        for (int i = 0; i < var.size() - 3; i += 3) {
            var.segment(i, 2) = T * var.segment(i, 2);
        }
        var.tail(3) << 1.0 + s1, 0, 1.0 + s2;

        solve(var, verbose);
        
        Eigen::VectorXd grad;
        m_objFunc(var, &grad, nullptr);
        if (grad.norm() < m_sTol) {
            saveStretchData(file, var, angle);
        } else {
            const std::vector<int> tmp = m_fixedVar;
            m_fixedVar = {s1Idx, s1Idx+1, s1Idx+2};
            Newton::optSolve(m_objFunc, var, m_sMaxIter, m_sTol, verbose, &m_Uproj);
            m_objFunc(var, &grad, nullptr);
            m_fixedVar = tmp;
            if (grad.norm() < m_sTol) {
                saveStretchData(file, var, angle);
            }
        }

        // reset weight
        setWeight(m_sInitWeight);
    };

    
    std::cout << "Uniaxial stretching..." << std::endl;
    m_fixedVar = { s1Idx };
    constexpr double ds = 0.01;
    int nAngles = 20, nStretch = 16, nCompress = 8;

    for (int i = 0; i < nAngles; ++i) {
        const double angle = i * M_PI / nAngles;
        setDir(angle);

        for (int j = 1; j < nStretch + 1; ++j) {
            const double sj = ds * j;
            sim(sj, 0, angle);
        }
        
        for (int j = 1; j < nCompress + 1; ++j) {
            const double sj = -ds * j;
            sim(sj, 0, angle);
        }
    }

    
    std::cout << "Biaxial stretching..." << std::endl;
    m_fixedVar = { s1Idx, s1Idx+2 };
    nAngles = 10, nStretch = 8, nCompress = 4;
    for (int i = 0; i < nAngles; ++i) {
        const double angle = i * M_PI_2 / nAngles;
        setDir(angle);

        for (int j = 1; j < nStretch + 1; ++j) {
            const double sj = ds * j;
            for (int k = 1; k < nStretch + 1; ++k) {
                const double sk = ds * k;
                sim(sj, sk, angle);
            }
        }
        
        for (int j = 1; j < nStretch + 1; ++j) {
            const double sj = ds * j;
            for (int k = 1; k < nCompress + 1; ++k) {
                const double sk = -ds * k;
                sim(sj, sk, angle);
            }
        }

        for (int j = 1; j < nCompress + 1; ++j) {
            const double sj = -ds * j;
            for (int k = 1; k < nStretch + 1; ++k) {
                const double sk = ds * k;
                sim(sj, sk, angle);
            }
        }

        for (int j = 1; j < nCompress + 1; ++j) {
            const double sj = -ds * j;
            for (int k = 1; k < nCompress + 1; ++k) {
                const double sk = -ds * k;
                sim(sj, sk, angle);
            }
        }
    }
}

void Homogenization::saveStretchData(std::ofstream &file, const Eigen::VectorXd& var, double angle) const
{
    // Green strain
    const Eigen::Vector3d pbc = var.tail(3);
    Eigen::Matrix2d F; F << pbc(0), pbc(1), pbc(1), pbc(2);
    Eigen::Matrix2d R; R << cos(angle), -sin(angle), sin(angle), cos(angle);
    F = R * F * R.transpose();
    const Eigen::Matrix2d E = 0.5 * (F.transpose() * F - Eigen::Matrix2d::Identity());
    Eigen::Vector3d EVoigt; EVoigt << E(0, 0), E(1, 1), 2.0 * E(0, 1);
    
    // Second Piola Kirchhoff stress
    const Eigen::MatrixX3d pos = (m_rProj * var).reshaped<Eigen::RowMajor>(m_rProj.rows()/3, 3);
    const Eigen::Matrix<double, 6, 1> SVoigt = Elasticity::secondPKStress(m_mesh, pos, m_restPos, m_lambda, m_mu) / m_area;

    // Elastic energy
    const double psi = Elasticity::NeoHookean(m_mesh, pos, m_restPos, m_lambda, m_mu) / m_area;

    const double s1 = pbc(0), s2 = m_fixedVar.size()==1 ? 0.0 : pbc(2);
    file << s1 << " " << s2 << " "  << angle << " "
        << EVoigt(0) << " " << EVoigt(1) << " " << EVoigt(2) << " " 
        << SVoigt(0) << " " << SVoigt(1) << " " << SVoigt(5) << " "
        << psi << std::endl;
}

void Homogenization::bend(std::ofstream &file, bool verbose)
{
    file << "# k angle k00 k11 k01 m00 m11 m01 Psi" << std::endl;

    Eigen::VectorXd var = m_restVar;
    Eigen::VectorXd grad;

    constexpr double dk = 10;
    constexpr int nAngles = 25, nBend = 15;

    std::cout << "Uniaxial bending..." << std::endl;
    for (int i = 0; i < nAngles; ++i) {
        m_angle = i * M_PI / nAngles;

        var = m_restVar;
        setWeight(m_sInitWeight);
        for (int j = 1; j < nBend + 1; ++j) {
            m_curvature = dk * j;
            solve(var, verbose);
            m_objFunc(var, &grad, nullptr);
            if (grad.norm() < m_sTol) {
                saveBendData(file, var);
            }
        }
        
        var = m_restVar;
        setWeight(m_sInitWeight);
        for (int j = 1; j < nBend + 1; ++j) {
            m_curvature = -dk * j;
            solve(var, verbose);
            m_objFunc(var, &grad, nullptr);
            if (grad.norm() < m_sTol) {
                saveBendData(file, var);
            }
        }
    }
}

void Homogenization::saveBendData(std::ofstream &file, const Eigen::VectorXd& var) const
{
    Eigen::Matrix2d R; R << cos(m_angle), -sin(m_angle), sin(m_angle), cos(m_angle);
    Eigen::Matrix2d S; S << m_curvature, 0, 0, 0;
    S = R * S * R.transpose();
    const Eigen::Vector3d kVoigt(S(0,0), S(1,1), 2.0 * S(0,1));
    
    const Eigen::MatrixX3d w = (m_rProj * var).reshaped<Eigen::RowMajor>(m_rProj.rows()/3, 3);
    const Eigen::MatrixX3d pos = CylindricalProj::proj(w, m_curvature, m_angle);
    const Eigen::Vector3d mVoigt = Elasticity::moment(m_mesh, pos, m_restPos, m_lambda, m_mu) / m_area;
    const double psi = Elasticity::NeoHookean(m_mesh, pos, m_restPos, m_lambda, m_mu) / m_area;

    file << m_curvature << " "  << m_angle << " "
        << kVoigt(0) << " " << kVoigt(1) << " " << kVoigt(2) << " " 
        << mVoigt(0) << " " << mVoigt(1) << " " << mVoigt(2) << " "
        << psi << std::endl;
}


Eigen::SparseMatrix<double> Homogenization::WoodburyMatrix(const Eigen::SparseMatrix<double>& massMatrix)
{
    const int n = massMatrix.rows();
    const Eigen::VectorXd e = massMatrix.transpose() * Eigen::VectorXd::Ones(n);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * n);
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
            triplets.emplace_back(3*i+d, d, e(i));
        }
    }
    Eigen::SparseMatrix<double> U(3 * n, 3);
    U.setFromTriplets(triplets.begin(), triplets.end());
    return U;
}
