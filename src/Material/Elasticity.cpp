#include <Material/Elasticity.h>
#include <Element/Tet10.h>

#include <Eigen/Eigenvalues> 

Eigen::Matrix<double, 6, 6> Elasticity::stiffnessTensor(double E, double nu)
{
    auto [lambda, mu] = lameParameters(E, nu);
    Eigen::Matrix<double, 6, 6> C;
    C << 2.0 * mu + lambda, lambda, lambda, 0.0, 0.0, 0.0,
        lambda, 2.0 * mu + lambda, lambda, 0.0, 0.0, 0.0,
        lambda, lambda,  2.0 * mu + lambda, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, mu, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, mu, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, mu;
    return C;
}

double Elasticity::StVK(
    const Tet10& mesh,
    const Eigen::MatrixX3d& x, 
    const Eigen::MatrixX3d& X, 
    const Eigen::Matrix<double, 6, 6>& C,
    Eigen::VectorXd* gradient,
    Eigen::SparseMatrix<double>* hessian)
{
    double e = 0.0;
    const int nTets = mesh.nTets();
    if (gradient != nullptr) {
        gradient->setZero(X.size());
    }
    std::vector<Eigen::Triplet<double>> triplets;
    if (hessian != nullptr) {
        hessian->resize(X.size(), X.size());
        triplets.reserve(nTets * 30 * 30);
    }


    const std::vector<double> vs = mesh.volumes(X); 
    for (int i = 0; i < nTets; ++i) {
        const Eigen::Array<int, 10, 1> ind = mesh.indices(i);
        Eigen::Matrix<double, 30, 1> gi; gi.setZero();
        Eigen::Matrix<double, 30, 30> hi; hi.setZero();
        for (const auto &[wj, qj] : mesh.m_sQuadrature) {
            
            Eigen::Matrix<double, 6, 30> dEdx;
            std::array<Eigen::Matrix<double, 30, 30>, 6> d2Edx2;
            const Eigen::Matrix<double, 6, 1> E = mesh.GreenStrainVoigt(
                i, qj, x, X, (gradient!=nullptr || hessian!=nullptr) ? &dEdx : nullptr, hessian!=nullptr ? &d2Edx2 : nullptr);

            const double Psi = 0.5 * E.transpose() * C * E;
            
            if (gradient != nullptr) {
                const Eigen::Matrix<double, 6, 1> dPsidE = C * E;
                const Eigen::Matrix<double, 30, 1> dPsidx = dEdx.transpose() * dPsidE;
                gi += dPsidx * vs[i] * wj;
            }

            if (hessian != nullptr) {
                const Eigen::Matrix<double, 6, 6> &d2PsidE2 = C;
                const Eigen::Matrix<double, 6, 1> dPsidE = C * E;
                const Eigen::Matrix<double, 30, 30> d2Psidx2 = dEdx.transpose() * d2PsidE2 * dEdx + 
                    dPsidE(0) * d2Edx2[0] + dPsidE(1) * d2Edx2[1] + dPsidE(2) * d2Edx2[2] + 
                    dPsidE(3) * d2Edx2[3] + dPsidE(4) * d2Edx2[4] + dPsidE(5) * d2Edx2[5];
                hi += wj * vs[i] * d2Psidx2;
            }

            e += Psi * vs[i] * wj;
        }

        if (gradient != nullptr) {
            for (int j = 0; j < 10; ++j) {
                gradient->segment(3 * ind(j), 3) += gi.segment(3 * j, 3);
            }
        }

        if (hessian != nullptr) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    for (int d = 0; d < 3; ++d) {
                        for (int dd = 0; dd < 3; ++dd) {
                            triplets.emplace_back(3 * ind(j) + d, 3 * ind(k) + dd, hi(3 * j + d, 3 * k + dd));                           
                        }
                    }
                }
            }
        }
    }

    if (hessian != nullptr) {
        hessian->setFromTriplets(triplets.begin(), triplets.end());
    }
    
    return e;
}

Eigen::Matrix<double, 6, 1> Elasticity::secondPKStress(
    const Tet10& mesh, 
    const Eigen::MatrixX3d& x,
    const Eigen::MatrixX3d& X, 
    const Eigen::Matrix<double, 6, 6>& C)
{
    const std::vector<double> vs = mesh.volumes(X); 
    Eigen::Matrix<double, 6, 1> S; S.setZero();
    for (int i = 0; i < mesh.nTets(); ++i) {
        for (const auto &[wj, qj] : mesh.m_sQuadrature) {
            const Eigen::Matrix<double, 6, 1> E = mesh.GreenStrainVoigt(i, qj, x, X);
            S += C * E * vs[i] * wj;
        }
    }
    return S;
}

Eigen::Vector3d Elasticity::moment(
    const Tet10& mesh, 
    const Eigen::MatrixX3d& x,
    const Eigen::MatrixX3d& X, 
    const Eigen::Matrix<double, 6, 6>& C)
{
    // Assume midsurface is at z=0
    Eigen::Vector3d mTotal; mTotal.setZero();
    const std::vector<double> vs = mesh.volumes(X); 
    for (int i = 0; i < mesh.nTets(); ++i) {
        for (const auto &[wj, qj] : mesh.m_sQuadrature) {
            const Eigen::Matrix<double, 6, 1> E = mesh.GreenStrainVoigt(i, qj, x, X);
            const Eigen::Matrix<double, 6, 1> S = C * E;
            const Eigen::Vector3d S2D(S(0), S(1), S(5));

            const double z = mesh.restZ(i, qj, X); 

            mTotal += wj * vs[i] * S2D * z;
        }
    }
    return mTotal;
}


std::array<double, 2> Elasticity::lameParameters(double E, double nu)
{
    const double lambda = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu);
    const double mu = E / 2.0 / (1.0 + nu);
    return {lambda, mu};
}

double Elasticity::NeoHookean(
    const Tet10& mesh,
    const Eigen::MatrixX3d& x, 
    const Eigen::MatrixX3d& X, 
    double lambda,
    double mu,
    Eigen::VectorXd* gradient,
    Eigen::SparseMatrix<double>* hessian)
{
    double e = 0.0;
    const int nTets = mesh.nTets();
    if (gradient != nullptr) {
        gradient->setZero(X.size());
    }
    std::vector<Eigen::Triplet<double>> triplets;
    if (hessian != nullptr) {
        hessian->resize(X.size(), X.size());
        triplets.reserve(nTets * 30 * 30);
    }


    const std::vector<double> vs = mesh.volumes(X); 
    for (int i = 0; i < nTets; ++i) {
        const Eigen::Array<int, 10, 1> ind = mesh.indices(i);
        Eigen::Matrix<double, 30, 1> gi; gi.setZero();
        Eigen::Matrix<double, 30, 30> hi; hi.setZero();
        for (const auto &[wj, qj] : mesh.m_sQuadrature) {
            Eigen::Matrix<double, 9, 30> dFdx;
            const Eigen::Matrix3d F = mesh.deformationGradient(
                i, qj, x, X, (gradient!=nullptr || hessian!=nullptr) ? &dFdx : nullptr);
            
            const double J = F.determinant();
            const double Psi = 0.5 * mu * ((F.transpose() * F).trace() - 3) - mu * log(J) + 0.5 * lambda * log(J) * log(J);

            
            if (gradient != nullptr) {
                const Eigen::Matrix3d dPsidF = mu * (F - F.transpose().inverse()) + lambda * log(J) * F.transpose().inverse();
                const Eigen::Matrix<double, 30, 1> dPsidx = dFdx.transpose() * dPsidF.reshaped(9,1);
                gi += dPsidx * vs[i] * wj;
            }

            if (hessian != nullptr) {
                const Eigen::Matrix3d Finv = F.inverse();
                const Eigen::Matrix3d FinvT = Finv.transpose();
                Eigen::Matrix<double,9,9> dFinvTdF; 
                for (int m = 0; m < 3; ++m) {
                    for (int n = 0; n < 3; ++n) {
                        for (int k = 0; k < 3; ++k) {
                            for (int l = 0; l < 3; ++l) {
                                dFinvTdF(3*k+l, m+3*n) = - Finv(k,m) * Finv(n,l);
                            }
                        }
                    }
                }

                const Eigen::Matrix<double, 9, 1> FinvTv = FinvT.reshaped(); 
                Eigen::Matrix<double, 9, 9> dlnJFinvTdF;
                for (int j = 0; j < 9; ++j) {
                    dlnJFinvTdF.row(j) = log(J) * dFinvTdF.row(j) + FinvTv(j) * FinvTv.transpose();
                }

                const Eigen::Matrix<double, 9, 9> d2PsidF2 = mu * (Eigen::Matrix<double,9,9>::Identity() - dFinvTdF) + lambda * dlnJFinvTdF;
                const Eigen::Matrix<double, 30, 30> d2Psidx2 = dFdx.transpose() * d2PsidF2 * dFdx;
                hi += wj * vs[i] * d2Psidx2;
            }

            e += Psi * vs[i] * wj;
        }

        if (gradient != nullptr) {
            for (int j = 0; j < 10; ++j) {
                gradient->segment(3 * ind(j), 3) += gi.segment(3 * j, 3);
            }
        }

        if (hessian != nullptr) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    for (int d = 0; d < 3; ++d) {
                        for (int dd = 0; dd < 3; ++dd) {
                            triplets.emplace_back(3 * ind(j) + d, 3 * ind(k) + dd, hi(3 * j + d, 3 * k + dd));                           
                        }
                    }
                }
            }
        }
    }

    if (hessian != nullptr) {
        hessian->setFromTriplets(triplets.begin(), triplets.end());
    }
    
    return e;
}

Eigen::Matrix<double, 6, 1> Elasticity::secondPKStress(
    const Tet10& mesh, 
    const Eigen::MatrixX3d& x,
    const Eigen::MatrixX3d& X, 
    double lambda,
    double mu)
{
    const std::vector<double> vs = mesh.volumes(X); 
    Eigen::Matrix<double, 6, 1> S; S.setZero();
    for (int i = 0; i < mesh.nTets(); ++i) {
        for (const auto &[wj, qj] : mesh.m_sQuadrature) {
            const Eigen::Matrix3d F = mesh.deformationGradient(i, qj, x, X);
            const double J = F.determinant();
            const Eigen::Matrix3d FTFinv = (F.transpose() * F).inverse();
            const Eigen::Matrix3d Si = mu * (Eigen::Matrix3d::Identity() - FTFinv) + lambda * log(J) * FTFinv;

            Eigen::Matrix<double, 6, 1> Svi;
            Svi << Si(0, 0), Si(1, 1), Si(2, 2), Si(1, 2), Si(0, 2), Si(0, 1);
            S += Svi * vs[i] * wj;
        }
    }
    return S;
}

Eigen::Vector3d Elasticity::moment(
    const Tet10& mesh, 
    const Eigen::MatrixX3d& x,
    const Eigen::MatrixX3d& X, 
    double lambda,
    double mu)
{
    // Assume midsurface is at z=0
    Eigen::Vector3d mTotal; mTotal.setZero();
    const std::vector<double> vs = mesh.volumes(X); 
    for (int i = 0; i < mesh.nTets(); ++i) {
        for (const auto &[wj, qj] : mesh.m_sQuadrature) {
            const Eigen::Matrix3d F = mesh.deformationGradient(i, qj, x, X);
            const double J = F.determinant();
            const Eigen::Matrix3d FTFinv = (F.transpose() * F).inverse();
            const Eigen::Matrix3d Si = mu * (Eigen::Matrix3d::Identity() - FTFinv) + lambda * log(J) * FTFinv;

            const Eigen::Vector3d S2D(Si(0, 0), Si(1, 1), Si(0, 1));

            const double z = mesh.restZ(i, qj, X); 

            mTotal += wj * vs[i] * S2D * z;
        }
    }
    return mTotal;
}