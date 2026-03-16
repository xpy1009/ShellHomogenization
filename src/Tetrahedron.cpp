#include <Tetrahedron.h>
#include <ElasticityLib.h>

#include <igl/boundary_facets.h>
#include <oneapi/tbb.h>
#include <functional>
#include <exception>

Tetrahedron::Tetrahedron(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T) : undeformed(V), deformed(V), tets(T)
{
    if (V.cols() != 3) {
        throw std::runtime_error("Tets not in 3D?");
    }
    if (T.cols() != 10) {
        throw std::runtime_error("Only quadratic tets!");
    }

    // build linearized boundary faces for visualization
    Eigen::MatrixXi _T = tets.leftCols(4);
    Eigen::MatrixXi F;
    Eigen::VectorXi J, K;
    igl::boundary_facets(_T, F, J, K);
    faces.resize(4*F.rows(),3);
    for (int i = 0; i < J.size(); i++) {
        const int tid = J(i), vid = K(i);
        switch (vid) {
        case 0:
            faces.middleRows(4*i,4) << tets(tid,1),tets(tid,5),tets(tid,8),
                                    tets(tid,5),tets(tid,2),tets(tid,9),
                                    tets(tid,8),tets(tid,9),tets(tid,3),
                                    tets(tid,8),tets(tid,5),tets(tid,9);
            break;
        case 1:
            faces.middleRows(4*i,4) << tets(tid,0),tets(tid,7),tets(tid,6),
                                    tets(tid,3),tets(tid,9),tets(tid,7),
                                    tets(tid,9),tets(tid,2),tets(tid,6),
                                    tets(tid,7),tets(tid,9),tets(tid,6);
            break;
        case 2:
            faces.middleRows(4*i,4) << tets(tid,0),tets(tid,4),tets(tid,7),
                                    tets(tid,4),tets(tid,1),tets(tid,8),
                                    tets(tid,7),tets(tid,8),tets(tid,3),
                                    tets(tid,4),tets(tid,8),tets(tid,7);
            break;
        case 3:
            faces.middleRows(4*i,4) << tets(tid,4),tets(tid,0),tets(tid,6),
                                    tets(tid,1),tets(tid,4),tets(tid,5),
                                    tets(tid,5),tets(tid,6),tets(tid,2),
                                    tets(tid,4),tets(tid,6),tets(tid,5);
            break;
        default:
            throw std::runtime_error("Face for tets");
        }
    }

    // https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_et1.html#b66e328lmm
    // Zienkiewicz, O.C. and Taylor, R.L., The finite element method.McGraw-Hill Book Company (1989).
    qs = {{0.25, {0.58541020, 0.13819660, 0.13819660}},
        {0.25, {0.13819660, 0.58541020, 0.13819660}},
        {0.25, {0.13819660, 0.13819660, 0.58541020}},
        {0.25, {0.13819660, 0.13819660, 0.13819660}}};

    // double E = 1e9, nu = 0.4;
    double E = 1e6, nu = 0.4;
    m_lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    m_mu = E / 2 / (1 + nu);
    C << 2 * m_mu + m_lambda, m_lambda, m_lambda, 0, 0, 0,
        m_lambda, 2 * m_mu + m_lambda, m_lambda, 0, 0, 0,
        m_lambda, m_lambda,  2 * m_mu + m_lambda, 0, 0, 0,
        0, 0, 0, m_mu, 0, 0,
        0, 0, 0, 0, m_mu, 0,
        0, 0, 0, 0, 0, m_mu;

}

Eigen::Matrix<double, 10, 1> Tetrahedron::basisFunction(const Eigen::Vector3d &q) const
{
    Eigen::Matrix<double, 10, 1> N;
    N << (1 - q(0) - q(1) - q(2)) * (2 * (1 - q(0) - q(1) - q(2)) - 1), 
        q(0) * (2 * q(0) - 1), 
        q(1) * (2 * q(1) - 1),
        q(2) * (2 * q(2) - 1),
        4 * q(0) * (1 - q(0) - q(1) - q(2)),
        4 * q(0) * q(1),
        4 * q(1) * (1 - q(0) - q(1) - q(2)),
        4 * q(2) * (1 - q(0) - q(1) - q(2)),
        4 * q(0) * q(2),
        4 * q(1) * q(2);
    return N;
}

Eigen::Matrix<double, 10, 3> Tetrahedron::basisFunctionGradient(const Eigen::Vector3d &q) const
{
    Eigen::Matrix<double, 10, 3> dNdq;
    dNdq << 4*(q(0)+q(1)+q(2))-3, 4*(q(0)+q(1)+q(2))-3, 4*(q(0)+q(1)+q(2))-3, 
            4*q(0)-1, 0, 0, 
            0, 4*q(1)-1, 0, 
            0, 0, 4*q(2)-1,
            4*(1-2*q(0)-q(1)-q(2)), -4*q(0), -4*q(0),
            4*q(1), 4*q(0), 0,
            -4*q(1), 4*(1-q(0)-2*q(1)-q(2)), -4*q(1),
            -4*q(2), -4*q(2), 4*(1-q(0)-q(1)-2*q(2)),
            4*q(2), 0, 4*q(0),
            0, 4*q(2), 4*q(1);
    return dNdq;
}

double Tetrahedron::energy() const
{
    auto elementEnergy = [&](int i) {
        double e = 0.0;
        const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
        const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::placeholders::all);
        const Eigen::Matrix<double, 10, 3> xi = deformed(ind, Eigen::placeholders::all);
        for (int j = 0; j < qs.size(); ++j)
        {
            const double& wj = qs[j].first;
            const Eigen::Vector3d& qj = qs[j].second;
            const Eigen::Matrix<double, 10, 3> dNdqj = basisFunctionGradient(qj);

            const Eigen::Matrix3d dXdq = Xi.transpose() * dNdqj;
            const Eigen::Matrix3d dxdq = xi.transpose() * dNdqj;
            const Eigen::Matrix3d F = dxdq * dXdq.inverse();

            double Psi;
            if (m_model == Model::StVK) {
                const Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
                Eigen::Matrix<double, 6, 1> E_voigt; 
                E_voigt << E(0,0), E(1,1), E(2,2), 2*E(1,2), 2*E(0,2), 2*E(0,1);
                Psi = 0.5 * E_voigt.transpose() * C * E_voigt;
            } else if (m_model == Model::NeoHookean) {
                const double J = F.determinant();
                Psi = 0.5 * m_mu * ((F.transpose() * F).trace() - 3) - m_mu * log(J) + 0.5 * m_lambda * log(J) * log(J);
            }

            const double volume = abs(dXdq.determinant()) / 6.0;
            e += Psi * volume * wj;
        }
        return e;
    };

    return oneapi::tbb::parallel_reduce(
        oneapi::tbb::blocked_range<int>(0, tets.rows()), 0.0,
        [&](const oneapi::tbb::blocked_range<int>& r, double running_total) {
            for(int i = r.begin(); i != r.end(); ++i) {
                running_total += elementEnergy(i);
            }
            return running_total;
        }, std::plus<double>()
    );
}

Eigen::VectorXd Tetrahedron::gradient() const
{
    std::vector<Eigen::Matrix<double, 30, 1>> gradients(tets.rows(), Eigen::Matrix<double, 30, 1>::Zero());
    auto elementGradient = [&](int i) {
        const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
        const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::placeholders::all);
        const Eigen::Matrix<double, 10, 3> xi = deformed(ind, Eigen::placeholders::all);

        Eigen::Matrix<double, 30, 1> gi; gi.setZero();
        for (int j = 0; j < qs.size(); ++j)
        {
            const double& wj = qs[j].first;
            const Eigen::Vector3d& qj = qs[j].second;
            const Eigen::Matrix<double, 10, 3> dNdqj = basisFunctionGradient(qj);

            const Eigen::Matrix3d dXdqj = Xi.transpose() * dNdqj;
            const Eigen::Matrix3d dxdqj = xi.transpose() * dNdqj;
            const Eigen::Matrix<double, 10, 3> B = dNdqj * dXdqj.inverse();
            const Eigen::Matrix<double, 9, 30> dFdx = ElasticityLib::deformationGradient3DGradient(B);

            const Eigen::Matrix3d F = dxdqj * dXdqj.inverse();
            const double volume = abs(dXdqj.determinant()) / 6.0;

            Eigen::Matrix<double, 30, 1> dedx;
            if (m_model == Model::StVK) {
                const Eigen::Matrix<double,6,9> dEdF = ElasticityLib::GreenStrainGradient(F);
                const Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
                Eigen::Matrix<double,6,1> E_voigt; 
                E_voigt << E(0,0), E(1,1), E(2,2), 2*E(1,2), 2*E(0,2), 2*E(0,1);
                const Eigen::Matrix<double, 6, 1> dPsidE = C * E_voigt;
                dedx = volume * dFdx.transpose() * dEdF.transpose() * dPsidE;
            } else if (m_model == Model::NeoHookean) {
                const double J = F.determinant();
                const Eigen::Matrix3d dPsidF = m_mu * (F - F.transpose().inverse()) + m_lambda * log(J) * F.transpose().inverse();
                dedx = volume * dFdx.transpose() * dPsidF.reshaped(9,1);
            }

            gi += wj * dedx;
        }
        gradients[i] = gi;
    };

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<int>(0, tets.rows()),
        [&](const oneapi::tbb::blocked_range<int>& r) {
            for(int i = r.begin(); i != r.end(); ++i) {
                elementGradient(i);
            }
        }
    );

    // assemble
    Eigen::VectorXd gradient(deformed.size()); gradient.setZero();
    for (int i = 0; i < tets.rows(); ++i) {
        for (int j = 0; j < tets.cols(); ++j) {
            gradient.segment(3*tets(i,j), 3) += gradients[i].segment(3*j, 3);
        }
    }

    return gradient;
}

Eigen::SparseMatrix<double> Tetrahedron::hessian() const
{
    std::vector<Eigen::Matrix<double, 30, 30>> hessians(tets.rows(), Eigen::Matrix<double, 30, 30>::Zero());
    auto elementHessian = [&](int i) {
        const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
        const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::placeholders::all);
        const Eigen::Matrix<double, 10, 3> xi = deformed(ind, Eigen::placeholders::all);

        Eigen::Matrix<double, 30, 30> hi; hi.setZero();
        for (int j = 0; j < qs.size(); ++j)
        {
            const double& wj = qs[j].first;
            const Eigen::Vector3d& qj = qs[j].second;
            const Eigen::Matrix<double, 10, 3> dNdqj = basisFunctionGradient(qj);

            const Eigen::Matrix3d dXdq = Xi.transpose() * dNdqj;
            const Eigen::Matrix3d dxdq = xi.transpose() * dNdqj;
            const Eigen::Matrix<double, 10, 3> B = dNdqj * dXdq.inverse();
            const Eigen::Matrix<double, 9, 30> dFdx = ElasticityLib::deformationGradient3DGradient(B);

            const Eigen::Matrix3d F = dxdq * dXdq.inverse();

            Eigen::Matrix<double,9,9> d2PsidF2;
            if (m_model == Model::StVK) {
                const Eigen::Matrix<double,6,9> dEdF = ElasticityLib::GreenStrainGradient(F);
                Eigen::Matrix<double,9,9> d2E1dF2, d2E2dF2, d2E3dF2, d2E4dF2, d2E5dF2, d2E6dF2;
                d2E1dF2.setZero(); d2E1dF2(0,0) = d2E1dF2(1,1) = d2E1dF2(2,2) = 1;
                d2E2dF2.setZero(); d2E2dF2(3,3) = d2E2dF2(4,4) = d2E2dF2(5,5) = 1;
                d2E3dF2.setZero(); d2E3dF2(6,6) = d2E3dF2(7,7) = d2E3dF2(8,8) = 1;
                d2E4dF2.setZero(); d2E4dF2(3,6) = d2E4dF2(6,3) = d2E4dF2(4,7) = d2E4dF2(7,4) = d2E4dF2(5,8) = d2E4dF2(8,5) = 1;
                d2E5dF2.setZero(); d2E5dF2(0,6) = d2E5dF2(6,0) = d2E5dF2(1,7) = d2E5dF2(7,1) = d2E5dF2(2,8) = d2E5dF2(8,2) = 1;
                d2E6dF2.setZero(); d2E6dF2(0,3) = d2E6dF2(3,0) = d2E6dF2(1,4) = d2E6dF2(4,1) = d2E6dF2(2,5) = d2E6dF2(5,2) = 1;

                const Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
                Eigen::Matrix<double,6,1> E_voigt; 
                E_voigt << E(0,0), E(1,1), E(2,2), 2*E(1,2), 2*E(0,2), 2*E(0,1);
                const Eigen::Matrix<double,6,1> dPsidE = C * E_voigt;
                const Eigen::Matrix<double,6,6> d2PsidE2 = C;

                d2PsidF2 = dEdF.transpose() * d2PsidE2 * dEdF + 
                        dPsidE(0)*d2E1dF2 + dPsidE(1)*d2E2dF2 + dPsidE(2)*d2E3dF2 + 
                        dPsidE(3)*d2E4dF2 + dPsidE(4)*d2E5dF2 + dPsidE(5)*d2E6dF2;
            } else if (m_model == Model::NeoHookean) {
                const Eigen::Matrix3d Finv = F.inverse();
                const Eigen::Matrix3d FinvT = F.inverse().transpose();
                const double J = F.determinant();
                Eigen::Matrix<double,9,9> dFinvTdF; 
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {
                            for (int l = 0; l < 3; l++) {
                                dFinvTdF(3*k+l, i+3*j) = - Finv(k,i) * Finv(j,l);
                            }
                        }
                    }
                }

                const Eigen::Matrix<double,9,1> FinvTv = FinvT.reshaped(); 
                Eigen::Matrix<double,9,9> dlnJFinvTdF;
                for (int i = 0; i < 9; i++) {
                    dlnJFinvTdF.row(i) = log(J)*dFinvTdF.row(i) + FinvTv(i)*FinvTv.transpose();
                }

                d2PsidF2 = m_mu * (Eigen::Matrix<double,9,9>::Identity() - dFinvTdF) + m_lambda * dlnJFinvTdF;
            }

            const double volume = abs(dXdq.determinant()) / 6.0;
            const Eigen::Matrix<double, 30, 30> d2edx2 = wj * volume * dFdx.transpose() * d2PsidF2 * dFdx;
            hi += d2edx2;
        }
        hessians[i] = hi;
    };

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<int>(0, tets.rows()),
        [&](const oneapi::tbb::blocked_range<int>& r) {
            for(int i = r.begin(); i != r.end(); ++i) {
                elementHessian(i);
            }
        }
    );

    // assemble
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(tets.rows() * tets.cols() * tets.cols() * 3 * 3);
    for (int i = 0; i < tets.rows(); ++i) {
        for (int j = 0; j < tets.cols(); ++j) {
            for (int k = 0; k < tets.cols(); ++k) {
                for (int d = 0; d < 3; ++d) {
                    for (int dd = 0; dd < 3; ++dd) {
                        triplets.emplace_back(3*tets(i,j)+d, 3*tets(i,k)+dd, hessians[i](3*j+d, 3*k+dd));                           
                    }
                }
            }
        }
    }

    Eigen::SparseMatrix<double> hessian(deformed.size(), deformed.size());
    hessian.setFromTriplets(triplets.begin(), triplets.end());

    return hessian;
}

Eigen::SparseMatrix<double> Tetrahedron::massMatrix() const
{
    // NUMPDE 2.7.5 Local Computations 2.7.5.5
    Eigen::Matrix<double, 10, 10> eleMat; 
    eleMat << 6,1,1,1,-4,-6,-4,-4,-6,-6,
            1,6,1,1,-4,-4,-6,-6,-4,-6,
            1,1,6,1,-6,-4,-4,-6,-6,-4,
            1,1,1,6,-6,-6,-6,-4,-4,-4,
            -4,-4,-6,-6,32,16,16,16,16,8,
            -6,-4,-4,-6,16,32,16,8,16,16,
            -4,-6,-4,-6,16,16,32,16,8,16,
            -4,-6,-6,-4,16,8,16,32,16,16,
            -6,-4,-6,-4,16,16,8,16,32,16,
            -6,-6,-4,-4,8,16,16,16,16,32;
    eleMat /= 420.0;

    std::vector<Eigen::Triplet<double>> triplets;
    const int nV = tets.cols();
    triplets.reserve(tets.rows() * nV * nV);
    for (int i = 0; i < tets.rows(); ++i) {
        const double vi = elementVolume(i);
        Eigen::Matrix<double, 10, 10> eleMat_i = vi * eleMat;
        
        for (int j = 0; j < nV; ++j) {
            for (int k = 0; k < nV; ++k) {
                triplets.emplace_back(tets(i,j), tets(i,k), eleMat_i(j,k));   
            }
        }
    }
    const int n = undeformed.rows();
    Eigen::SparseMatrix<double> massMat(n, n);
    massMat.setFromTriplets(triplets.begin(), triplets.end());

    double volume = 0.0;
    for(int i = 0; i < tets.rows(); ++i) {
        volume += elementVolume(i);
    }
    massMat /= volume;

    // check
    Eigen::VectorXd _V(n); _V.setOnes();
    const double value = (_V.transpose() * massMat * _V).value();
    if (abs(value - 1) > 1e-8) {
        throw std::runtime_error("Tet mass matrix check fails.");
    }
    
    return massMat;
}

double Tetrahedron::elementVolume(int i) const
{
    const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
    const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::placeholders::all);
    const Eigen::RowVector3d X0 = Xi.row(0), X1 = Xi.row(1), X2 = Xi.row(2), X3 = Xi.row(3);

    Eigen::Matrix3d F;
    F << X1 - X0, X2 - X0, X3 - X0;

    return abs(F.determinant()) / 6.0;
}

Eigen::Vector3d Tetrahedron::moment() const
{
    Eigen::Vector3d mTotal; mTotal.setZero();
    for (int i = 0; i < tets.rows(); ++i)
    {
        const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
        const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::placeholders::all);
        const Eigen::Matrix<double, 10, 3> xi = deformed(ind, Eigen::placeholders::all);

        for (const auto &[wj, qj] : qs) {
            const Eigen::Matrix<double, 10, 3> dNdqj = basisFunctionGradient(qj);
            const Eigen::Matrix3d dXdq = Xi.transpose() * dNdqj;
            const Eigen::Matrix3d dxdq = xi.transpose() * dNdqj;
            const Eigen::Matrix3d F = dxdq * dXdq.inverse();
            const Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
            Eigen::Matrix<double, 6, 1> E_voigt; E_voigt << E(0,0), E(1,1), E(2,2), 2*E(1,2), 2*E(0,2), 2*E(0,1);
            const Eigen::Matrix<double, 6, 1> sigmaVoigt = C * E_voigt;
            const Eigen::Vector3d sigma2D(sigmaVoigt(0), sigmaVoigt(1), sigmaVoigt(5));

            const Eigen::Matrix<double, 10, 1> Nj = basisFunction(qj);
            const Eigen::Vector3d Vj = Xi.transpose() * Nj;
            const double z = Vj.z(); 

            const double volume = abs(dXdq.determinant()) / 6.0;

            mTotal += sigma2D * z * volume * wj;
        }
    }
    return mTotal;
}

double Tetrahedron::volume() const
{
    double volume = 0.0;
    for(int i = 0; i < tets.rows(); ++i) {
        volume += elementVolume(i);
    }
    return volume;
}
