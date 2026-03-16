#include <Material/Shell.h>
#include <Material/ShellHelper.h>

#include <oneapi/tbb.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/doublearea.h>

Shell::Shell(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F) : m_v(V), m_F(F), m_V(V.leftCols(2))
{
    constexpr double eps = 1e-8;
    if (abs(V.col(2).maxCoeff() - V.col(2).minCoeff()) > eps) {
        throw std::runtime_error("Only flat triangles!");
    }

    igl::doublearea(m_V, m_F, m_A);
    m_A *= 0.5;

    // build flaps
    Eigen::MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    m_flaps.resize(F.rows());
    m_ts.resize(F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        const Eigen::RowVector2d p0 = m_V.row(F(i, 0));
        const Eigen::RowVector2d p1 = m_V.row(F(i, 1));
        const Eigen::RowVector2d p2 = m_V.row(F(i, 2));
        const std::array<Eigen::RowVector2d, 3> vs = {p2 - p1, p0 - p2, p1 - p0};
        for (int j = 0; j < 3; ++j) {
            const int eid = (j + 1) % 3;
            const int fid = TT(i, eid);
            if (fid < 0) {
                m_flaps[i][j] = F(i, j);
                m_ts[i][j].setZero();
            } else {
                m_flaps[i][j] = F(fid, (TTi(i, eid) + 2) % 3);
                m_ts[i][j] = Eigen::Vector2d(vs[j].y(), -vs[j].x()).normalized(); 
            }
        }
    }

    // set materials
    m_h = 1e-1;
    constexpr double E = 1e6, nu = 0.4;
    m_C << 1.0, nu, 0.0, 
            nu, 1.0, 0.0, 
            0.0, 0.0, 0.5 * (1.0 - nu);
    m_C *= E / (1.0 - nu * nu);
    m_B = m_C * m_h * m_h * m_h / 12.0;
}

double Shell::energyStretch(int i, Eigen::Matrix<double, 9, 1>* grad, Eigen::Matrix<double, 9, 9>* hess) const
{
    const Eigen::Array<int, 3, 1> ind = m_F.row(i).array();
    const Eigen::Matrix<double, 3, 2>& Xi = m_V(ind, Eigen::indexing::all);
    const Eigen::Matrix3d& xi = m_v(ind, Eigen::indexing::all);
    const double Vi = m_A(i) * m_h;
    return ShellHelper::stretchingEnergy(xi, Xi, m_C, Vi, grad, hess);
}

double Shell::energyBend(int i, Eigen::Matrix<double, 18, 1>*grad, Eigen::Matrix<double, 18, 18>*hess) const
{
    const Eigen::Array<int, 3, 1> face = m_F.row(i).array();
    const Eigen::Matrix3d& pi = m_v(face, Eigen::indexing::all);
    const std::array<int, 3>& flap = m_flaps[i];
    const Eigen::Matrix3d& qi = m_v(flap, Eigen::indexing::all);
    const std::array<Eigen::Vector2d, 3>& ti = m_ts[i];
    return ShellHelper::bendingEnergy(pi, qi, ti, m_A(i), m_B, grad, hess);
}

double Shell::energy() const
{
    // double e = 0.0;
    // for (int i = 0; i < m_F.rows(); ++i) {
    //     e += energyStretch(i) + energyBend(i);
    // }
    // return e;

    return oneapi::tbb::parallel_reduce(
        oneapi::tbb::blocked_range<int>(0, m_F.rows()), 0.0,
        [&](const oneapi::tbb::blocked_range<int>& r, double running_total) {
            for(int i = r.begin(); i != r.end(); ++i) {
                running_total += energyStretch(i) + energyBend(i);
            }
            return running_total;
        }, std::plus<double>()
    );
}

Eigen::VectorXd Shell::gradient() const
{
    std::vector<Eigen::Matrix<double, 9, 1>> gradStretch(m_F.rows(), Eigen::Matrix<double, 9, 1>::Zero());
    std::vector<Eigen::Matrix<double, 18, 1>> gradBend(m_F.rows(), Eigen::Matrix<double, 18, 1>::Zero());

    // for (int i = 0; i < m_F.rows(); ++i) {
    //     energyStretch(i, &gradStretch[i]);
    //     energyBend(i, &gradBend[i]);
    // }

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<int>(0, m_F.rows()),
        [&](const oneapi::tbb::blocked_range<int>& r) {
            for(int i = r.begin(); i != r.end(); ++i) {
                energyStretch(i, &gradStretch[i]);
                energyBend(i, &gradBend[i]);
            }
        }
    );

    // assemble
    Eigen::VectorXd gradient(m_v.size()); gradient.setZero();
    for (int i = 0; i < m_F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            gradient.segment(3 * m_F(i, j), 3) += gradStretch[i].segment(3 * j, 3);
            gradient.segment(3 * m_F(i, j), 3) += gradBend[i].segment(3 * j, 3);
            gradient.segment(3 * m_flaps[i][j], 3) += gradBend[i].segment(9 + 3 * j, 3);

        }
    }

    return gradient;
}

Eigen::SparseMatrix<double> Shell::hessian() const
{
    std::vector<Eigen::Matrix<double, 9, 9>> hessStretch(m_F.rows(), Eigen::Matrix<double, 9, 9>::Zero());
    std::vector<Eigen::Matrix<double, 18, 18>> hessBend(m_F.rows(), Eigen::Matrix<double, 18, 18>::Zero());

    // for (int i = 0; i < m_F.rows(); ++i) {
    //     energyStretch(i, nullptr, &hessStretch[i]);
    //     energyBend(i, nullptr, &hessBend[i]);
    // }

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<int>(0, m_F.rows()),
        [&](const oneapi::tbb::blocked_range<int>& r) {
            for(int i = r.begin(); i != r.end(); ++i) {
                energyStretch(i, nullptr, &hessStretch[i]);
                energyBend(i, nullptr, &hessBend[i]);
            }
        }
    );

    // assemble
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(m_F.rows() * m_F.cols() * m_F.cols() * 3 * 3 * 5);
    for (int i = 0; i < m_F.rows(); ++i) {
        for (int j = 0; j < m_F.cols(); ++j) {
            for (int k = 0; k < m_F.cols(); ++k) {
                for (int d = 0; d < 3; ++d) {
                    for (int dd = 0; dd < 3; ++dd) {
                        triplets.emplace_back(3 * m_F(i, j) + d, 3 * m_F(i, k) + dd, hessStretch[i](3 * j + d, 3 * k + dd));                           
                        triplets.emplace_back(3 * m_F(i, j) + d, 3 * m_F(i, k) + dd, hessBend[i](3 * j + d, 3 * k + dd));                           
                        triplets.emplace_back(3 * m_F(i, j) + d, 3 * m_flaps[i][k] + dd, hessBend[i](3 * j + d, 9 + 3 * k + dd));                           
                        triplets.emplace_back(3 * m_flaps[i][j] + d, 3 * m_F(i, k) + dd, hessBend[i](9 + 3 * j + d, 3 * k + dd));                           
                        triplets.emplace_back(3 * m_flaps[i][j] + d, 3 * m_flaps[i][k] + dd, hessBend[i](9 + 3 * j + d, 9 + 3 * k + dd));                           
                    }
                }
            }
        }
    }

    Eigen::SparseMatrix<double> hessian(m_v.size(), m_v.size());
    hessian.setFromTriplets(triplets.begin(), triplets.end());

    return hessian;
}

// Eigen::SparseMatrix<double> Tetrahedron::massMatrix() const
// {
//     // NUMPDE 2.7.5 Local Computations 2.7.5.5
//     Eigen::Matrix<double, 10, 10> eleMat; 
//     eleMat << 6,1,1,1,-4,-6,-4,-4,-6,-6,
//             1,6,1,1,-4,-4,-6,-6,-4,-6,
//             1,1,6,1,-6,-4,-4,-6,-6,-4,
//             1,1,1,6,-6,-6,-6,-4,-4,-4,
//             -4,-4,-6,-6,32,16,16,16,16,8,
//             -6,-4,-4,-6,16,32,16,8,16,16,
//             -4,-6,-4,-6,16,16,32,16,8,16,
//             -4,-6,-6,-4,16,8,16,32,16,16,
//             -6,-4,-6,-4,16,16,8,16,32,16,
//             -6,-6,-4,-4,8,16,16,16,16,32;
//     eleMat /= 420.0;

//     std::vector<Eigen::Triplet<double>> triplets;
//     const int nV = tets.cols();
//     triplets.reserve(tets.rows() * nV * nV);
//     for (int i = 0; i < tets.rows(); ++i) {
//         const double vi = elementVolume(i);
//         Eigen::Matrix<double, 10, 10> eleMat_i = vi * eleMat;
        
//         for (int j = 0; j < nV; ++j) {
//             for (int k = 0; k < nV; ++k) {
//                 triplets.emplace_back(tets(i,j), tets(i,k), eleMat_i(j,k));   
//             }
//         }
//     }
//     const int n = undeformed.rows();
//     Eigen::SparseMatrix<double> massMat(n, n);
//     massMat.setFromTriplets(triplets.begin(), triplets.end());

//     double volume = 0.0;
//     for(int i = 0; i < tets.rows(); ++i) {
//         volume += elementVolume(i);
//     }
//     massMat /= volume;

//     // check
//     Eigen::VectorXd _V(n); _V.setOnes();
//     const double value = (_V.transpose() * massMat * _V).value();
//     if (abs(value - 1) > 1e-8) {
//         throw std::runtime_error("Tet mass matrix check fails.");
//     }
    
//     return massMat;
// }

// double Tetrahedron::elementVolume(int i) const
// {
//     const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
//     const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::all);
//     const Eigen::RowVector3d X0 = Xi.row(0), X1 = Xi.row(1), X2 = Xi.row(2), X3 = Xi.row(3);

//     Eigen::Matrix3d F;
//     F << X1 - X0, X2 - X0, X3 - X0;

//     return abs(F.determinant()) / 6.0;
// }

// Eigen::Vector3d Tetrahedron::moment() const
// {
//     Eigen::Vector3d mTotal; mTotal.setZero();
//     for (int i = 0; i < tets.rows(); ++i)
//     {
//         const Eigen::Array<int, 10, 1> ind = tets.row(i).array();
//         const Eigen::Matrix<double, 10, 3> Xi = undeformed(ind, Eigen::all);
//         const Eigen::Matrix<double, 10, 3> xi = deformed(ind, Eigen::all);

//         for (const auto &[wj, qj] : qs) {
//             const Eigen::Matrix<double, 10, 3> dNdqj = basisFunctionGradient(qj);
//             const Eigen::Matrix3d dXdq = Xi.transpose() * dNdqj;
//             const Eigen::Matrix3d dxdq = xi.transpose() * dNdqj;
//             const Eigen::Matrix3d F = dxdq * dXdq.inverse();
//             const Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
//             Eigen::Matrix<double, 6, 1> E_voigt; E_voigt << E(0,0), E(1,1), E(2,2), 2*E(1,2), 2*E(0,2), 2*E(0,1);
//             const Eigen::Matrix<double, 6, 1> sigmaVoigt = C * E_voigt;
//             const Eigen::Vector3d sigma2D(sigmaVoigt(0), sigmaVoigt(1), sigmaVoigt(5));

//             const Eigen::Matrix<double, 10, 1> Nj = basisFunction(qj);
//             const Eigen::Vector3d Vj = Xi.transpose() * Nj;
//             const double z = Vj.z(); 

//             const double volume = abs(dXdq.determinant()) / 6.0;

//             mTotal += sigma2D * z * volume * wj;
//         }
//     }
//     return mTotal;
// }

// double Tetrahedron::volume() const
// {
//     double volume = 0.0;
//     for(int i = 0; i < tets.rows(); ++i) {
//         volume += elementVolume(i);
//     }
//     return volume;
// }

// std::array<Eigen::VectorXd, 2> DiscreteShell::principalCurvature() const
// {
//     Eigen::VectorXd k1s(m_faces.rows()), k2s(m_faces.rows());
//     for (int i = 0; i < m_faces.rows(); ++i) {
//         const std::array<Eigen::RowVector3d, 6> ps = {  m_deformed.row(m_faces(i, 0)),
//                                                         m_deformed.row(m_faces(i, 1)), 
//                                                         m_deformed.row(m_faces(i, 2)), 
//                                                         m_deformed.row(m_flaps(i, 0)), 
//                                                         m_deformed.row(m_flaps(i, 1)), 
//                                                         m_deformed.row(m_flaps(i, 2)) };
//         const std::array<Eigen::RowVector2d, 3> _ps = { m_undeformed.row(m_faces(i, 0)), 
//                                                         m_undeformed.row(m_faces(i, 1)), 
//                                                         m_undeformed.row(m_faces(i, 2)) };

//         const Eigen::Vector3d kv = GeometryLib::shapeOperator(ps, _ps, m_isBoundary[i]);
//         Eigen::Matrix2d kt;
//         kt << kv(0), 0.5 * kv(2), kv(1), 0.5 * kv(2);
//         Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(kt);
//         const Eigen::Vector2d ki = eigensolver.eigenvalues();
//         k1s(i) = ki(0);
//         k2s(i) = ki(1);
//     }
//     return {k1s, k2s};
// }