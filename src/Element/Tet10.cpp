#include <Element/Tet10.h>

#include <Eigen/Dense>
#include <igl/boundary_facets.h>


Tet10::Tet10(const Eigen::Matrix<int, -1, 10> &T) : m_T(T)
{
    // Build linearized boundary faces for visualization
    const Eigen::MatrixX4i tet4 = m_T.leftCols(4);
    Eigen::MatrixXi F;
    Eigen::VectorXi J, K;
    igl::boundary_facets(tet4, F, J, K);
    m_F.resize(4 * F.rows(), 3);
    for (int i = 0; i < J.size(); ++i) {
        const int tid = J(i), vid = K(i);
        switch (vid) {
        case 0:
            m_F.middleRows(4 * i, 4) << m_T(tid, 1), m_T(tid, 5), m_T(tid, 8),
                                        m_T(tid, 5), m_T(tid, 2), m_T(tid, 9),
                                        m_T(tid, 8), m_T(tid, 9), m_T(tid, 3),
                                        m_T(tid, 8), m_T(tid, 5), m_T(tid, 9);
            break;
        case 1:
            m_F.middleRows(4 * i, 4) << m_T(tid, 0), m_T(tid, 7), m_T(tid, 6),
                                        m_T(tid, 3), m_T(tid, 9), m_T(tid, 7),
                                        m_T(tid, 9), m_T(tid, 2), m_T(tid, 6),
                                        m_T(tid, 7), m_T(tid, 9), m_T(tid, 6);
            break;
        case 2:
            m_F.middleRows(4 * i, 4) << m_T(tid, 0), m_T(tid, 4), m_T(tid, 7),
                                        m_T(tid, 4), m_T(tid, 1), m_T(tid, 8),
                                        m_T(tid, 7), m_T(tid, 8), m_T(tid, 3),
                                        m_T(tid, 4), m_T(tid, 8), m_T(tid, 7);
            break;
        case 3:
            m_F.middleRows(4 * i, 4) << m_T(tid, 4), m_T(tid, 0), m_T(tid, 6),
                                        m_T(tid, 1), m_T(tid, 4), m_T(tid, 5),
                                        m_T(tid, 5), m_T(tid, 6), m_T(tid, 2),
                                        m_T(tid, 4), m_T(tid, 6), m_T(tid, 5);
            break;
        default:
            throw std::runtime_error("Fail to build faces for tets");
        }
    }

}

std::vector<double> Tet10::volumes(const Eigen::MatrixX3d& X) const
{
    std::vector<double> vs(m_T.rows());
    for (int i = 0; i < m_T.rows(); ++i) {
        const Eigen::Matrix<double, 10, 3> Xi = X(m_T.row(i), Eigen::placeholders::all);
        // undeformed volumes of quadratic and linear elements should be the same
        const Eigen::RowVector3d X0 = Xi.row(0), X1 = Xi.row(1), X2 = Xi.row(2), X3 = Xi.row(3);
        Eigen::Matrix3d F;
        F << X1 - X0, X2 - X0, X3 - X0;
        vs[i] = std::abs(F.determinant()) / 6.0;        
    }
    return vs;
}

double Tet10::restZ(int i, const std::array<double, 3>& q, const Eigen::MatrixX3d& X) const
{
    const Eigen::Matrix<double, 10, 1> N = basisFunction(q);
    const Eigen::Matrix<double, 10, 3> Xi = X(m_T.row(i), Eigen::placeholders::all);
    const Eigen::Vector3d V = Xi.transpose() * N;
    return V.z(); 
}

Eigen::Matrix<double, 10, 1> Tet10::basisFunction(const std::array<double, 3> &q, Eigen::Matrix<double, 10, 3>* dNdq)
{
    Eigen::Matrix<double, 10, 1> N;
    N << (1.0 - q[0] - q[1] - q[2]) * (2.0 * (1.0 - q[0] - q[1] - q[2]) - 1.0), 
        q[0] * (2.0 * q[0] - 1.0), 
        q[1] * (2.0 * q[1] - 1.0),
        q[2] * (2.0 * q[2] - 1.0),
        4.0 * q[0] * (1.0 - q[0] - q[1] - q[2]),
        4.0 * q[0] * q[1],
        4.0 * q[1] * (1.0 - q[0] - q[1] - q[2]),
        4.0 * q[2] * (1.0 - q[0] - q[1] - q[2]),
        4.0 * q[0] * q[2],
        4.0 * q[1] * q[2];

    if (dNdq != nullptr) {
        *dNdq << 4.0 * (q[0] + q[1] + q[2]) - 3.0, 4.0 * (q[0] + q[1] + q[2]) - 3.0, 4.0 * (q[0] + q[1] + q[2]) - 3.0, 
                4.0 * q[0] - 1.0, 0.0, 0.0, 
                0.0, 4.0 * q[1] - 1.0, 0.0, 
                0.0, 0.0, 4.0 * q[2] - 1.0,
                4.0 * (1.0 - 2.0 * q[0] - q[1] - q[2]), -4.0 * q[0], -4.0 * q[0],
                4.0 * q[1], 4.0 * q[0], 0.0,
                -4.0 * q[1], 4.0 * (1.0 - q[0] - 2.0 * q[1] - q[2]), -4.0 * q[1],
                -4.0 * q[2], -4.0 * q[2], 4.0 * (1.0 - q[0] - q[1] - 2.0 * q[2]),
                4.0 * q[2], 0.0, 4.0 * q[0],
                0.0, 4.0 * q[2], 4.0 * q[1];
    }
    return N;
}

Eigen::SparseMatrix<double> Tet10::massMatrix(const Eigen::MatrixX3d& V) const
{
    constexpr int N = 10;
    // NUMPDE 2.7.5 Local Computations 2.7.5.5
    Eigen::Matrix<double, N, N> eleMat; 
    eleMat << 6.0, 1.0, 1.0, 1.0, -4.0, -6.0, -4.0, -4.0, -6.0, -6.0, 
              1.0, 6.0, 1.0, 1.0, -4.0, -4.0, -6.0, -6.0, -4.0, -6.0, 
              1.0, 1.0, 6.0, 1.0, -6.0, -4.0, -4.0, -6.0, -6.0, -4.0, 
              1.0, 1.0, 1.0, 6.0, -6.0, -6.0, -6.0, -4.0, -4.0, -4.0, 
              -4.0, -4.0, -6.0, -6.0, 32.0, 16.0, 16.0, 16.0, 16.0, 8.0, 
              -6.0, -4.0, -4.0, -6.0, 16.0, 32.0, 16.0, 8.0, 16.0, 16.0, 
              -4.0, -6.0, -4.0, -6.0, 16.0, 16.0, 32.0, 16.0, 8.0, 16.0, 
              -4.0, -6.0, -6.0, -4.0, 16.0, 8.0, 16.0, 32.0, 16.0, 16.0, 
              -6.0, -4.0, -6.0, -4.0, 16.0, 16.0, 8.0, 16.0, 32.0, 16.0, 
              -6.0, -6.0, -4.0, -4.0, 8.0, 16.0, 16.0, 16.0, 16.0, 32.0;
    eleMat /= 420.0;

    const int nT = m_T.rows();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nT * N * N);
    double totalV = 0.0;
    const std::vector<double> vs = volumes(V);
    for (int i = 0; i < nT; ++i) {
        totalV += vs[i];
        const Eigen::Matrix<double, N, N> eleMati = vs[i] * eleMat;
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                triplets.emplace_back(m_T(i, j), m_T(i ,k), eleMati(j, k));   
            }
        }
    }
    const int nV = V.rows();
    Eigen::SparseMatrix<double> massMat(nV, nV);
    massMat.setFromTriplets(triplets.begin(), triplets.end());
    massMat /= totalV;

    // check
    Eigen::VectorXd _V(nV); _V.setOnes();
    const double value = (_V.transpose() * massMat * _V).value();
    if (abs(value - 1.0) > 1e-8) {
        throw std::runtime_error("Tet mass matrix check fails.");
    }
    
    return massMat;
}

Eigen::Matrix3d Tet10::deformationGradient(
    int i,
    const std::array<double, 3>& q,
    const Eigen::MatrixX3d& x,
    const Eigen::MatrixX3d& X,
    Eigen::Matrix<double, 9, 30>* dFdxi) const
{
    const Eigen::Matrix<double, 10, 3> xi = x(m_T.row(i), Eigen::placeholders::all);
    const Eigen::Matrix<double, 10, 3> Xi = X(m_T.row(i), Eigen::placeholders::all);

    Eigen::Matrix<double, 10, 3> dNdq;
    basisFunction(q, &dNdq);
    const Eigen::Matrix3d dXdq = Xi.transpose() * dNdq;
    const Eigen::Matrix3d dxdq = xi.transpose() * dNdq;
    const Eigen::Matrix3d F = dxdq * dXdq.inverse();

    if (dFdxi != nullptr) {
        Eigen::Matrix3d I; 
        I.setIdentity();
        const Eigen::Matrix<double, 10, 3> B = dNdq * dXdq.inverse();
        for (int i = 0; i < B.rows(); ++i) {
            for (int j = 0; j < B.cols(); ++j) {
                dFdxi->block(3 * j, 3 * i, 3, 3) = B(i, j) * I;
            }
        }
    }

    return F;
}

Eigen::Matrix<double, 6, 1> Tet10::GreenStrainVoigt(
    int i,
    const std::array<double, 3>& q,
    const Eigen::MatrixX3d& x,
    const Eigen::MatrixX3d& X,
    Eigen::Matrix<double, 6, 30>* dEdxi,
    std::array<Eigen::Matrix<double, 30, 30>, 6>* d2Edxi2) const
{
    Eigen::Matrix<double, 9, 30> dFdxi;
    const Eigen::Matrix3d F = deformationGradient(i, q, x, X, (dEdxi!=nullptr || d2Edxi2!=nullptr) ? &dFdxi : nullptr);

    const Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity());
    Eigen::Matrix<double, 6, 1> E_voigt; 
    E_voigt << E(0, 0), E(1, 1), E(2, 2), 2.0 * E(1, 2), 2.0 * E(0, 2), 2.0 * E(0, 1);

    if (dEdxi != nullptr) {
        Eigen::Matrix<double, 6, 9> dEdF;
        dEdF << F(0,0), F(1,0), F(2,0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, F(0,1), F(1,1), F(2,1), 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, F(0,2), F(1,2), F(2,2),
                0.0, 0.0, 0.0, F(0,2), F(1,2), F(2,2), F(0,1), F(1,1), F(2,1),
                F(0,2), F(1,2), F(2,2), 0.0, 0.0, 0.0, F(0,0), F(1,0), F(2,0),
                F(0,1), F(1,1), F(2,1), F(0,0), F(1,0), F(2,0), 0.0, 0.0, 0.0;
        *dEdxi = dEdF * dFdxi;
    }

    if (d2Edxi2 != nullptr) {
        std::array<Eigen::Matrix<double, 9, 9>, 6> d2EdF2;
        d2EdF2[0].setZero(); d2EdF2[0](0, 0) = d2EdF2[0](1, 1) = d2EdF2[0](2, 2) = 1.0;
        d2EdF2[1].setZero(); d2EdF2[1](3, 3) = d2EdF2[1](4, 4) = d2EdF2[1](5, 5) = 1.0;
        d2EdF2[2].setZero(); d2EdF2[2](6, 6) = d2EdF2[2](7, 7) = d2EdF2[2](8, 8) = 1.0;
        d2EdF2[3].setZero(); d2EdF2[3](3, 6) = d2EdF2[3](6, 3) = d2EdF2[3](4, 7) = d2EdF2[3](7, 4) = d2EdF2[3](5, 8) = d2EdF2[3](8, 5) = 1.0;
        d2EdF2[4].setZero(); d2EdF2[4](0, 6) = d2EdF2[4](6, 0) = d2EdF2[4](1, 7) = d2EdF2[4](7, 1) = d2EdF2[4](2, 8) = d2EdF2[4](8, 2) = 1.0;
        d2EdF2[5].setZero(); d2EdF2[5](0, 3) = d2EdF2[5](3, 0) = d2EdF2[5](1, 4) = d2EdF2[5](4, 1) = d2EdF2[5](2, 5) = d2EdF2[5](5, 2) = 1.0;
        for (size_t j = 0; j < 6; ++j) {
            d2Edxi2->at(j) = dFdxi.transpose() * d2EdF2[j] * dFdxi;
        }
        
    }
    return E_voigt;
}
