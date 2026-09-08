#ifndef ELEMENT_TET10_H
#define ELEMENT_TET10_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

// Quadratic tetrahedra
class Tet10
{    
public:
    Tet10(const Eigen::Matrix<int, -1, 10>& T);

    int nTets() const { return static_cast<int>(m_T.rows()); }

    Eigen::Array<int, 10, 1> indices(int i) const { return m_T.row(i).array(); }

    std::vector<double> volumes(const Eigen::MatrixX3d& X) const;

    double restZ(int i, const std::array<double, 3>& q, const Eigen::MatrixX3d& X) const;

    static Eigen::Matrix<double, 10, 1> basisFunction(
        const std::array<double, 3>& q, 
        Eigen::Matrix<double, 10, 3>* dNdq=nullptr);

    Eigen::SparseMatrix<double> massMatrix(const Eigen::MatrixX3d& V) const;

    Eigen::Matrix3d deformationGradient(
        int i,
        const std::array<double, 3>& q,
        const Eigen::MatrixX3d& x,
        const Eigen::MatrixX3d& X,
        Eigen::Matrix<double, 9, 30>* dFdxi=nullptr) const;

    Eigen::Matrix<double, 6, 1> GreenStrainVoigt(
        int i,
        const std::array<double, 3>& q,
        const Eigen::MatrixX3d& x,
        const Eigen::MatrixX3d& X,
        Eigen::Matrix<double, 6, 30>* dEdxi=nullptr,
        std::array<Eigen::Matrix<double, 30, 30>, 6>* d2EdF2=nullptr) const;


    const Eigen::Matrix<int, -1, 10> &m_T; // tetrahedra list
    Eigen::MatrixX3i m_F; // faces for visualization

    // quadrature points
    // https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_et1.html#b66e328lmm
    // Zienkiewicz, O.C. and Taylor, R.L., The finite element method.McGraw-Hill Book Company (1989).
    static constexpr std::array<std::pair<double, std::array<double, 3>>, 4> m_sQuadrature = 
        std::to_array<std::pair<double, std::array<double, 3>>>({
            {0.25, {0.58541020, 0.13819660, 0.13819660}},
            {0.25, {0.13819660, 0.58541020, 0.13819660}},
            {0.25, {0.13819660, 0.13819660, 0.58541020}},
            {0.25, {0.13819660, 0.13819660, 0.13819660}}
        });
};

#endif