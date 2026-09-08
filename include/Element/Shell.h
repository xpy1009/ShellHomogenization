#ifndef ELEMENT_SHELL_H
#define ELEMENT_SHELL_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class Shell
{
public:
    Shell(const Eigen::MatrixX3i& F);

    int nFaces() const { return m_F.rows(); }

    Eigen::Array3i indices(int i) const { return m_F.row(i).array(); }

    std::array<int, 3> oppositeIndices(int i) const { return m_flaps[i]; }

    Eigen::VectorXd areas(const Eigen::MatrixX2d& X) const;

    Eigen::Matrix<double, 3, 2> deformationGradient(
        int i,
        const Eigen::MatrixX2d& X,
        const Eigen::MatrixX3d& x, 
        Eigen::Matrix<double, 6, 9>* dFdxi=nullptr) const;

    Eigen::Vector3d GreenStrainVoigt(
        int i,
        const Eigen::MatrixX2d& X,
        const Eigen::MatrixX3d& x, 
        Eigen::Matrix<double, 3, 9>* dEdF=nullptr,
        std::array<Eigen::Matrix<double, 9, 9>, 3>* d2EdF2=nullptr) const;


    // based on libshell https://github.com/evouga/libshell and 
    // triangle-averaged shape operator from Computing Discrete Shape Operators on General Meshes
    static Eigen::Matrix3d crossMatrix(const Eigen::RowVector3d& v);
    static double height(
        const Eigen::RowVector3d &x0, 
        const Eigen::RowVector3d &x1, 
        const Eigen::RowVector3d &x2,
        Eigen::Matrix<double, 9, 1>* dhdx=nullptr,
        Eigen::Matrix<double, 9, 9>* d2hdx2=nullptr);
    static double angle(
        const Eigen::RowVector3d& v, 
        const Eigen::RowVector3d& w, 
        const Eigen::RowVector3d& axis,
        Eigen::Matrix<double, 9, 1>* dadv=nullptr,
        Eigen::Matrix<double, 9, 9>* d2adv2=nullptr);
    static double exteriorDihedralAngle(
        const Eigen::RowVector3d &x0, 
        const Eigen::RowVector3d &x1, 
        const Eigen::RowVector3d &x2, 
        const Eigen::RowVector3d &x3,
        Eigen::Matrix<double, 12, 1>* dtdx=nullptr,
        Eigen::Matrix<double, 12, 12>* d2tdx2=nullptr);
    Eigen::Vector3d shapeOperator(
        int idx,
        const Eigen::MatrixX2d& X,
        const Eigen::MatrixX3d& x, 
        Eigen::Matrix<double, 3, 18>* dkdxi=nullptr,
        std::array<Eigen::Matrix<double, 18, 18>, 3>* d2kdxi2=nullptr) const;

    // to decide strains for native scale simulation
    std::array<double, 2> maxStrain(const Eigen::MatrixX2d& X, const Eigen::MatrixX3d& x) const;
    std::array<double, 2> maxCurvature(const Eigen::MatrixX2d& X, const Eigen::MatrixX3d& x) const;


    const Eigen::MatrixX3i& m_F;
    std::vector<std::array<int, 3>> m_flaps;
};

#endif