#ifndef SHELL_HELPER_H
#define SHELL_HELPER_H

#include <Eigen/Core>

class ShellHelper
{
public:
    ShellHelper() = delete;

    // Deformation gradient for 3D triangles is only well defined if it is 3x2: map a 2D undeformed vector to 3D deformed vector
    // Only works for initially flat surfaces. Not sure how to enable orthogonal materials for curved surfaces.
    static Eigen::Matrix<double, 3, 2> deformationGradient(const Eigen::Matrix3d& x, 
                                                           const Eigen::Matrix<double, 3, 2>& X, 
                                                           Eigen::Matrix<double, 6, 9>* dFdx=nullptr);
    static Eigen::Vector3d GreenStrainVoigt(const Eigen::Matrix3d& x, 
                                            const Eigen::Matrix<double, 3, 2>& X, 
                                            Eigen::Matrix<double, 3, 9>* dEdx=nullptr,
                                            std::array<Eigen::Matrix<double, 9, 9>, 3>* d2Edx2=nullptr);
    static double stretchingEnergy(const Eigen::Matrix3d& x, 
                                   const Eigen::Matrix<double, 3, 2>& X, 
                                   const Eigen::Matrix3d& C, double V,
                                   Eigen::Matrix<double, 9, 1>* dedx=nullptr,
                                   Eigen::Matrix<double, 9, 9>* d2edx2=nullptr);


    static Eigen::Matrix3d crossMatrix(const Eigen::RowVector3d& v);
    /*
    * From libshell:
    * Signed angle between two vectors, as measured on the oriented plane with normal parallel to the given axis (which
    * must be perpendicular to both vectors).
    * Derivatives are with respect to v, w, a. Note that the derivative with respect to a is always zero (but the function
    * fill a 1x9 vector for consistency with the Hessian, which *does* have non-zero blocks with respect to a.
    */
    static double angle(const Eigen::RowVector3d& v, 
                        const Eigen::RowVector3d& w, 
                        const Eigen::RowVector3d& axis,
                        Eigen::Matrix<double, 9, 1>* dadv=nullptr,
                        Eigen::Matrix<double, 9, 9>* d2adv2=nullptr);
    /*
    *        x2
    *       /  \
    *     x0 -- x1
    *       \  /
    *        x3
    * "Discrete Quadratic Curvature Energies" 
    * https://github.com/evouga/libshell/blob/master/src/SecondFundamentalForm/MidedgeAngleThetaFormulation.cpp#L11
    */
    static double exteriorDihedralAngle(const Eigen::RowVector3d &x0, 
                                        const Eigen::RowVector3d &x1, 
                                        const Eigen::RowVector3d &x2, 
                                        const Eigen::RowVector3d &x3,
                                        Eigen::Matrix<double, 12, 1>* dtdx=nullptr,
                                        Eigen::Matrix<double, 12, 12>* d2tdx2=nullptr);
    static double height(const Eigen::RowVector3d &x0, 
                         const Eigen::RowVector3d &x1, 
                         const Eigen::RowVector3d &x2,
                         Eigen::Matrix<double, 9, 1>* dhdx=nullptr,
                         Eigen::Matrix<double, 9, 9>* d2hdx2=nullptr);
    // discrete shape operator (https://cims.nyu.edu/gcl/papers/grinspun2006cds.pdf)
    static Eigen::Vector3d shapeOperator(const Eigen::Matrix3d &ps, 
                                         const Eigen::Matrix3d &qs, 
                                         const std::array<Eigen::Vector2d, 3> &ts,
                                         Eigen::Matrix<double, 3, 18>* dkdx=nullptr,
                                         std::array<Eigen::Matrix<double, 18, 18>, 3>* d2kdx2=nullptr);
    static double bendingEnergy(const Eigen::Matrix3d &ps, 
                                const Eigen::Matrix3d &qs, 
                                const std::array<Eigen::Vector2d, 3> &ts,
                                double A, const Eigen::Matrix3d& B,
                                Eigen::Matrix<double, 18, 1>* dedx=nullptr,
                                Eigen::Matrix<double, 18, 18>* d2edx2=nullptr);

};

#endif