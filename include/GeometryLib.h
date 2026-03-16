#ifndef GEOMETRY_LIB_H
#define GEOMETRY_LIB_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class GeometryLib
{
public:
    GeometryLib() = delete;

    // project (x,y,z) to cylinder: ((k*z+1)/k*sin(kx), y, ((k*z+1)*cos(kx)-1)/k)
    static Eigen::MatrixX3d cylindricalProject(const Eigen::MatrixX3d& w, double k, double alpha);
    static Eigen::SparseMatrix<double> cylindricalJacobian(const Eigen::MatrixX3d& w, double k, double alpha);
    static Eigen::SparseMatrix<double> cylindricalHessian(const Eigen::MatrixX3d& w, double k, double alpha, const Eigen::VectorXd& gradient);

    // // triangle
    // static double area(const Eigen::Matrix<double,3,2>& X);

    // // Deformation gradient for 3D triangles is only well defined if it is 3x2: map a 2D undeformed vector to 3D deformed vector
    // // Only works for initially flat surfaces. Not sure how to enable orthogonal materials for curved surfaces.
    // static Eigen::Matrix<double,3,2> deformationGradient(const Eigen::Matrix3d& x, 
    //                                                      const Eigen::Matrix<double,3,2>& X, 
    //                                                      Eigen::Matrix<double,6,9>* dFdx=nullptr);
    
    // // shape operator
    // // static double area(const Eigen::RowVector3d &p0, const Eigen::RowVector3d &p1, const Eigen::RowVector3d &p2);
    // static double height(const Eigen::RowVector3d &x0, 
    //                      const Eigen::RowVector3d &x1, 
    //                      const Eigen::RowVector3d &x2,
    //                      Eigen::Matrix<double,9,1>* dhdx=nullptr);
    // static Eigen::Matrix3d crossMatrix(const Eigen::RowVector3d& v);
    // /*
    // *        x2
    // *       /  \
    // *     x0 -- x1
    // *       \  /
    // *        x3
    // * "Discrete Quadratic Curvature Energies" 
    // * https://github.com/evouga/libshell/blob/master/src/SecondFundamentalForm/MidedgeAngleThetaFormulation.cpp#L11
    // */
    // static double exteriorDihedralAngle(const Eigen::RowVector3d &x0, 
    //                                     const Eigen::RowVector3d &x1, 
    //                                     const Eigen::RowVector3d &x2, 
    //                                     const Eigen::RowVector3d &x3,
    //                                     Eigen::Matrix<double,12,1>* dtdx=nullptr,
    //                                     Eigen::Matrix<double,12,12>* d2tdx2=nullptr);
    // /*
    // * From libshell:
    // * Signed angle between two vectors, as measured on the oriented plane with normal parallel to the given axis (which
    // * must be perpendicular to both vectors).
    // * Derivatives are with respect to v, w, a. Note that the derivative with respect to a is always zero (but the function
    // * fill a 1x9 vector for consistency with the Hessian, which *does* have non-zero blocks with respect to a.
    // */
    // static double angle(const Eigen::RowVector3d& v, 
    //                     const Eigen::RowVector3d& w, 
    //                     const Eigen::RowVector3d& axis,
    //                     Eigen::Matrix<double,9,1>* dadv=nullptr,
    //                     Eigen::Matrix<double,9,9>* d2adv2=nullptr);
    // static Eigen::Vector3d shapeOperator(const std::array<Eigen::RowVector3d, 6> &ps, 
    //                                      const std::array<Eigen::RowVector2d, 3> &ts,
    //                                      Eigen::Matrix<double,3,18>* dkdx=nullptr);

};

#endif