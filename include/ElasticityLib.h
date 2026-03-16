#ifndef ELASTICITY_LIB_H
#define ELASTICITY_LIB_H

#include <Eigen/Core>

// result is in column major
class ElasticityLib
{    
public:
    ElasticityLib() = delete;
    // (Multi-Layer Thick Shells: Supplemental Document)
    static Eigen::MatrixXd deformationGradient3DGradient(const Eigen::MatrixXd &B);
    static Eigen::Matrix<double,6,9> GreenStrainGradient(const Eigen::Matrix3d &F);

    // // for triangles
    // static Eigen::Vector3d GreenStrainVoigt(const Eigen::Matrix3d& x, 
    //                                         const Eigen::Matrix<double,3,2>& X, 
    //                                         Eigen::Matrix<double,3,9>* dEdx=nullptr,
    //                                         std::array<Eigen::Matrix<double,9,9>,3>* d2Edx2=nullptr);
    // // StVK
    // static double shellStVK(const Eigen::Matrix3d& x, 
    //                         const Eigen::Matrix<double,3,2>& X, 
    //                         const Eigen::Matrix3d& C, double h,
    //                         Eigen::Matrix<double,9,1>* dedx=nullptr,
    //                         Eigen::Matrix<double,9,9>* d2edx2=nullptr);
};

#endif