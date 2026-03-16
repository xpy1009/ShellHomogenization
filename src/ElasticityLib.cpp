#include <ElasticityLib.h>
#include <GeometryLib.h>
#include <iostream>

Eigen::MatrixXd ElasticityLib::deformationGradient3DGradient(const Eigen::MatrixXd &B)
{
    Eigen::Matrix3d I; 
    I.setIdentity();
    Eigen::MatrixXd dFdx(9, 3 * B.rows());
    for (int i = 0; i < B.rows(); ++i) {
        for (int j = 0; j < B.cols(); ++j) {
            dFdx.block(3*j, 3*i, 3, 3) = B(i, j) * I;
        }
    }
    return dFdx;
}

Eigen::Matrix<double,6,9> ElasticityLib::GreenStrainGradient(const Eigen::Matrix3d &F)
{
    Eigen::Matrix<double,6,9> dEdF;
    dEdF << F(0,0),F(1,0),F(2,0),0,0,0,0,0,0,
            0,0,0,F(0,1),F(1,1),F(2,1),0,0,0,
            0,0,0,0,0,0,F(0,2),F(1,2),F(2,2),
            0,0,0,F(0,2),F(1,2),F(2,2),F(0,1),F(1,1),F(2,1),
            F(0,2),F(1,2),F(2,2),0,0,0,F(0,0),F(1,0),F(2,0),
            F(0,1),F(1,1),F(2,1),F(0,0),F(1,0),F(2,0),0,0,0;
    return dEdF;
}

// Eigen::Vector3d ElasticityLib::GreenStrainVoigt(const Eigen::Matrix3d& x, 
//                                                 const Eigen::Matrix<double,3,2>& X,
//                                                 Eigen::Matrix<double,3,9>* dEdx,
//                                                 std::array<Eigen::Matrix<double,9,9>,3>* d2Edx2)
// {
//     Eigen::Matrix<double,6,9> dFdx;
// 	const Eigen::Matrix<double,3,2> F = GeometryLib::deformationGradient(x, X, (dEdx==nullptr && d2Edx2==nullptr) ? nullptr : &dFdx);
// 	const Eigen::Matrix2d E = 0.5 * (F.transpose() * F - Eigen::Matrix2d::Identity()); // Green strain

//     if (dEdx != nullptr) {
//         // F in column major
//         Eigen::Matrix<double,3,6> dEdF;
//         dEdF << F.col(0).transpose(), 0,0,0,
//                 0,0,0, F.col(1).transpose(),
//                 F.col(1).transpose(), F.col(0).transpose();
//         *dEdx = dEdF * dFdx;
//     }

//     if (d2Edx2 != nullptr) {
//         // dE1dF = I 0 \\ 0 0 ...
//         const Eigen::Matrix<double,3,9> tmp1 = dFdx.topRows(3);
//         const Eigen::Matrix<double,3,9> tmp2 = dFdx.bottomRows(3);
//         d2Edx2->at(0) = tmp1.transpose() * tmp1;
//         d2Edx2->at(1) = tmp2.transpose() * tmp2;
//         d2Edx2->at(2) = tmp1.transpose() * tmp2 + tmp2.transpose() * tmp1;
//     }

//     return Eigen::Vector3d(E(0,0), E(1,1), 2.0 * E(0,1)); // Voigt notation
// }

// double ElasticityLib::shellStVK(const Eigen::Matrix3d& x, 
//                                 const Eigen::Matrix<double,3,2>& X, 
//                                 const Eigen::Matrix3d& C, double h,
//                                 Eigen::Matrix<double,9,1>* dedx,
//                                 Eigen::Matrix<double,9,9>* d2edx2)
// {
//     Eigen::Matrix<double,3,9> dEdx;
//     std::array<Eigen::Matrix<double,9,9>,3> d2Edx2;
//     const Eigen::Vector3d EVoigt = GreenStrainVoigt(x, X, 
//                                                     dedx==nullptr && d2edx2==nullptr ? nullptr : &dEdx,
//                                                     d2edx2==nullptr ? nullptr : &d2Edx2);
// 	const double psi = 0.5 * EVoigt.transpose() * C * EVoigt; // energy density
//     const double V_ = GeometryLib::area(X) * h; // undeformed volume
//     const double energy = psi * V_;

//     if (dedx != nullptr) {
//         const Eigen::Vector3d dedE = V_ * C * EVoigt;
//         *dedx = dEdx.transpose() * dedE;
//     }

//     if (d2edx2 != nullptr) {
//         d2edx2->setZero();
//         // dEdx d2PsidE2 dEdx
//         for (int i = 0; i < 3; ++i) {
//             for (int j = 0; j < 3; ++j) {
//                 *d2edx2 += C(i,j) * dEdx.row(i).transpose() * dEdx.row(j);
//             }
//         }
        
//         // dPsidE d2Edx2
//         const Eigen::Vector3d dPsidE = C * EVoigt;
//         for (int i = 0; i < 3; ++i) {
//             *d2edx2 += dPsidE(i) * d2Edx2[i];
//         }

//         // dedPsi
//         *d2edx2 *= V_;
//     }

// 	return energy;
// }
