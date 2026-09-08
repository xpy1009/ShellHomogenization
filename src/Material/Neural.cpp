#include <Material/Neural.h>
#include <Element/Shell.h>

double Neural::elasticEnergy(
        const Shell& mesh,
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX2d& X, 
        const torch::jit::Module& stretchModel,
        const torch::jit::Module& bendModel,
        Eigen::VectorXd* gradient,
        Eigen::SparseMatrix<double>* hessian)
{
    double energy = 0.0;

    Eigen::VectorXd stretchGrad;
    Eigen::SparseMatrix<double> stretchHess;
    energy += stretchingEnergy(
        mesh, x, X, stretchModel, gradient!=nullptr ? &stretchGrad : nullptr, hessian!=nullptr ? &stretchHess : nullptr);

    Eigen::VectorXd bendGrad;
    Eigen::SparseMatrix<double> bendHess;
    energy += bendingEnergy(
        mesh, x, X, bendModel, gradient!=nullptr ? &bendGrad : nullptr, hessian!=nullptr ? &bendHess : nullptr);
    
    if (gradient != nullptr) {
        *gradient = stretchGrad + bendGrad;
    }

    if (hessian != nullptr) {
        *hessian = stretchHess + bendHess;
    }

    return energy;
}

double Neural::stretchingEnergy(
    const Shell& mesh,
    const Eigen::MatrixX3d& x, 
    const Eigen::MatrixX2d& X, 
    const torch::jit::Module& model,
    Eigen::VectorXd* gradient,
    Eigen::SparseMatrix<double>* hessian) 
{
    const int nFaces = mesh.nFaces();
    std::vector<Eigen::Matrix<double, 3, 9>> dEdx;
    std::vector<std::array<Eigen::Matrix<double, 9, 9>, 3>> d2Edx2;
    if (gradient != nullptr) {
        gradient->setZero(x.size());
        dEdx.resize(nFaces);
    }
    std::vector<Eigen::Triplet<double>> triplets;
    if (hessian != nullptr) {
        hessian->resize(x.size(), x.size());
        triplets.reserve(nFaces * 81);
        d2Edx2.resize(nFaces);
    }


    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> strains(nFaces, 3);
    Eigen::RowVectorXd As = mesh.areas(X).transpose();

    for (int i = 0; i < nFaces; ++i) {
        const Eigen::Vector3d E = mesh.GreenStrainVoigt(
            i, X, x, (gradient!=nullptr || hessian!=nullptr) ? &dEdx[i] : nullptr, hessian!=nullptr ? &d2Edx2[i] : nullptr);

        strains.row(i) = E.transpose();
    }

    const torch::TensorOptions options(torch::kFloat64);
    const torch::Tensor epsilon = torch::from_blob(strains.data(), {nFaces, 3}, options);
    const std::vector<torch::jit::IValue> epsilon_at = {epsilon};

    const at::Tensor psi = model.get_method("energy")(epsilon_at).toTensor();
    const at::Tensor a = at::from_blob(As.data(), {nFaces, 1}, options);
    const double energy = (psi * a).sum().item<double>();


    if (gradient != nullptr || hessian != nullptr) {
        const at::Tensor dedE_at = model.get_method("gradient")(epsilon_at).toTensor();
        const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> dedE_eigen(dedE_at.data_ptr<double>(), dedE_at.size(0), 3);

        if (gradient != nullptr) {
            for(int i = 0; i < nFaces; ++i) {
                const Eigen::Matrix<double, 9, 1> gi = As[i] * (dedE_eigen.row(i) * dEdx[i]).transpose();
                const Eigen::Array3i ind = mesh.indices(i);
                for (int j = 0; j < 3; ++j) {
                    gradient->segment(3 * ind(j), 3) += gi.segment(3 * j, 3);
                }
            }
        } 

        if (hessian != nullptr) {
            const at::Tensor d2edE2_at = model.get_method("hessian")(epsilon_at).toTensor();
            for(int i = 0; i < nFaces; ++i) {
                const Eigen::Array3i ind = mesh.indices(i);
                const Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> d2edE2(d2edE2_at[i].data_ptr<double>(), 3, 3);
                Eigen::Matrix<double, 9, 9> hi = dEdx[i].transpose() * d2edE2 * dEdx[i];
                for (int j = 0; j < 3; ++j) {
                    hi += d2Edx2[i][j] * dedE_eigen(i, j);
                }
                hi *= As[i];

                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int d = 0; d < 3; ++d) {
                            for (int dd = 0; dd < 3; ++dd) {
                                triplets.emplace_back(3 * ind(j) + d, 3 * ind(k) + dd, hi(3 * j + d, 3 * k + dd));                           
                            }
                        }
                    }
                }
            }
            hessian->setFromTriplets(triplets.begin(), triplets.end());
        }
    }

    return energy;
}

double Neural::bendingEnergy(
    const Shell& mesh,
    const Eigen::MatrixX3d& x, 
    const Eigen::MatrixX2d& X, 
    const torch::jit::Module& model,
    Eigen::VectorXd* gradient,
    Eigen::SparseMatrix<double>* hessian) 
{
    const int nFaces = mesh.nFaces();
    std::vector<Eigen::Matrix<double, 3, 18>> dkdx;
    std::vector<std::array<Eigen::Matrix<double, 18, 18>, 3>> d2kdx2;
    if (gradient != nullptr) {
        gradient->setZero(x.size());
        dkdx.resize(nFaces);
    }
    std::vector<Eigen::Triplet<double>> triplets;
    if (hessian != nullptr) {
        hessian->resize(x.size(), x.size());
        triplets.reserve(nFaces * 324);
        d2kdx2.resize(nFaces);
    }


    Eigen::RowVectorXd As = mesh.areas(X).transpose(); 
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> ks(nFaces, 3);
    for (int i = 0; i < nFaces; ++i) {
        ks.row(i) = mesh.shapeOperator(i, X, x, (gradient!=nullptr || hessian!=nullptr) ? &dkdx[i] : nullptr, hessian!=nullptr ? &d2kdx2[i] : nullptr).transpose();
    }


    const torch::TensorOptions options(torch::kFloat64);
    const torch::Tensor epsilon = torch::from_blob(ks.data(), {nFaces, 3}, options);
    const std::vector<torch::jit::IValue> epsilon_at = {epsilon};

    const at::Tensor psi = model.get_method("energy")(epsilon_at).toTensor();
    const at::Tensor a = at::from_blob(As.data(), {nFaces, 1}, options);
    const double energy = (psi * a).sum().item<double>();


    if (gradient != nullptr || hessian != nullptr) {
        const at::Tensor dedk_at = model.get_method("gradient")(epsilon_at).toTensor();
        const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> dedk_eigen(dedk_at.data_ptr<double>(), nFaces, 3);

        if (gradient != nullptr) {
            for(int i = 0; i < nFaces; ++i) {
                const Eigen::Matrix<double, 18, 1> gi = As[i] * (dedk_eigen.row(i) * dkdx[i]).transpose();
                const Eigen::Array3i ind = mesh.indices(i);
                const std::array<int, 3> oppositeInd = mesh.oppositeIndices(i);
                for (int j = 0; j < 3; ++j) {
                    gradient->segment(3 * ind(j), 3) += gi.segment(3 * j, 3);
                    if (oppositeInd[j] > -1) {
                        gradient->segment(3 * oppositeInd[j], 3) += gi.segment(9 + 3 * j, 3);
                    }
                }
            }
        } 


        if (hessian != nullptr) {
            const at::Tensor d2edk2_at = model.get_method("hessian")(epsilon_at).toTensor();
            for(int i = 0; i < nFaces; ++i) {
                const Eigen::Array3i ind = mesh.indices(i);
                std::array<int, 3> oppositeInd = mesh.oppositeIndices(i);
                for (int &oid : oppositeInd) {
                    oid = std::max(oid, 0);
                }

                const Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> d2edk2(d2edk2_at[i].data_ptr<double>(), 3, 3);
                Eigen::Matrix<double, 18, 18> hi = dkdx[i].transpose() * d2edk2 * dkdx[i];
                for (int j = 0; j < 3; ++j) {
                    hi += d2kdx2[i][j] * dedk_eigen(i, j);
                }
                hi *= As[i];

                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int d = 0; d < 3; ++d) {
                            for (int dd = 0; dd < 3; ++dd) {
                                triplets.emplace_back(3 * ind(j) + d, 3 * ind(k) + dd, hi(3 * j + d, 3 * k + dd)); 
                                triplets.emplace_back(3 * ind(j) + d, 3 * oppositeInd[k] + dd, hi(3 * j + d, 9 + 3 * k + dd));                           
                                triplets.emplace_back(3 * oppositeInd[j] + d, 3 * ind(k) + dd, hi(9 + 3 * j + d, 3 * k + dd));                           
                                triplets.emplace_back(3 * oppositeInd[j] + d, 3 * oppositeInd[k] + dd, hi(9 + 3 * j + d, 9 + 3 * k + dd));                           
                            }
                        }
                    }
                }
            }

            hessian->setFromTriplets(triplets.begin(), triplets.end());
        }
    }
    
    return energy;
}

