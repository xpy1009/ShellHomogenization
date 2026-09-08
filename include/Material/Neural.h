#ifndef MATERIAL_NEURAL_H
#define MATERIAL_NEURAL_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <torch/script.h>

class Shell;

class Neural
{
public:
    Neural() = delete;

    static double elasticEnergy(
        const Shell& mesh,
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX2d& X, 
        const torch::jit::Module& stretchModel,
        const torch::jit::Module& bendModel,
        Eigen::VectorXd* gradient=nullptr,
        Eigen::SparseMatrix<double>* hessian=nullptr);

    static double stretchingEnergy(
        const Shell& mesh,
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX2d& X, 
        const torch::jit::Module& model,
        Eigen::VectorXd* gradient=nullptr,
        Eigen::SparseMatrix<double>* hessian=nullptr);

    static double bendingEnergy(
        const Shell& mesh,
        const Eigen::MatrixX3d& x, 
        const Eigen::MatrixX2d& X, 
        const torch::jit::Module& model,
        Eigen::VectorXd* gradient=nullptr,
        Eigen::SparseMatrix<double>* hessian=nullptr);
};

#endif