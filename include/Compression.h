#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <Eigen/Core>
#include <torch/script.h>

class Tet10;
class Shell;

class Compression
{
public:
    Compression() = delete;

    static Eigen::MatrixX3d run(
        const Tet10& mesh, 
        const Eigen::MatrixX3d& restPos, 
        double E,
        double nu,
        std::ofstream &file,
        bool verbose=false);

    static Eigen::MatrixX3d run(
        const Shell& mesh, 
        const Eigen::MatrixX2d& restPos, 
        const torch::jit::Module& stretchModel,
        const torch::jit::Module& bendModel,
        std::ofstream &file,
        bool verbose=false);

    static std::array<double, 2> fit(const Eigen::MatrixX3d& pos);


    static constexpr int m_sMaxIter = 100;
    static constexpr double m_sTol = 1e-5;
    static constexpr double m_sEps = 1e-5;
};

#endif