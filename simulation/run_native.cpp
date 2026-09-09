#include <iostream>
#include <fstream>

#include <Mesh/MeshLib.h>
#include <Element/Tet10.h>
#include <Homogenization.h>
#include <Compression.h>

#include <igl/writeOBJ.h>


int main(int argc, char **argv)
{
    const auto start = std::chrono::high_resolution_clock::now();

    // tiling pattern
    constexpr int IH = 16;
    const std::vector<double> params = {0.11};

    // material parameters
    constexpr double height = 4e-3;
    constexpr double inflateWidth = 1e-3;
    // elasticity
    constexpr double E = 1e6;
    constexpr double nu = 0.4;


    constexpr bool verbose = false;
    // homogenization
    {
        MeshLib::Tet10 T;
        Eigen::SparseMatrix<double> proj;
        Eigen::VectorXd var;
        std::array<double, 4> trans;
        MeshLib::structuredSheet(IH, params, inflateWidth, height, T, proj, var, trans);
        Eigen::MatrixX3d restPos = (proj * var).reshaped<Eigen::RowMajor>(var.size()/3, 3);
        Tet10 mesh(T);

        Homogenization sim(mesh, proj, var, trans, E, nu);
        std::ofstream sFile("../data/stretching.txt");
        sim.stretch(sFile, verbose);
        std::ofstream bFile("../data/bending.txt");
        sim.bend(bFile, verbose);
    }


    // compression
    for (int i = 0; i < 3; ++i) {
        const int n = 2 * i + 3;
        const std::array<int, 2> size = {10, n};
        MeshLib::Tet10 T;
        Eigen::MatrixX3d restPos;
        MeshLib::structuredSheet(IH, params, inflateWidth, size, height, restPos, T);
        Tet10 mesh(T);

        std::cout << "Compressing 9 x " + std::to_string(n) + " tiles in "
                  << restPos.col(0).maxCoeff() - restPos.col(0).minCoeff() << "m x " 
                  << restPos.col(1).maxCoeff() - restPos.col(1).minCoeff() << "m..." << std::endl;

        std::ofstream cFile("../data/N" + std::to_string(n) + ".txt");
        const Eigen::MatrixX3d pos = Compression::run(mesh, restPos, E, nu, cFile, verbose);
        igl::writeOBJ("../data/N" + std::to_string(n) + ".obj", pos, mesh.m_F);
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> diff = end - start;
    std::cout << diff << std::endl;
    
    return 0;
}