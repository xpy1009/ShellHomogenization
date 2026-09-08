#ifndef MESH_LIB_H
#define MESH_LIB_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class MeshLib
{    
public:
    MeshLib() = delete;
    
    using Tet10 = Eigen::Matrix<int, Eigen::Dynamic, 10>;

    static void getQuadTets(Eigen::MatrixX3d &vertices, Tet10 &tets);

    static void rectangle(double dx, double dy, Eigen::MatrixX2d &vertices, Eigen::MatrixX3i &faces, bool preview=false);

    // multiple tiles for native compression
    static void structuredSheet(
        int IH, 
        const std::vector<double>& params,
        double inflateWidth,
        const std::array<int, 2>& size,
        double height,
        Eigen::MatrixX3d &vertices, 
        Tet10 &tets,
        bool preview=false);

    // periodic tet mesh for homogenization
    static void structuredSheet(
        int IH,
        const std::vector<double>& params,
        double inflateWidth,
        double height,
        Tet10 &tets, 
        Eigen::SparseMatrix<double>& proj, 
        Eigen::VectorXd& reduced,
        std::array<double, 4>& trans,
        bool preview=false);
    
    static std::vector<std::pair<int,int>> setPeriodic(const std::vector<std::pair<double, double>>& offset);

    static void getReduced(
        const std::vector<std::pair<int, int>>& sDimTags, 
        const Eigen::MatrixXd& vertices,  
        Eigen::SparseMatrix<double>& proj, 
        Eigen::VectorXd& reduced);
};

#endif