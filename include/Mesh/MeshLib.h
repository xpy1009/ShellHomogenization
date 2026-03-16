#ifndef MESH_LIB_H
#define MESH_LIB_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class MeshLib
{    
public:
    MeshLib() = delete;
    static void getQuadTets(Eigen::MatrixXd &vertices, Eigen::MatrixXi &tets);
    static void generateBox(Eigen::MatrixXd &vertices, Eigen::MatrixXi &tets);
    static void rectangle(Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces);

    // From tactile/demo/psdemo.cpp
    // Return a list of 2D points for Clipper2
    static void constructTiling(int idx, Eigen::Matrix2d& bbox, std::vector<std::vector<double>>& tiles);

    // Inflate path
    static std::vector<std::vector<double>> inflatePaths(const std::vector<std::vector<double>>& lists);

    static void structuredSheet(Eigen::MatrixXd &vertices, Eigen::MatrixXi &tets);

    // Build periodic mesh for homogenization
    static void structuredSheetPeriodic(Eigen::MatrixXi &tets, Eigen::SparseMatrix<double>& paramsToCoords, Eigen::VectorXd& params, Eigen::Matrix2d& bbox, std::vector<int>* bdyIdx1=nullptr, std::vector<int>* bdyIdx2=nullptr);
    static std::vector<std::pair<int,int>> setPeriodic(const std::vector<std::pair<int, int>>& dimTags, double dx, double dy);
    static void buildReduced(const Eigen::MatrixXd& vertices, 
        const std::vector<std::pair<int, int>>& dimTags1, 
        const std::vector<std::pair<int, int>>& dimTags2, 
        Eigen::SparseMatrix<double>& reducedToFull, Eigen::VectorXd& reduced, 
        std::vector<int>* bdyIdx1, std::vector<int>* bdyIdx2);

    // analysis
    static void structuredSheetMidsurface(Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces);
};

#endif