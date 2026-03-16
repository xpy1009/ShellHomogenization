#ifndef SHELL_MESH_H
#define SHELL_MESH_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

class ShellMesh
{
public:
    ShellMesh() = delete;
    static void getVF(Eigen::MatrixX3d &vertices, Eigen::MatrixX3i &faces);
    static void biMatShell(Eigen::MatrixX3d &vertices, Eigen::MatrixX3i &faces);
    static void biMatShellPeriodic(Eigen::MatrixX3i &faces, 
                                Eigen::SparseMatrix<double>& paramsToCoords, 
                                Eigen::VectorXd& params, 
                                Eigen::Matrix2d& bbox, 
                                std::vector<int>* bdyIdx1=nullptr, 
                                std::vector<int>* bdyIdx2=nullptr);
};

#endif