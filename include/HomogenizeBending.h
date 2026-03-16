#ifndef HOMOGENIZE_BENDING_H
#define HOMOGENIZE_BENDING_H

#include <Objective.h>
#include <Tetrahedron.h>

class HomogenizeBending : public Objective
{
public:
    HomogenizeBending();
    void visualize();
    
    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;
    Eigen::SparseMatrix<double> woodburyHessian() const;

    bool step(int iter);
    bool sim();
    void generateData();
    void generateDataHYLC();

    std::unique_ptr<Tetrahedron> m_tets;

    Eigen::VectorXd m_params, m_paramsUndeformed;
    Eigen::MatrixXd m_coords;
    Eigen::SparseMatrix<double> m_paramsToCoords;
    Eigen::Matrix2d m_bbox;

    double m_weight = 1e6;
    Eigen::SparseMatrix<double> m_massMatrix;
    Eigen::SparseMatrix<double> m_woodU;

    double m_k = 0.0;
    double m_angle = 0.0;

    double m_tol = 1e-5; // tolerance
    int m_maxIter = 50; // max iteration
    std::vector<int> m_dirichletIdx;
};

#endif