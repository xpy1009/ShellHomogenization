#ifndef HOMOGENIZE_STRETCHING_H
#define HOMOGENIZE_STRETCHING_H

#include <Objective.h>
#include <Tetrahedron.h>

class HomogenizeStretching : public Objective
{
public:
    HomogenizeStretching();
    void visualize();
    
    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;
    Eigen::SparseMatrix<double> woodburyHessian() const;

    bool step(int iter);
    bool sim();
    void generateData();
    void updateMap();

    std::unique_ptr<Tetrahedron> m_tets;

    Eigen::VectorXd m_params, m_paramsUndeformed;
    // Eigen::MatrixXd m_coords;
    Eigen::SparseMatrix<double> m_paramsToCoords;
    Eigen::Matrix2d m_bbox;
    std::vector<int> m_bdyIdx1, m_bdyIdx2;

    double m_weight = 1e6;
    Eigen::SparseMatrix<double> m_massMatrix;
    Eigen::SparseMatrix<double> m_woodU;

    const int m_dim;
    double m_angle = 0.0;

    double m_tol = 1e-5; // tolerance
    int m_maxIter = 0; // max iteration
    std::vector<int> m_dirichletIdx;
};

#endif