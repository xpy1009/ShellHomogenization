#ifndef SIM_DRAPE_H
#define SIM_DRAPE_H

#include <Objective.h>
#include <Tetrahedron.h>

class SimDrape : public Objective
{
public:
    SimDrape();
    void visualize();
    
    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;

    bool step(int iter);
    void sim();

    void guess();

    std::unique_ptr<Tetrahedron> m_tets;
    Eigen::SparseMatrix<double> m_massMatrix;
    double rho_g = 4e3 * 9.8;

    double tol = 1e-5; // tolerance
    int max_iter = 500; // max iteration
    std::vector<int> dirichlet_idx;
};

#endif