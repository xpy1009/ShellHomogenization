#ifndef SIM_PRESS_H
#define SIM_PRESS_H

#include <Objective.h>
#include <Tetrahedron.h>

class SimPress : public Objective
{
public:
    SimPress();
    void visualize();
    void guess();
    
    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;

    bool step(int iter);
    void sim();

    std::unique_ptr<Tetrahedron> m_tets;

    double tol = 1e-5; // tolerance
    int max_iter = 500; // max iteration
    std::vector<int> dirichlet_idx;
};

#endif