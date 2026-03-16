#ifndef SIM_TWIST_H
#define SIM_TWIST_H

#include <Objective.h>
#include <Tetrahedron.h>

class SimTwist : public Objective
{
public:
    SimTwist();
    void visualize();
    
    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;

    bool step(int iter);
    void sim();

    std::unique_ptr<Tetrahedron> tets;

    double tol = 1e-5; // tolerance
    int max_iter = 500; // max iteration
    std::vector<int> dirichlet_idx;
};

#endif