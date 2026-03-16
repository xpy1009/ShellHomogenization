#ifndef SIM_SHELL_H
#define SIM_SHELL_H

#include <Objective.h>
#include <Material/Shell.h>

class SimShell : public Objective
{
public:
    SimShell();
    void visualize();
    
    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;

    // bool step(int iter);

    std::unique_ptr<Shell> m_shell;

    double m_tol = 1e-5; // tolerance
    int m_maxIter = 500; // max iteration
    std::vector<int> m_dirichletIdx;
};

#endif