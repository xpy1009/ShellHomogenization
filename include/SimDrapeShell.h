// #ifndef SIM_DRAPE_SHELL_H
// #define SIM_DRAPE_SHELL_H

// #include <Objective.h>
// #include <DiscreteShell.h>

// class SimDrapeShell : public Objective
// {
// public:
//     SimDrapeShell();
//     void visualize();
    
//     void update(const Eigen::VectorXd &dx) override;
//     double energy() const override;
//     Eigen::VectorXd gradient() const override;
//     Eigen::SparseMatrix<double> hessian() const override;

//     bool step(int iter);
//     void sim();

//     std::unique_ptr<DiscreteShell> m_shell;
//     // Eigen::SparseMatrix<double> m_massMatrix;
//     // double rho = 1e6;

//     // double tol = 1e-5; // tolerance
//     // int max_iter = 500; // max iteration
//     // std::vector<int> dirichlet_idx;
// };

// #endif