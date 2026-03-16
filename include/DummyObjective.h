#ifndef DUMMY_OBJECTIVE_H
#define DUMMY_OBJECTIVE_H

#include <Objective.h>
#include <Material/Shell.h>

class DummyObjective : public Objective
{
public:
    DummyObjective();

    void update(const Eigen::VectorXd &dx) override;
    double energy() const override;
    Eigen::VectorXd gradient() const override;
    Eigen::SparseMatrix<double> hessian() const override;

    std::unique_ptr<Shell> m_pShell;
    Eigen::MatrixXd x, X;
    Eigen::MatrixXd ts;
    Eigen::Matrix3d B;
    double A = 0.5;
};

#endif