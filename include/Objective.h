#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

/**
 * \brief Base class of all objectives.
 *
 * Used for derivative tests
 */
class Objective
{    
public:
    virtual void update(const Eigen::VectorXd &dx) = 0;
    virtual double energy() const = 0;
    virtual Eigen::VectorXd gradient() const = 0;
    virtual Eigen::SparseMatrix<double> hessian() const = 0;
};
#endif