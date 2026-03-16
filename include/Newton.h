#ifndef NEWTON_H
#define NEWTON_H

#include <Objective.h>

class Newton
{
public:
    Newton() = delete;
    static Eigen::VectorXd solve(Eigen::SparseMatrix<double> &A, const Eigen::VectorXd& b);

    // https://en.wikipedia.org/wiki/Woodbury_matrix_identity
    // U = V, C = I
    static Eigen::VectorXd solve(Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &U, const Eigen::VectorXd& b);

    static void solveOpt(Objective& obj, double maxIter, double tol);
};


#endif