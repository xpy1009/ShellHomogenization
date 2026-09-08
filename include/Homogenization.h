#ifndef HOMOGENIZATION_H
#define HOMOGENIZATION_H

#include <Element/Tet10.h>

class Homogenization
{
public:
    Homogenization(
        const Tet10& mesh, 
        const Eigen::SparseMatrix<double>& proj, 
        const Eigen::VectorXd& var, 
        const std::array<double, 4>& trans,
        double E,
        double nu);

    // set stretch direction
    void setDir(double angle);
    // set weight for constraints
    void setWeight(double weight);
    // solve while satisfying constraints
    void solve(Eigen::VectorXd& var, bool verbose=false);

    void stretch(std::ofstream &file, bool verbose=false);
    void saveStretchData(std::ofstream &file, const Eigen::VectorXd& var, double angle) const;

    void bend(std::ofstream &file, bool verbose=false);
    void saveBendData(std::ofstream &file, const Eigen::VectorXd& var) const;


    // apply soft constraints on center positions via Woodbury matrix identity
    static Eigen::SparseMatrix<double> WoodburyMatrix(const Eigen::SparseMatrix<double>& massMatrix);

    
    const Tet10& m_mesh;
    const Eigen::SparseMatrix<double>& m_proj;
    Eigen::VectorXd m_restVar;

    double m_area;
    Eigen::MatrixX3d m_restPos;
    Eigen::SparseMatrix<double> m_massMat;

    // projection matrix with varying stretch directions
    Eigen::SparseMatrix<double> m_rProj;

    // lame parameters
    double m_lambda, m_mu;

    // for bending
    double m_curvature = 0.0, m_angle = 0.0;

    // for constraint
    double m_weight;
    Eigen::SparseMatrix<double> m_U, m_Uproj;
    static constexpr double m_sConstraintTol = 1e-6;
    static constexpr double m_sInitWeight = 1e4;
    static constexpr double m_sMaxWeight = 1e12;

    // for optimization
    std::vector<int> m_fixedVar;
    std::function<double(const Eigen::VectorXd &, Eigen::VectorXd *, Eigen::SparseMatrix<double> *)> m_objFunc;
    static constexpr int m_sMaxIter = 50;
    static constexpr double m_sTol = 1e-5;
};

#endif