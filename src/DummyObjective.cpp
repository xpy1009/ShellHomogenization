#include <DummyObjective.h>
#include <GeometryLib.h>
#include <ElasticityLib.h>
#include <Material/ShellHelper.h>
#include <iostream>

DummyObjective::DummyObjective()
{
    // X.resize(6,3);
    // X << 0,0,0, 1,0,0, 0,1,0, 1,1,0, -1,1,0, 1,-1,0;
    // Eigen::MatrixXd x2d = X.leftCols(2);
	// const Eigen::RowVector2d& _p1 = x2d.row(0);
	// const Eigen::RowVector2d& _p2 = x2d.row(1);
	// const Eigen::RowVector2d& _p3 = x2d.row(2);
    // const Eigen::RowVector2d _v1 = _p3 - _p2;
    // const Eigen::RowVector2d _v2 = _p1 - _p3;
    // const Eigen::RowVector2d _v3 = _p2 - _p1;
    // Eigen::RowVector2d _t1(_v1.y(), -_v1.x()); 
    // Eigen::RowVector2d _t2(_v2.y(), -_v2.x()); 
    // Eigen::RowVector2d _t3(_v3.y(), -_v3.x()); 
    // _t1.normalize();
    // _t2.normalize();
    // _t3.normalize();
    // ts.resize(3,2);
    // ts << _t1, _t2, _t3;
    // X = GeometryLib::cylindricalProject(X, 1, M_PI_4);
    // // std::cout << Geometry2D::shapeOperator(X.topRows(3), X.bottomRows(3), ts) << std::endl;
    // B << 1,6,5,
    //     6,2,4,
    //     5,4,3;
    
    // X.resize(3,3);
    // X << -1,0,0, 1,0,0, 0,1,0;

    Eigen::MatrixXd V(4,3);
    V << 0,0,0, 1,0,0, 0,1,0, 1,1,0;
    Eigen::MatrixXi F(2,3);
    F << 0,1,2, 3,2,1;
    m_pShell = std::make_unique<Shell>(V, F);
    m_pShell->m_v = GeometryLib::cylindricalProject(V, 1, M_PI_4);
    // m_pShell->m_deformed *= 2;
}

void DummyObjective::update(const Eigen::VectorXd &dx)
{
    m_pShell->m_v += dx.reshaped<Eigen::RowMajor>(4,3);
    // X += dx.reshaped<Eigen::RowMajor>(6,3);
    // std::cout << Geometry2D::exteriorDihedralAngle(X.row(0), X.row(1), X.row(2), X.row(3)) << std::endl;
}

double DummyObjective::energy() const
{
    // auto e = Geometry2D::exteriorDihedralAngle(X.row(0), X.row(1), X.row(2), X.row(3));
    // auto e = Geometry2D::height(X.row(0), X.row(1), X.row(2));
    // auto e = Geometry2D::shapeOperator(X.topRows(3), X.bottomRows(3), ts)[0];
    // auto e = Geometry2D::bendingEnergy(X.topRows(3), X.bottomRows(3), ts, A, B);
    auto e = m_pShell->energy();
    return e;
    return 0;
}

Eigen::VectorXd DummyObjective::gradient() const
{
    // Eigen::Matrix<double,12,1> grad;
    // Geometry2D::exteriorDihedralAngle(X.row(0), X.row(1), X.row(2), X.row(3), &grad);
    // Eigen::Matrix<double,9,1> grad;
    // Geometry2D::height(X.row(0), X.row(1), X.row(2), &grad);
    // Eigen::Matrix<double,3,18> tmp;
    // Geometry2D::shapeOperator(X.topRows(3), X.bottomRows(3), ts, &tmp);
    // Eigen::Matrix<double,18,1> grad = tmp.row(0).transpose();
    // Eigen::Matrix<double,18,1> grad;
    // Geometry2D::bendingEnergy(X.topRows(3), X.bottomRows(3), ts, A, B, &grad);
    Eigen::Matrix<double,12,1> grad = m_pShell->gradient();
    return grad;
}

Eigen::SparseMatrix<double> DummyObjective::hessian() const
{
    // return m_pShell->hessian();
    // Eigen::Matrix<double,9,9> hess;
    // Geometry2D::height(X.row(0), X.row(1), X.row(2), nullptr, &hess);
    // Eigen::Matrix<double,12,12> hess;
    // Geometry2D::exteriorDihedralAngle(X.row(0), X.row(1), X.row(2), X.row(3), nullptr, &hess);
    // std::array<Eigen::Matrix<double,18,18>,3> tmp;
    // Geometry2D::shapeOperator(X.topRows(3), X.bottomRows(3), ts, nullptr, &tmp);
    // Eigen::Matrix<double,18,18> hess = tmp[0];
    // Eigen::Matrix<double,18,18> hess;
    // Geometry2D::bendingEnergy(X.topRows(3), X.bottomRows(3), ts, A, B, nullptr, &hess);
    // return hess.sparseView();
    return m_pShell->hessian();
}