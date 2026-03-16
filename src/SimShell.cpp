#include <SimShell.h>
#include <Mesh/MeshLib.h>
#include <GeometryLib.h>
#include <Newton.h>
#include <Mesh/ShellMesh.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

SimShell::SimShell()
{
    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    Eigen::Matrix2d bbox;
    Eigen::SparseMatrix<double> paramsToCoords;
    Eigen::VectorXd params;
    std::vector<int> bdyIdx1, bdyIdx2;
    MeshLib::rectangle(V, F);
    // ShellMesh::biMatShellPeriodic(F, paramsToCoords, params, bbox, &bdyIdx1, &bdyIdx2);
    m_shell = std::make_unique<Shell>(V, F);

    // double xmin = V.col(0).minCoeff();
    // double xmax = V.col(0).maxCoeff();
    // double ymin = V.col(1).minCoeff();
    // double ymax = V.col(1).maxCoeff();

    // for (int i = 0; i < V.rows(); i++)
    // {
    //     if ((abs(V(i, 0) - xmin) < 1e-5 || abs(V(i, 0) - xmax) < 1e-5) && 
    //         (abs(V(i, 1) - ymin) < 1e-5 || abs(V(i, 1) - ymax) < 1e-5)) {
    //         m_dirichletIdx.push_back(3*i);
    //         m_dirichletIdx.push_back(3*i+1);
    //         m_dirichletIdx.push_back(3*i+2);
    //     }
    // }
    // m_shell->m_v = GeometryLib::cylindricalProject(m_shell->m_v, 2, M_PI_4);
}

void SimShell::visualize()
{
    // Initialize Polyscope
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
    polyscope::init();

    // Register the mesh with Polyscope
    polyscope::SurfaceMesh* mesh = polyscope::registerSurfaceMesh("input mesh", m_shell->m_v, m_shell->m_F);
    mesh->setEdgeWidth(1.0);
    mesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    mesh->setBackFaceColor({0.2, 0.2, 0.2});

    Eigen::MatrixX3d points(m_dirichletIdx.size()/3, 3);
    for (int i = 0; i < m_dirichletIdx.size(); i+=3)
    {
        points.row(i/3) = m_shell->m_v.row(m_dirichletIdx[i]/3);
    }
    polyscope::registerPointCloud("points", points);

    polyscope::state::userCallback = [&]() 
    {
        ImGui::PushItemWidth(100);
        if (ImGui::Button("Step")) {
            Newton::solveOpt(*this, 1, m_tol);
            mesh->updateVertexPositions(m_shell->m_v);
        }
    };

    // Show the GUI
    polyscope::show();
}

void SimShell::update(const Eigen::VectorXd &dx)
{
    m_shell->m_v += dx.reshaped<Eigen::RowMajor>(dx.size()/3, 3);
}

double SimShell::energy() const
{
    double e = m_shell->energy();
    return e;
}

Eigen::VectorXd SimShell::gradient() const 
{
    Eigen::VectorXd g = m_shell->gradient();
    g(m_dirichletIdx).setZero();
    return g;
}

Eigen::SparseMatrix<double> SimShell::hessian() const
{
    Eigen::SparseMatrix<double> h = m_shell->hessian();
    for (const int& idx : m_dirichletIdx) {
        h.row(idx) *= 0.0;
        h.col(idx) *= 0.0;
        h.coeffRef(idx,idx) = 1.0;
    }
    return h;
}
