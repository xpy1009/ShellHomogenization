#include <HomogenizeBending.h>
#include <Mesh/MeshLib.h>
#include <Newton.h>
#include <GeometryLib.h>

#include <Eigen/Dense>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <fstream>

HomogenizeBending::HomogenizeBending()
{
    Eigen::MatrixXi T;
    MeshLib::structuredSheetPeriodic(T, m_paramsToCoords, m_params, m_bbox);

    // scale such that one tile is 1 cm^2
    m_params *= 1e-2;
    m_paramsToCoords.rightCols(3) *= 1e-2;
    m_bbox *= 1e-2;
    m_params.tail(3) << 1,0,1;

    const Eigen::VectorXd coords = m_paramsToCoords * m_params;
    const Eigen::MatrixXd V = coords.reshaped<Eigen::RowMajor>(coords.size()/3, 3);
    m_tets = std::make_unique<Tetrahedron>(V, T);

    m_massMatrix = m_tets->massMatrix();

    // Woodbury matrix U
    const int n = V.rows();
    Eigen::VectorXd e(n); e.setOnes();
    e = m_massMatrix.transpose() * e;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * n);
    for (int i = 0; i < n; i++) {
        for (int d = 0; d < 3; d++) {
            triplets.emplace_back(3*i+d, d, e(i));
        }
    }
    m_woodU.resize(3 * n, 3);
    m_woodU.setFromTriplets(triplets.begin(), triplets.end());
    m_woodU = m_paramsToCoords.transpose() * m_woodU;

    Eigen::VectorXd dx(m_params.size()); dx.setZero();
    update(dx);
    sim();
    m_paramsUndeformed = m_params;
}

void HomogenizeBending::visualize()
{    
    // Initialize Polyscope
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
    polyscope::init();

    // Register the mesh with Polyscope
    polyscope::SurfaceMesh* mesh = polyscope::registerSurfaceMesh("input mesh", m_tets->deformed, m_tets->faces);
    mesh->setEdgeWidth(1.0);
    mesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    mesh->setBackFaceColor({0.2, 0.2, 0.2});

    // Eigen::MatrixXd vd(dirichlet_idx.size()/3, 3);
    // for (int i = 0; i < dirichlet_idx.size()/3; i++)
    // {
    //     int vid = dirichlet_idx[3*i]/3;
    //     vd.row(i) = tets->deformed.row(vid);
    // }
    // polyscope::PointCloud* psCloud = polyscope::registerPointCloud("origin", vd);

    polyscope::state::userCallback = [&]() 
    {
        ImGui::PushItemWidth(100);
        if (ImGui::Button("Step")) {
            step(0);
            mesh->updateVertexPositions(m_tets->deformed);
        }
    };

    // Show the GUI
    polyscope::show();
}

void HomogenizeBending::update(const Eigen::VectorXd &dx)
{
    if (dx.size() != m_params.size()) {
        throw std::runtime_error("# params neq # dx?");
    }

    m_params += dx;
    const Eigen::VectorXd coords = m_paramsToCoords * m_params;
    m_coords = coords.reshaped<Eigen::RowMajor>(coords.size()/3, 3);
    m_tets->deformed = GeometryLib::cylindricalProject(m_coords, m_k, m_angle);
}

double HomogenizeBending::energy() const
{
    // elasticity
    double e = m_tets->energy();

    // constraint
    const Eigen::RowVector3d cs = (m_massMatrix * m_coords).colwise().sum();
    e += 0.5 * m_weight * cs.squaredNorm();

    return e;
}

Eigen::VectorXd HomogenizeBending::gradient() const
{
    // elasticity
    Eigen::VectorXd gTets = m_tets->gradient();
    // map to coords before cylindrical projection
    const Eigen::SparseMatrix<double> dxdw = GeometryLib::cylindricalJacobian(m_coords, m_k, m_angle);
    gTets = dxdw.transpose() * gTets;

    // constraint
    const Eigen::RowVector3d cs = (m_massMatrix * m_coords).colwise().sum();
    Eigen::MatrixXd e(m_coords.rows(), m_coords.cols()); e.rowwise() = cs;
    const Eigen::VectorXd gConstraints = m_weight * (m_massMatrix * e).reshaped<Eigen::RowMajor>();

    // full
    const Eigen::VectorXd gFull = gConstraints + gTets;
    // full to reduced
    Eigen::VectorXd grad = m_paramsToCoords.transpose() * gFull;

    // fixed dofs
    for(const int& idx : m_dirichletIdx) {
        grad(idx) = 0.0;
    }

    return grad;
}

Eigen::SparseMatrix<double> HomogenizeBending::hessian() const
{
    Eigen::SparseMatrix<double> hess = woodburyHessian();

    hess += m_weight * m_woodU * m_woodU.transpose();

    return hess;
}

Eigen::SparseMatrix<double> HomogenizeBending::woodburyHessian() const
{
    // elasticity
    Eigen::SparseMatrix<double> hTets = m_tets->hessian();
    // map to coords before cylindrical projection
    const Eigen::VectorXd gTets = m_tets->gradient();
    const Eigen::SparseMatrix<double> dxdw = GeometryLib::cylindricalJacobian(m_coords, m_k, m_angle);
    const Eigen::SparseMatrix<double> d2xdw2 = GeometryLib::cylindricalHessian(m_coords, m_k, m_angle, gTets);
    hTets = dxdw.transpose() * hTets * dxdw + d2xdw2;

    // full to reduced
    Eigen::SparseMatrix<double> hess = m_paramsToCoords.transpose() * hTets * m_paramsToCoords;
   
    for (const int& idx : m_dirichletIdx) {
        hess.row(idx) *= 0.0;
        hess.col(idx) *= 0.0;
        hess.coeffRef(idx,idx) = 1.0;
    }
    return hess;
}

bool HomogenizeBending::step(int iter)
{
    const Eigen::VectorXd grad = gradient();
    
    const double grad_norm = grad.norm();
    std::cout << "[NEWTON] iter " << iter << "/" << m_maxIter << " gradient norm: " << grad_norm << " tol: " << m_tol << std::endl;

    if (iter == m_maxIter) {
        std::cout << m_params.tail(3).transpose() << " " << m_angle * M_1_PI << std::endl;
        // visualize();
        return true;
        // throw std::runtime_error("Exceed max iteration?!!");
    }
    if (grad_norm < m_tol) {
        // check if constraints are satisfied
        const Eigen::RowVector3d cs = (m_massMatrix * m_coords).colwise().sum().cwiseAbs();
        if (cs.maxCoeff() > 1e-10 && m_weight < 1e12) {
            std::cout << "Constraints: " << cs << " not statisfied, increase weight" << std::endl;
            m_weight *= 1e2;
            return false;
        }
        if (cs.maxCoeff() > 1e-8 && m_weight == 1e12) {
            std::cout << m_k << " " << m_angle << " " << cs << std::endl;
            std::exit(0);
        }
        return true;
    }

    // linear solve
    Eigen::SparseMatrix<double> K = woodburyHessian();
    Eigen::SparseMatrix<double> U = sqrt(m_weight) * m_woodU;
    for (const auto& i : m_dirichletIdx) {
        U.row(i) *= 0.0;
    }
    const Eigen::VectorXd dir = Newton::solve(K, U, -grad);

    // line search
    const double E0 = energy();
    double step_size = 1;
    int cnt = 0;

    while (true) {
        Eigen::VectorXd dx = step_size * dir;
        update(dx);

        const double E1 = energy();

        if (E1 - E0 < 0 || cnt > 15) {
            if (cnt > 15) {
                std::cout << "line search max" << std::endl;
            }
            break;
        }

        dx *= -1;
        update(dx);
        step_size *= 0.5;
        cnt++;
    }

    return false;
}

bool HomogenizeBending::sim()
{
    int iter = 0;
    while (!step(iter++));
    return iter >= m_maxIter ? false : true;
}

void HomogenizeBending::generateData()
{
    Eigen::VectorXd dx(m_params.size()); dx.setZero();
    std::ofstream file("bending.txt");
    auto saveData = [&]() {
        Eigen::Matrix2d R; R << cos(m_angle), -sin(m_angle), sin(m_angle), cos(m_angle);
        Eigen::Matrix2d S; S << m_k, 0, 0, 0;
        S = R * S * R.transpose();
        const Eigen::Vector3d kVoigt(S(0,0), S(1,1), 2*S(0,1));
        const double area = abs(m_bbox.determinant());
        const Eigen::Vector3d mVoigt = m_tets->moment() / area;
        const double psi = m_tets->energy() / area;
        file << m_k << " "  << m_angle << " "
            << kVoigt(0) << " " << kVoigt(1) << " " << kVoigt(2) << " " 
            << mVoigt(0) << " " << mVoigt(1) << " " << mVoigt(2) << " "
            << psi << std::endl;
    };

    const int an = 20;
    const int kn = 20;
    const double kStep = 10;
    for (int ai = 0; ai < an; ++ai) {
        m_angle = ai * M_PI / an;

        m_params = m_paramsUndeformed;
        m_weight = 1e6;
        for (int ki = 0; ki < kn; ++ki) {
            m_k = kStep * (ki + 1);
            update(dx);
            if (sim())
                saveData();
        }

        m_params = m_paramsUndeformed;
        m_weight = 1e6;
        for (int ki = 0; ki < kn; ++ki) {
            m_k = -kStep * (ki + 1);
            update(dx);
            if (sim())
                saveData();
        }
    }
    file.close();

}

void HomogenizeBending::generateDataHYLC()
{
    const double area = abs(m_bbox.determinant());
    const int pbcStart = m_params.size() - 3;
    m_dirichletIdx = {pbcStart, pbcStart + 1, pbcStart + 2};
    Eigen::VectorXd dx(m_params.size()); dx.setZero();

    constexpr int kn = 100;
    constexpr double kStep = 2;
    constexpr double eps = 1e-3;

    auto saveData = [&](std::ofstream& file, double step) {
        m_params = m_paramsUndeformed;
        m_weight = 1e6;
        for (int ki = 0; ki < kn; ++ki) {
            m_k = step * (ki + 1);
            update(dx);
            sim();
            double e = m_tets->energy() / area;

            m_k += eps;
            update(dx);
            sim();
            double e_ = m_tets->energy() / area;
            double deds = (e_ - e) / eps;
            file << step * (ki + 1) << " " << e << " " << deds << std::endl;
        }
    };

    std::ofstream fs3("s3.txt");
    m_angle = 0.0;
    saveData(fs3, kStep);
    saveData(fs3, -kStep);
    fs3.close();

    std::ofstream fs4("s4.txt");
    m_angle = M_PI_2;
    saveData(fs4, kStep);
    saveData(fs4, -kStep);
    fs4.close();

}
