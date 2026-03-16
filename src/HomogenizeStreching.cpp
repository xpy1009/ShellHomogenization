#include <HomogenizeStretching.h>
#include <Mesh/MeshLib.h>
#include <Newton.h>
#include <GeometryLib.h>

#include <Eigen/Dense>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <fstream>

HomogenizeStretching::HomogenizeStretching(): m_dim(2)
{
    Eigen::MatrixXi T;
    MeshLib::structuredSheetPeriodic(T, m_paramsToCoords, m_params, m_bbox, &m_bdyIdx1, &m_bdyIdx2);

    // scale such that one tile is 1 cm^2
    m_params *= 1e-2;
    m_paramsToCoords.rightCols(3) *= 1e-2;
    m_bbox *= 1e-2;
    m_params.tail(3) << 1,0,1;

    const Eigen::VectorXd coords = m_paramsToCoords * m_params;
    const Eigen::MatrixXd V = coords.reshaped<Eigen::RowMajor>(coords.size()/3, 3);
    m_tets = std::make_unique<Tetrahedron>(V, T);
    m_tets->m_model = Tetrahedron::Model::NeoHookean;

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

    m_maxIter = 50;
    sim();

    // forbid buckle
    if (m_dim == 2) {
        int nv = (m_params.size() - 3) / 3;
        Eigen::VectorXd params(2 * nv + 4);
        params.head(2 * nv) = m_params.head(3 * nv).reshaped(3, nv).topRows(2).reshaped();
        params.tail(4) << 1, 1,0,1;
        std::vector<Eigen::Triplet<double>> tri;
        int zid = params.size() - 4;
        for (int i = 0; i < nv; ++i)
        {
            tri.push_back(Eigen::Triplet<double>(3*i, 2*i, 1));
            tri.push_back(Eigen::Triplet<double>(3*i+1, 2*i+1, 1));
            tri.push_back(Eigen::Triplet<double>(3*i+2, zid, m_params(3*i+2)));
        }
        for (int i = 0; i < 3; ++i)
        {
            tri.push_back(Eigen::Triplet<double>(3*nv+i, 2*nv+1+i, 1));
        }
        Eigen::SparseMatrix<double> p2dTo3d(m_params.size(), params.size());
        p2dTo3d.setFromTriplets(tri.begin(), tri.end());
        m_params = params;
        m_paramsToCoords = m_paramsToCoords * p2dTo3d;
        m_woodU = p2dTo3d.transpose() * m_woodU;

    }

    m_paramsUndeformed = m_params;

    // m_params.tail(3) << 2,0,1;
    // m_angle = M_PI_2;
    // updateMap();
    // Eigen::VectorXd dx(m_params.size()); dx.setZero();
    // update(dx);
}

void HomogenizeStretching::visualize()
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
        if (ImGui::Button("Stretch")) {
            m_params(m_params.size()-3) += 0.05;
            Eigen::VectorXd dx(m_params.size());
            dx.setZero();
            update(dx);
            mesh->updateVertexPositions(m_tets->deformed);
        }
    };

    // Show the GUI
    polyscope::show();
}

void HomogenizeStretching::update(const Eigen::VectorXd &dx)
{
    if (dx.size() != m_params.size()) {
        throw std::runtime_error("# params neq # dx?");
    }

    m_params += dx;
    const Eigen::VectorXd coords = m_paramsToCoords * m_params;
    m_tets->deformed = coords.reshaped<Eigen::RowMajor>(coords.size()/3, 3);
}

double HomogenizeStretching::energy() const
{
    // elasticity
    double e = m_tets->energy();

    // constraint
    const Eigen::MatrixXd& V = m_tets->deformed;
    const Eigen::RowVector3d cs = (m_massMatrix * V).colwise().sum();
    e += 0.5 * m_weight * cs.squaredNorm();

    return e;
}

Eigen::VectorXd HomogenizeStretching::gradient() const
{
    // elasticity
    Eigen::VectorXd gTets = m_tets->gradient();

    // constraint
    const Eigen::MatrixXd& V = m_tets->deformed;
    const Eigen::RowVector3d cs = (m_massMatrix * V).colwise().sum();
    Eigen::MatrixXd e(V.rows(), V.cols()); e.rowwise() = cs;
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

Eigen::SparseMatrix<double> HomogenizeStretching::hessian() const
{
    Eigen::SparseMatrix<double> hess = woodburyHessian();

    hess += m_weight * m_woodU * m_woodU.transpose();

    return hess;
}

Eigen::SparseMatrix<double> HomogenizeStretching::woodburyHessian() const
{
    // elasticity
    Eigen::SparseMatrix<double> hTets = m_tets->hessian();

    // full to reduced
    Eigen::SparseMatrix<double> hess = m_paramsToCoords.transpose() * hTets * m_paramsToCoords;
   
    for (const int& idx : m_dirichletIdx) {
        hess.row(idx) *= 0.0;
        hess.col(idx) *= 0.0;
        hess.coeffRef(idx,idx) = 1.0;
    }
    return hess;
}

bool HomogenizeStretching::step(int iter)
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
        const Eigen::MatrixXd& V = m_tets->deformed;
        const Eigen::RowVector3d cs = (m_massMatrix * V).colwise().sum().cwiseAbs();
        if (cs.maxCoeff() > 1e-10 && m_weight < 1e12) {
            std::cout << "Constraints: " << cs << " not statisfied, increase weight" << std::endl;
            m_weight *= 1e2;
            return false;
        }
        if (cs.maxCoeff() > 1e-10 && m_weight == 1e12) {
            throw std::runtime_error("constraints not satisfied?");
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

bool HomogenizeStretching::sim()
{
    int iter = 0;
    while (!step(iter++));
    return iter >= m_maxIter ? false : true;
}

void HomogenizeStretching::generateData()
{
    Eigen::VectorXd dx(m_params.size()); dx.setZero();
    std::ofstream file("stretching.txt");
    auto saveData = [&]() {
        // Green strain
        const Eigen::Vector3d pbc = m_params.tail(3);
        Eigen::Matrix2d F; F << pbc(0), pbc(1), pbc(1), pbc(2);
        Eigen::Matrix2d R; R << cos(m_angle), -sin(m_angle), sin(m_angle), cos(m_angle);
        F = R * F * R.transpose();
        const Eigen::Matrix2d E = 0.5 * (F.transpose() * F - Eigen::Matrix2d::Identity());
        Eigen::Vector3d EVoigt; EVoigt << E(0,0), E(1,1), 2*E(0,1);
        
        // boundary forces
        Eigen::Vector2d f1, f2; f1.setZero(); f2.setZero();
        const Eigen::VectorXd grad = m_tets->gradient();
        for (const int& vid : m_bdyIdx1) {
            f1 += grad.segment(3 * vid, 2);
        }
        for (const int& vid : m_bdyIdx2) {
            f2 += grad.segment(3 * vid, 2);
        }
        // boundary vector is the offset of another boundary
        Eigen::Vector2d b1 = m_bbox.col(1), b2 = m_bbox.col(0);
        b1 = F * b1, b2 = F * b2;
        // force density
        f1 /= b1.norm();
        f2 /= b2.norm();
        // normal of boundary
        Eigen::Vector2d n1(b1.y(), -b1.x()), n2(-b2.y(), b2.x());
        n1.normalize(); n2.normalize();
        if (n1.dot(b2) < 0) { n1 *= -1; }
        if (n2.dot(b1) < 0) { n2 *= -1; }
        // Cauchy stress
        Eigen::Matrix2d f, n; f << f1, f2; n << n1, n2;
        const Eigen::Matrix2d sigma = f * n.inverse();
        // Second Piola Kirchhoff stress -- Green Lagrange strain / First Piola Kirchhoff stress -- deformation gradient
        const double J = F.determinant();
        const Eigen::Matrix2d S = J * (F.inverse() * sigma * F.inverse().transpose());
        Eigen::Vector3d SVoigt; SVoigt << S(0,0), S(1,1), S(0,1);

        const double area = abs(m_bbox.determinant());
        const double psi = m_tets->energy() / area;

        const double s1 = pbc(0), s2 = pbc(2);
        file << s1 << " " << s2 << " "  << m_angle << " "
            << EVoigt(0) << " " << EVoigt(1) << " " << EVoigt(2) << " " 
            << SVoigt(0) << " " << SVoigt(1) << " " << SVoigt(2) << " "
            << psi << std::endl;
    };

    const int nv = m_dim == 2 ? m_params.size() - 4 : m_params.size() - 3;;
    const int s1Idx = m_params.size() - 3;
    const int s2Idx = m_params.size() - 1;

    auto generateUni = [&](int ns, double step) {
        m_params = m_paramsUndeformed;
        m_weight = 1e6;
        for (int si = 0; si < ns; ++si) {
            const double s = step * (si + 1);

            Eigen::Matrix2d T;
            T << (1 + s) / m_params(s1Idx), 0, 0, 1;
            Eigen::Matrix2d R; R << cos(m_angle), -sin(m_angle), sin(m_angle), cos(m_angle);
            T = R * T * R.transpose();
            for (int i = 0; i < nv; i += m_dim) {
                m_params.segment(i, 2) = T * m_params.segment(i, 2);
            }
            m_params(s1Idx) = 1 + s;
            update(dx);

            if (sim())
                saveData();
            else {
                // for (int i = 0; i < nv; i += 2*m_dim) {
                //     m_params(i) += 1e-5;
                // }
                // update(dx);
                // visualize();
                m_dirichletIdx = {s1Idx, s1Idx+1, s1Idx+2};
                if (sim())
                    saveData();
                m_dirichletIdx = {s1Idx};
            }
        }
    };

    const double sStep = 0.05;
    // uniaxial stretching
    const int nUniAngle = 20;
    const int nUniStretch = 10;
    const int nUniCompress = 6;
    m_dirichletIdx = {s1Idx};
    
    for (int ai = 0; ai < nUniAngle; ++ai) {
        m_angle = ai * M_PI / nUniAngle;
        updateMap();

        // uni compress
        generateUni(nUniCompress, -sStep);
        // uni stretch
        generateUni(nUniStretch, sStep);
    }

    auto generateBi = [&](int n1, int n2, double step1, double step2) {
        m_params = m_paramsUndeformed;
        m_weight = 1e6;
        for (int si1 = 0; si1 < n1; ++si1) {
            const double s1 = step1 * (si1 + 1);

            Eigen::Matrix2d T;
            T.setIdentity();
            T(0,0) = (1 + s1) / m_params(s1Idx);
            m_params(s1Idx) = 1 + s1;

            for (int si2 = 0; si2 < n2 + 1; ++si2)
            {
                const double s2 = step2 * si2;
                T(1,1) = (1 + s2) / m_params(s2Idx);
                m_params(s2Idx) = 1 + s2;
                Eigen::Matrix2d R; R << cos(m_angle), -sin(m_angle), sin(m_angle), cos(m_angle);
                T = R * T * R.transpose();
                for (int i = 0; i < nv; i += m_dim) {
                    m_params.segment(i, 2) = T * m_params.segment(i, 2);
                }
                T.setIdentity();

                update(dx);
                if (sim())
                    saveData();
                else {
                    m_dirichletIdx = {s1Idx, s1Idx+1, s1Idx+2};
                    if (sim())
                        saveData();
                    m_dirichletIdx = {s1Idx, s1Idx+2};
                }

            }
        }
    };

    // biaxial stretching
    const int nBiAngle = 10;
    const int nBiStretch = 4;
    const int nBiCompress = 2;
    m_dirichletIdx = {s1Idx, s2Idx};
    for (int ai = 0; ai < nBiAngle; ++ai) {
        m_angle = ai * M_PI_2 / nBiAngle;
        updateMap();

        // compress - compress
        generateBi(nBiCompress, nBiCompress, -sStep, -sStep);
        // stretch - compress
        generateBi(nBiStretch, nBiCompress, sStep, -sStep);
        // compress - stretch
        generateBi(nBiCompress, nBiStretch, -sStep, sStep);
        // stretch - stretch
        generateBi(nBiStretch, nBiStretch, sStep, sStep);
    }
    file.close();

}

void HomogenizeStretching::updateMap()
{
    const int j = m_paramsToCoords.cols() - 3;
    m_paramsToCoords.rightCols(3) *= 0.0;
    const double dx1 = m_bbox(0, 0), dy1 = m_bbox(1, 0);
    for (const int& vid : m_bdyIdx1) {
        const int i = 3 * vid;

        m_paramsToCoords.coeffRef(i, j) += cos(m_angle) * cos(m_angle) * dx1 + cos(m_angle) * sin(m_angle) * dy1;
        m_paramsToCoords.coeffRef(i, j+1) += (cos(m_angle) * cos(m_angle) - sin(m_angle) * sin(m_angle)) * dy1 - 2 * cos(m_angle) * sin(m_angle) * dx1;
        m_paramsToCoords.coeffRef(i, j+2) += sin(m_angle) * sin(m_angle) * dx1 - cos(m_angle) * sin(m_angle) * dy1;
        m_paramsToCoords.coeffRef(i+1, j) += cos(m_angle) * sin(m_angle) * dx1 + sin(m_angle) * sin(m_angle) * dy1;
        m_paramsToCoords.coeffRef(i+1, j+1) += (cos(m_angle) * cos(m_angle) - sin(m_angle) * sin(m_angle)) * dx1 + 2 * cos(m_angle) * sin(m_angle) * dy1;
        m_paramsToCoords.coeffRef(i+1, j+2) += cos(m_angle) * cos(m_angle) * dy1 - cos(m_angle) * sin(m_angle) * dx1;
    }
    const double dx2 = m_bbox(0, 1), dy2 = m_bbox(1, 1);
    for (const int& vid : m_bdyIdx2) {
        const int i = 3 * vid;

        m_paramsToCoords.coeffRef(i, j) += cos(m_angle) * cos(m_angle) * dx2 + cos(m_angle) * sin(m_angle) * dy2;
        m_paramsToCoords.coeffRef(i, j+1) += (cos(m_angle) * cos(m_angle) - sin(m_angle) * sin(m_angle)) * dy2 - 2 * cos(m_angle) * sin(m_angle) * dx2;
        m_paramsToCoords.coeffRef(i, j+2) += sin(m_angle) * sin(m_angle) * dx2 - cos(m_angle) * sin(m_angle) * dy2;
        m_paramsToCoords.coeffRef(i+1, j) += cos(m_angle) * sin(m_angle) * dx2 + sin(m_angle) * sin(m_angle) * dy2;
        m_paramsToCoords.coeffRef(i+1, j+1) += (cos(m_angle) * cos(m_angle) - sin(m_angle) * sin(m_angle)) * dx2 + 2 * cos(m_angle) * sin(m_angle) * dy2;
        m_paramsToCoords.coeffRef(i+1, j+2) += cos(m_angle) * cos(m_angle) * dy2 - cos(m_angle) * sin(m_angle) * dx2;
    }
}