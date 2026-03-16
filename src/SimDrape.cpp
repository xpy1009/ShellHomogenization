// #include <SimDrape.h>
// #include <MeshLib.h>
// #include <Newton.h>

// #include <Eigen/Dense>
// #include <igl/writeOBJ.h>
// #include <igl/readOBJ.h>
// #include <igl/point_mesh_squared_distance.h>
// #include <igl/biharmonic_coordinates.h>

// #include <polyscope/polyscope.h>
// #include <polyscope/surface_mesh.h>
// #include <polyscope/point_cloud.h>

// SimDrape::SimDrape()
// {
//     Eigen::MatrixXd V;
//     Eigen::MatrixXi T;
//     MeshLib::structuredSheet(V, T);

//     // MeshLib::generateBox(V, T);

//     m_tets = std::make_unique<Tetrahedron>(V, T);
//     m_massMatrix = m_tets->massMatrix() * m_tets->volume();

//     const Eigen::MatrixXd& coords =  m_tets->deformed;
//     double x = coords.col(0).maxCoeff(), y = coords.col(1).maxCoeff();
//     for (int i = 0; i < coords.rows(); i++) {
//         // if (abs(coords(i,2)) < 1e-4 && ((abs(coords(i,0)) < 5e-3 && abs(coords(i,1)) < 5e-3) || (abs(coords(i,0)-x) < 5e-3 && abs(coords(i,1)-y) < 5e-3)))  {
//         if (abs(coords(i,2)) < 1e-4 && abs(coords(i,0)-x/2) < x/6+1e-5 && abs(coords(i,1)-y/2) < y/6+1e-5)  {
//             dirichlet_idx.push_back(3*i);
//             dirichlet_idx.push_back(3*i+1);
//             dirichlet_idx.push_back(3*i+2);
//         } 
//     }
//     std::cout << std::endl << "Volume: " << m_tets->volume() << std::endl;

// }

// void SimDrape::visualize()
// {    
//     // Initialize Polyscope
//     polyscope::view::setUpDir(polyscope::UpDir::ZUp);
//     polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
//     polyscope::init();

//     // Register the mesh with Polyscope
//     polyscope::SurfaceMesh* mesh = polyscope::registerSurfaceMesh("input mesh", m_tets->deformed, m_tets->faces);
//     mesh->setEdgeWidth(1.0);
//     mesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
//     mesh->setBackFaceColor({0.2, 0.2, 0.2});

//     Eigen::MatrixXd vd(dirichlet_idx.size()/3, 3);
//     for (int i = 0; i < dirichlet_idx.size()/3; i++)
//     {
//         int vid = dirichlet_idx[3*i]/3;
//         vd.row(i) = m_tets->deformed.row(vid);
//     }
//     polyscope::PointCloud* psCloud = polyscope::registerPointCloud("origin", vd);

//     polyscope::state::userCallback = [&]() 
//     {
//         ImGui::PushItemWidth(100);
//         if (ImGui::Button("Optimize")) {
//             sim();
//             mesh->updateVertexPositions(m_tets->deformed);
//         }
//         if (ImGui::Button("Step")) {
//             step(0);
//             mesh->updateVertexPositions(m_tets->deformed);
//         }
//         if (ImGui::Button("Save")) {
//             igl::writeOBJ("tets.obj", m_tets->deformed, m_tets->faces);
//         }
//     };

//     // Show the GUI
//     polyscope::show();
// }

// void SimDrape::update(const Eigen::VectorXd &dx)
// {
//     m_tets->deformed += dx.reshaped<Eigen::RowMajor>(dx.size()/3, 3);
// }

// double SimDrape::energy() const
// {
//     double Ee = m_tets->energy();
//     double Eg = rho_g * (m_massMatrix * m_tets->deformed.col(2)).sum();
//     return Ee + Eg;
// }

// Eigen::VectorXd SimDrape::gradient() const
// {
//     Eigen::VectorXd Ge = m_tets->gradient();

//     Eigen::MatrixXd coords = m_tets->deformed;
//     coords.leftCols(2).setZero();
//     coords.col(2).setOnes();
//     Eigen::VectorXd Gg = rho_g * (m_massMatrix * coords).reshaped<Eigen::RowMajor>();
//     Eigen::VectorXd G = Ge + Gg;
//     for(const int& idx : dirichlet_idx) {
//         G(idx) = 0.0;
//     }
//     return G;
// }

// Eigen::SparseMatrix<double> SimDrape::hessian() const
// {
//     Eigen::SparseMatrix<double> h = m_tets->hessian();
//     for (const int& idx : dirichlet_idx) {
//         h.row(idx) *= 0.0;
//         h.col(idx) *= 0.0;
//         h.coeffRef(idx,idx) = 1.0;
//     }

//     return h;
// }

// bool SimDrape::step(int iter)
// {
//     Eigen::VectorXd grad = gradient();
    
//     double grad_norm = grad.norm();
//     std::cout << "[NEWTON] iter " << iter << "/" << max_iter << " gradient norm: " << grad_norm << " tol: " << tol << std::endl;

//     if (grad_norm < tol || iter == max_iter) {
//         return true;
//     }

//     // linear solve
//     Eigen::SparseMatrix<double> K = hessian();
//     const Eigen::VectorXd dir = Newton::solve(K, -grad);

//     // line search
//     const double E0 = energy();
//     double step_size = 1;
//     int cnt = 0;

//     while (true) {
//         Eigen::VectorXd dx = step_size * dir;
//         update(dx);

//         const double E1 = energy();

//         if (E1 - E0 < 0 || cnt > 15) {
//             if (cnt > 15) {
//                 std::cout << "line search max" << std::endl;
//             }
//             break;
//         }

//         dx *= -1;
//         update(dx);
//         step_size *= 0.5;
//         cnt++;
//     }

//     return false;
// }

// void SimDrape::sim()
// {
//     int iter = 0;
//     while (!step(iter++));
// }

// void SimDrape::guess()
// {
//     Eigen::MatrixXd V;
//     Eigen::MatrixXi F;
//     const Eigen::MatrixXd& Vn = m_tets->undeformed;
//     const Eigen::MatrixXi& Fn = m_tets->faces;
//     // mesh for midsurface with similar topology
//     MeshLib::structuredSheetMidsurface(V, F);
//     // find correspondence
//     Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(Vn.rows(), 0, Vn.rows()-1);
//     Eigen::VectorXd sqrD, I;
//     Eigen::MatrixXd C;
//     igl::point_mesh_squared_distance(V, Vn, Ele, sqrD, I, C);
//     std::vector<int> ind, indmid;
//     std::vector<std::vector<int>> S;
//     for (int i = 0; i < sqrD.size(); ++i) {
//         if (sqrD(i) < 1e-7) {
//             ind.push_back(I(i));
//             indmid.push_back(i);
//             S.push_back({i});
//         }
//     }
//     // deform based on biharmonic coordinates
//     Eigen::MatrixXd W;
//     igl::biharmonic_coordinates(V, F, S, W);
//     const Eigen::MatrixXd& Vnd = m_tets->deformed;
//     Eigen::MatrixXd Vd = W * Vnd(ind, Eigen::all);

//     // initial guess to eliminate multistability
//     Eigen::MatrixXd Vh;
//     Eigen::MatrixXi Fh;
//     igl::readOBJ("../data/grid40.obj", Vh, Fh);
//     Vh.col(0).array() += 0.1;
//     Vh *= 10;
//     Vh.col(0).array() *= V.col(0).maxCoeff();
//     Vh.col(1).array() *= V.col(1).maxCoeff();
//     igl::writeOBJ("1.obj", Vd, F);
//     // find barycentric coordinates
//     {
//         Eigen::VectorXd sqrD;
//         Eigen::VectorXi I;
//         Eigen::MatrixXd C;
//         igl::point_mesh_squared_distance(Vh, V, F, sqrD, I, C);
//         Eigen::MatrixXd L;
//         igl::barycentric_coordinates(Vh, V(F(I, 0), Eigen::all), V(F(I, 1), Eigen::all), V(F(I, 2), Eigen::all), L);
//         Vh = (Vd(F(I, 0), Eigen::all).array().colwise() * L.col(0).array()).matrix()
//             +(Vd(F(I, 1), Eigen::all).array().colwise() * L.col(1).array()).matrix()
//             +(Vd(F(I, 2), Eigen::all).array().colwise() * L.col(2).array()).matrix();
//     }
//     igl::writeOBJ("2.obj", Vh, Fh);
// }