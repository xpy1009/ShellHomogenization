// #include <SimPress.h>
// #include <MeshLib.h>
// #include <Newton.h>
// #include <GeometryLib.h>

// #include <Eigen/Dense>
// #include <igl/writeOBJ.h>
// #include <igl/readOBJ.h>
// #include <igl/point_mesh_squared_distance.h>
// #include <igl/biharmonic_coordinates.h>

// #include <polyscope/polyscope.h>
// #include <polyscope/surface_mesh.h>
// #include <polyscope/point_cloud.h>

// SimPress::SimPress()
// {
//     Eigen::MatrixXd V;
//     Eigen::MatrixXi T;
//     MeshLib::structuredSheet(V, T);

//     m_tets = std::make_unique<Tetrahedron>(V, T);

//     Eigen::MatrixXd& coords =  m_tets->deformed;
//     double min_x = coords.col(0).minCoeff(), max_x = coords.col(0).maxCoeff();
//     double min_y = coords.col(1).minCoeff(), max_y = coords.col(1).maxCoeff();
//     double mid_x = 0.5 * (min_x + max_x);

//     double h = 1e-2;
//     for (int i = 0; i < coords.rows(); i++)
//     {
//         if (abs(coords(i,0)-mid_x) < 7e-3 && abs(coords(i,2)) < 1e-5 && (abs(coords(i,1)-min_y) < 1e-5 || abs(coords(i,1)-max_y) < 1e-5)) 
//         {
//             dirichlet_idx.push_back(3*i);
//             dirichlet_idx.push_back(3*i+1);
//             dirichlet_idx.push_back(3*i+2);
//         }
//         coords(i, 2) += h * sin(coords(i, 1) / max_y * M_PI);
//         coords(i, 1) *= .9;
//     }
// }

// void SimPress::visualize()
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
//         if (ImGui::Button("Next")) {
//             m_tets->deformed.col(1) *= 0.9;
//             mesh->updateVertexPositions(m_tets->deformed);
//             for (int i = 0; i < dirichlet_idx.size()/3; i++)
//             {
//                 int vid = dirichlet_idx[3*i]/3;
//                 vd.row(i) = m_tets->deformed.row(vid);
//             }
//             psCloud->updatePointPositions(vd);
//         }
//         if (ImGui::Button("Save")) {
//             igl::writeOBJ("tets.obj", m_tets->deformed, m_tets->faces);
//         }
//     };

//     // Show the GUI
//     polyscope::show();
// }

// void SimPress::guess()
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
//     igl::writeOBJ("m.obj", Vd, F);
// }

// void SimPress::update(const Eigen::VectorXd &dx)
// {
//     m_tets->deformed += dx.reshaped<Eigen::RowMajor>(dx.size()/3, 3);
// }

// double SimPress::energy() const
// {
//     double e = m_tets->energy();
//     return e;
// }

// Eigen::VectorXd SimPress::gradient() const
// {
//     Eigen::VectorXd g = m_tets->gradient();
//     for(const int& idx : dirichlet_idx) {
//         g(idx) = 0.0;
//     }
//     return g;
// }

// Eigen::SparseMatrix<double> SimPress::hessian() const
// {
//     Eigen::SparseMatrix<double> h = m_tets->hessian();
//     for (const int& idx : dirichlet_idx) {
//         h.row(idx) *= 0.0;
//         h.col(idx) *= 0.0;
//         h.coeffRef(idx,idx) = 1.0;
//     }

//     return h;
// }

// bool SimPress::step(int iter)
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

// void SimPress::sim()
// {
//     int iter = 0;
//     while (!step(iter++));
// }

