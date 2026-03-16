// #include <SimDrapeShell.h>
// // #include <MeshLib.h>
// // #include <Newton.h>

// // #include <Eigen/Dense>
// // #include <igl/writeOBJ.h>
// #include <igl/readOBJ.h>

// // #include <polyscope/polyscope.h>
// // #include <polyscope/surface_mesh.h>
// // #include <polyscope/point_cloud.h>

// SimDrapeShell::SimDrapeShell()
// {
//     Eigen::MatrixXd V;
//     Eigen::MatrixXi F;
//     igl::readOBJ("../data/grid40.obj", V, F);
//     V.col(0).array() += 0.1;
//     V *= 10;

//     // MeshLib::generateBox(V, T);

//     m_shell = std::make_unique<DiscreteShell>(V, F);
//     // m_massMatrix = m_tets->massMatrix() * m_tets->volume();

//     // const Eigen::MatrixXd& coords =  m_tets->deformed;
//     // for (int i = 0; i < coords.rows(); i++) {
//     //     if (abs(coords(i,0)) < 1e-5)  {
//     //         dirichlet_idx.push_back(3*i);
//     //         dirichlet_idx.push_back(3*i+1);
//     //         dirichlet_idx.push_back(3*i+2);
//     //     } 
//     // }

// }

// void SimDrapeShell::visualize()
// {    
//     // // Initialize Polyscope
//     // polyscope::view::setUpDir(polyscope::UpDir::ZUp);
//     // polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
//     // polyscope::init();

//     // // Register the mesh with Polyscope
//     // polyscope::SurfaceMesh* mesh = polyscope::registerSurfaceMesh("input mesh", m_tets->deformed, m_tets->faces);
//     // mesh->setEdgeWidth(1.0);
//     // mesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
//     // mesh->setBackFaceColor({0.2, 0.2, 0.2});

//     // Eigen::MatrixXd vd(dirichlet_idx.size()/3, 3);
//     // for (int i = 0; i < dirichlet_idx.size()/3; i++)
//     // {
//     //     int vid = dirichlet_idx[3*i]/3;
//     //     vd.row(i) = m_tets->deformed.row(vid);
//     // }
//     // polyscope::PointCloud* psCloud = polyscope::registerPointCloud("origin", vd);

//     // polyscope::state::userCallback = [&]() 
//     // {
//     //     ImGui::PushItemWidth(100);
//     //     if (ImGui::Button("Optimize")) {
//     //         sim();
//     //         mesh->updateVertexPositions(m_tets->deformed);
//     //     }
//     //     if (ImGui::Button("Step")) {
//     //         step(0);
//     //         mesh->updateVertexPositions(m_tets->deformed);
//     //     }
//     //     if (ImGui::Button("Next")) {
//     //         m_tets->deformed.col(0) *= 0.9;
//     //         std::cout << m_tets->deformed.col(0).maxCoeff() << std::endl;
//     //         mesh->updateVertexPositions(m_tets->deformed);
//     //         for (int i = 0; i < dirichlet_idx.size()/3; i++)
//     //         {
//     //             int vid = dirichlet_idx[3*i]/3;
//     //             vd.row(i) = m_tets->deformed.row(vid);
//     //         }
//     //         psCloud->updatePointPositions(vd);
//     //     }
//     //     if (ImGui::Button("Save")) {
//     //         igl::writeOBJ("tets.obj", m_tets->deformed, m_tets->faces);
//     //     }
//     // };

//     // // Show the GUI
//     // polyscope::show();
// }

// void SimDrapeShell::update(const Eigen::VectorXd &dx)
// {
//     // m_tets->deformed += dx.reshaped<Eigen::RowMajor>(dx.size()/3, 3);
// }

// double SimDrapeShell::energy() const
// {
//     return 0;
//     // double Ee = m_tets->energy();
//     // double Eg = rho * (m_massMatrix * m_tets->deformed.col(2)).sum();
//     // return Ee + Eg;
// }

// Eigen::VectorXd SimDrapeShell::gradient() const
// {
//     // Eigen::VectorXd Ge = m_tets->gradient();

//     // Eigen::MatrixXd coords = m_tets->deformed;
//     // coords.leftCols(2).setZero();
//     // coords.col(2).setOnes();
//     // Eigen::VectorXd Gg = rho * (m_massMatrix * coords).reshaped<Eigen::RowMajor>();
//     // Eigen::VectorXd G = Ge + Gg;
//     // for(const int& idx : dirichlet_idx) {
//     //     G(idx) = 0.0;
//     // }
//     // return G;
//     return Eigen::VectorXd();
// }

// Eigen::SparseMatrix<double> SimDrapeShell::hessian() const
// {
//     Eigen::SparseMatrix<double> h;
//     //  = m_tets->hessian();
//     // for (const int& idx : dirichlet_idx) {
//     //     h.row(idx) *= 0.0;
//     //     h.col(idx) *= 0.0;
//     //     h.coeffRef(idx,idx) = 1.0;
//     // }

//     return h;
// }

// bool SimDrapeShell::step(int iter)
// {
//     // Eigen::VectorXd grad = gradient();
    
//     // double grad_norm = grad.norm();
//     // std::cout << "[NEWTON] iter " << iter << "/" << max_iter << " gradient norm: " << grad_norm << " tol: " << tol << std::endl;

//     // if (grad_norm < tol || iter == max_iter) {
//     //     return true;
//     // }

//     // // linear solve
//     // Eigen::SparseMatrix<double> K = hessian();
//     // const Eigen::VectorXd dir = Newton::solve(K, -grad);

//     // // line search
//     // const double E0 = energy();
//     // double step_size = 1;
//     // int cnt = 0;

//     // while (true) {
//     //     Eigen::VectorXd dx = step_size * dir;
//     //     update(dx);

//     //     const double E1 = energy();

//     //     if (E1 - E0 < 0 || cnt > 15) {
//     //         if (cnt > 15) {
//     //             std::cout << "line search max" << std::endl;
//     //         }
//     //         break;
//     //     }

//     //     dx *= -1;
//     //     update(dx);
//     //     step_size *= 0.5;
//     //     cnt++;
//     // }

//     return false;
// }

// void SimDrapeShell::sim()
// {
//     // int iter = 0;
//     // while (!step(iter++));
// }