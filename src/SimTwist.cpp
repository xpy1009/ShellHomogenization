// #include <SimTwist.h>
// #include <Mesh/MeshLib.h>
// #include <Newton.h>

// #include <Eigen/Dense>
// #include <igl/writeOBJ.h>
// #include <igl/readOBJ.h>

// #include <polyscope/polyscope.h>
// #include <polyscope/surface_mesh.h>
// #include <polyscope/point_cloud.h>

// SimTwist::SimTwist()
// {
//     Eigen::MatrixXd V;
//     Eigen::MatrixXi T;
//     MeshLib::structuredSheet(V, T);

//     // V *= 1e-2;
//     // V.col(2).array() -= V.col(2).mean();
//     // V.leftCols(2) = V.leftCols(2) * Eigen::Rotation2Dd(M_PI_2).matrix();

//     // GeometryLib::generateBox(V, T);

//     tets = std::make_unique<Tetrahedron>(V, T);

//     Eigen::MatrixXd& coords =  tets->deformed;
//     double min_x = coords.col(0).minCoeff(), max_x = coords.col(0).maxCoeff();
//     Eigen::Vector3d v = 0.5 * coords.col(1).maxCoeff() * Eigen::Vector3d::UnitY();
//     std::cout << v << std::endl << max_x << std::endl;
//     for (int i = 0; i < coords.rows(); i++)
//     {
//         if (abs(coords(i,0)-min_x) < 5e-3 || abs(coords(i,0)-max_x) < 5e-3) 
//         {
//             dirichlet_idx.push_back(3*i);
//             dirichlet_idx.push_back(3*i+1);
//             dirichlet_idx.push_back(3*i+2);
//         } 
//         double angle = coords(i,0) / max_x * M_PI_2;
//         auto R = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX());
//         coords.row(i) =( R * (coords.row(i).transpose() - v) + v).transpose();
//     }
//     coords.col(0) *= 0.9;
//     std::cout << tets->deformed.col(0).minCoeff() << " " << tets->deformed.col(0).maxCoeff() << " "
//         << tets->deformed.col(1).minCoeff() << " " << tets->deformed.col(1).maxCoeff() << std::endl;
    
//     // tets->deformed.col(1) += 0.1 * tets->deformed.col(0);
//     // std::cout <<  << std::endl;

//     // {
//     //     Eigen::MatrixXd V;
//     //     Eigen::MatrixXi F;
//     //     igl::readOBJ("tets.obj", V, F);
//     //     coords = V;
//     // }

// }

// void SimTwist::visualize()
// {    
//     // Initialize Polyscope
//     polyscope::view::setUpDir(polyscope::UpDir::ZUp);
//     polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
//     polyscope::init();

//     // Register the mesh with Polyscope
//     polyscope::SurfaceMesh* mesh = polyscope::registerSurfaceMesh("input mesh", tets->deformed, tets->faces);
//     mesh->setEdgeWidth(1.0);
//     mesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
//     mesh->setBackFaceColor({0.2, 0.2, 0.2});

//     Eigen::MatrixXd vd(dirichlet_idx.size()/3, 3);
//     for (int i = 0; i < dirichlet_idx.size()/3; i++)
//     {
//         int vid = dirichlet_idx[3*i]/3;
//         vd.row(i) = tets->deformed.row(vid);
//     }
//     polyscope::PointCloud* psCloud = polyscope::registerPointCloud("origin", vd);

//     polyscope::state::userCallback = [&]() 
//     {
//         ImGui::PushItemWidth(100);
//         if (ImGui::Button("Optimize")) {
//             sim();
//             mesh->updateVertexPositions(tets->deformed);
//         }
//         if (ImGui::Button("Step")) {
//             step(0);
//             mesh->updateVertexPositions(tets->deformed);
//         }
//         if (ImGui::Button("Next")) {
//             tets->deformed.col(0) *= 0.9;
//             std::cout << tets->deformed.col(0).maxCoeff() << std::endl;
//             mesh->updateVertexPositions(tets->deformed);
//             for (int i = 0; i < dirichlet_idx.size()/3; i++)
//             {
//                 int vid = dirichlet_idx[3*i]/3;
//                 vd.row(i) = tets->deformed.row(vid);
//             }
//             psCloud->updatePointPositions(vd);
//         }
//         if (ImGui::Button("Save")) {
//             igl::writeOBJ("tets.obj", tets->deformed, tets->faces);
//         }
//     };

//     // Show the GUI
//     polyscope::show();
// }

// void SimTwist::update(const Eigen::VectorXd &dx)
// {
//     tets->deformed += dx.reshaped<Eigen::RowMajor>(dx.size()/3, 3);
// }

// double SimTwist::energy() const
// {
//     double e = tets->energy();
//     return e;
// }

// Eigen::VectorXd SimTwist::gradient() const
// {
//     Eigen::VectorXd g = tets->gradient();
//     for(const int& idx : dirichlet_idx) {
//         g(idx) = 0.0;
//     }
//     return g;
// }

// Eigen::SparseMatrix<double> SimTwist::hessian() const
// {
//     Eigen::SparseMatrix<double> h = tets->hessian();
//     for (const int& idx : dirichlet_idx) {
//         h.row(idx) *= 0.0;
//         h.col(idx) *= 0.0;
//         h.coeffRef(idx,idx) = 1.0;
//     }

//     return h;
// }

// bool SimTwist::step(int iter)
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

// void SimTwist::sim()
// {
//     int iter = 0;
//     while (!step(iter++));
// }