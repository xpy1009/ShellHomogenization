#include <iostream>
#include <chrono>

#include <SimPress.h>
#include <SimTwist.h>
#include <SimDrape.h>
#include <HomogenizeBending.h>
#include <HomogenizeStretching.h>
#include <DerivativeTest.h>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include <DummyObjective.h>
#include <SimShell.h>

int main(int argc, char **argv)
{
    // DummyObjective obj;
    // DerivativeTest::checkGradient(obj);
    // DerivativeTest::checkHessian(obj);
    SimShell sim;
    sim.visualize();
    // DerivativeTest::checkGradient(sim);
    // DerivativeTest::checkHessian(sim);

    // HomogenizeStretching sim;
    // HomogenizeBending sim;
    // SimTwist sim;
    // SimDrape sim;
    // SimPress sim;

    // sim.generateData();
    // sim.generateDataHYLC();

    // DerivativeTest::checkGradient(sim);
    // DerivativeTest::checkHessian(sim);
    
    // sim.visualize();
    // sim.guess();

    // // Initialize Polyscope
    // polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    // polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
    // polyscope::init();

    // Eigen::MatrixXd Vn, V1, V2, Vmid;
    // Eigen::MatrixXi Fn, F1, F2, Fmid;
    // igl::readOBJ("o3.obj", V1, F1);
    // igl::readOBJ("n3.obj", V2, F2);
    // igl::readOBJ("m3.obj", Vmid, Fmid);
    // igl::readOBJ("3.obj", Vn, Fn);
    // Eigen::VectorXd sqrD;
    // Eigen::VectorXi I;
    // Eigen::MatrixXd C;
    // igl::point_mesh_squared_distance(V1, Vmid, Fmid, sqrD, I, C);
    // Eigen::VectorXd D1 = 100 * sqrD.cwiseSqrt();
    // igl::point_mesh_squared_distance(V2, Vmid, Fmid, sqrD, I, C);
    // Eigen::VectorXd D2 = 100 * sqrD.cwiseSqrt();
    // std::ofstream f1("do3.txt");
    // for (int i = 0; i < D1.size(); i++) {
    //     f1 << D1(i) << std::endl;
    // }
    // f1.close();
    // std::ofstream f2("dn3.txt");
    // for (int i = 0; i < D2.size(); i++) {
    //     f2 << D2(i) << std::endl;
    // }
    // f2.close();

    // // Register the mesh with Polyscope
    // polyscope::SurfaceMesh* m1 = polyscope::registerSurfaceMesh("1", V1, F1);
    // m1->setEdgeWidth(1.0);
    // m1->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    // m1->setBackFaceColor({0.2, 0.2, 0.2});
    // m1->addVertexScalarQuantity("dist", D1);

    // polyscope::SurfaceMesh* m2 = polyscope::registerSurfaceMesh("2", V2, F2);
    // m2->setEdgeWidth(1.0);
    // m2->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    // m2->setBackFaceColor({0.2, 0.2, 0.2});
    // m2->addVertexScalarQuantity("dist", D2);

    // polyscope::SurfaceMesh* m3 = polyscope::registerSurfaceMesh("3", Vn, Fn);
    // m3->setEdgeWidth(1.0);
    // m3->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    // m3->setBackFaceColor({0.2, 0.2, 0.2});
    
    // polyscope::SurfaceMesh* m4 = polyscope::registerSurfaceMesh("4", Vmid, Fmid);
    // m4->setEdgeWidth(1.0);
    // m4->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    // m4->setBackFaceColor({0.2, 0.2, 0.2});

    // polyscope::show();

    return 0;
}