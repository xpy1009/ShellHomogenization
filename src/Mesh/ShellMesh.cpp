#include <Mesh/ShellMesh.h>
#include <Mesh/MeshLib.h>
#include <iostream>
#include <gmsh.h>

constexpr int IH = 1;

void ShellMesh::getVF(Eigen::MatrixX3d &vertices, Eigen::MatrixX3i &faces)
{
    // vertices
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);
    vertices = Eigen::Map<Eigen::MatrixXd>(nodeCoords.data(), 3, nodeCoords.size()/3).transpose();

    // faces 
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, _elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, _elemNodeTags, 2);
    if (_elemNodeTags.size() != 1) {
        throw std::runtime_error("# Surface != 1");
    }

    const std::vector<size_t> &elemNodeTags = _elemNodeTags[0];
    constexpr int n = 3;
    const int nTris = elemNodeTags.size() / n;
    faces.resize(nTris, n);
    // https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    for (int i = 0; i < nTris; ++i) {
        faces.row(i) << elemNodeTags[n*i], elemNodeTags[n*i+1], elemNodeTags[n*i+2];
    }
    faces.array() -= 1.0;
}

void ShellMesh::biMatShell(Eigen::MatrixX3d &vertices, Eigen::MatrixX3i &faces)
{
    Eigen::Matrix2d bbox; 
    std::vector<double> tile;
    std::vector<std::vector<double>> tiles;
    bbox << -1, -1, 6, 6;
    MeshLib::constructTiling(IH, bbox, tiles);
    std::cout << "bbox: " << std::endl << bbox << std::endl;

    std::vector<std::vector<double>> paths = MeshLib::inflatePaths(tiles);

    gmsh::initialize();
    gmsh::model::add("biMatShell");
    namespace factory = gmsh::model::occ;

    // square sheet
    int nx = 6;
    int ny = 5;
    double dx = abs(nx * bbox(0,1));
    double dy = abs(ny * bbox(1,0));
    factory::addRectangle(0.0, 0.0, 0.0, dx, dy);

    // tiling polygons
    std::vector<std::pair<int, int>> tilingDimTags(paths.size());
    for (size_t i = 0; i < paths.size(); ++i) {
        std::vector<int> pids(paths[i].size() / 2);
        for (size_t j = 0; j < pids.size(); ++j) {
            pids[j] = factory::addPoint(paths[i][2 * j], paths[i][2 * j + 1], 0.0);
        }
        std::vector<int> ls(pids.size());
        for (size_t j = 0; j < ls.size(); ++j) {
            ls[j] = factory::addLine(pids[j], pids[(j+1)%pids.size()]);
        }
        const int clid = factory::addCurveLoop(ls);
        const int psid = factory::addPlaneSurface({clid});
        tilingDimTags[i] = std::make_pair(2, psid);
    }

    std::vector<std::pair<int, int>> ov1, ov2;
    std::vector<std::vector<std::pair<int, int>>> ovv1, ovv2;
    factory::intersect({{2, 1}}, tilingDimTags, ov1, ovv1, -1, false);
    factory::fragment({{2,1}}, ov1, ov2, ovv2);

    factory::synchronize();

    gmsh::option::setNumber("General.Verbosity", 2.0);

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::reverse(ov1);

    getVF(vertices, faces);

    gmsh::fltk::run();
    gmsh::finalize();
    // std::exit(0);
}

void ShellMesh::biMatShellPeriodic(Eigen::MatrixX3i &faces, 
                                Eigen::SparseMatrix<double>& paramsToCoords, 
                                Eigen::VectorXd& params, 
                                Eigen::Matrix2d& bbox, 
                                std::vector<int>* bdyIdx1, 
                                std::vector<int>* bdyIdx2)
{
    std::vector<std::vector<double>> tiles;
    bbox << -5,-5,5,5;
    MeshLib::constructTiling(IH, bbox, tiles);
    std::cout << "bbox: " << std::endl << bbox << std::endl;

    std::vector<std::vector<double>> paths = MeshLib::inflatePaths(tiles);

    gmsh::initialize();
    gmsh::model::add("biMatShellPeriodic");
    namespace factory = gmsh::model::occ;

    // periodic region
    double dx1 = bbox(0,0), dx2 = bbox(0,1), dy1 = bbox(1,0), dy2 = bbox(1,1);
    double dx = abs(dx1 + dx2);
    double dy = abs(dy1 + dy2);
    factory::addPoint(dx, dy, 0.0, 0.0, 1);
    factory::addPoint(dx1 + dx, dy1 + dy, 0.0, 0.0, 2);
    factory::addPoint(dx1 + dx2 + dx, dy1 + dy2 + dy, 0.0, 0.0, 3);
    factory::addPoint(dx2 + dx, dy2 + dy, 0.0, 0.0, 4);
    factory::addLine(1, 2, 1);
    factory::addLine(2, 3, 2);
    factory::addLine(3, 4, 3);
    factory::addLine(4, 1, 4);
    factory::addCurveLoop({1, 2, 3, 4}, 1);
    factory::addPlaneSurface({1}, 1);

    // tiling polygons
    std::vector<std::pair<int, int>> tilingDimTags(paths.size());
    for (size_t i = 0; i < paths.size(); ++i) {
        std::vector<int> pids(paths[i].size() / 2);
        for (size_t j = 0; j < pids.size(); ++j) {
            pids[j] = factory::addPoint(paths[i][2 * j], paths[i][2 * j + 1], 0.0);
        }
        std::vector<int> ls(pids.size());
        for (size_t j = 0; j < ls.size(); ++j) {
            ls[j] = factory::addLine(pids[j], pids[(j + 1) % pids.size()]);
        }
        const int clid = factory::addCurveLoop(ls);
        const int psid = factory::addPlaneSurface({clid});
        tilingDimTags[i] = std::make_pair(2, psid);
    }

    std::vector<std::pair<int, int>> ov1, ov2;
    std::vector<std::vector<std::pair<int, int>>> ovv1, ovv2;
    factory::intersect({{2, 1}}, tilingDimTags, ov1, ovv1, -1, false);
    factory::fragment({{2,1}}, ov1, ov2, ovv2);

    factory::synchronize();
    // gmsh::option::setNumber("Mesh.MeshSizeFactor", 2);

    // get surface
    std::vector<std::pair<int, int>> meshes, boundaryDimtags;
    gmsh::model::getEntities(meshes, 2);
    gmsh::model::getBoundary(meshes, boundaryDimtags, true, false);
    const std::vector<std::pair<int, int>> bdyDimTags1 = MeshLib::setPeriodic(boundaryDimtags, dx1, dy1);
    const std::vector<std::pair<int, int>> bdyDimTags2 = MeshLib::setPeriodic(boundaryDimtags, dx2, dy2);

    gmsh::option::setNumber("General.Verbosity", 2.0);

    gmsh::model::mesh::generate(2);

    Eigen::MatrixX3d vertices;
    getVF(vertices, faces);

    // build reduced coordinates
    MeshLib::buildReduced(vertices, bdyDimTags1, bdyDimTags2, paramsToCoords, params, bdyIdx1, bdyIdx2);

    gmsh::fltk::run();
    gmsh::finalize();
    std::exit(0);
}