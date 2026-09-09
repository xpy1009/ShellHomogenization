#include <Mesh/MeshLib.h>
#include <Mesh/Isohedral.h>

#include <gmsh.h>
#include <igl/remove_duplicate_vertices.h>
#include <iostream>


void MeshLib::getQuadTets(Eigen::MatrixX3d &vertices, Tet10 &tets)
{
    // vertices
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);
    vertices = Eigen::Map<Eigen::MatrixXd>(nodeCoords.data(), 3, nodeCoords.size()/3).transpose();

    // tets 
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, _elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, _elemNodeTags, 3);
    if (_elemNodeTags.size() < 1) {
        std::cout << "Not tets" << std::endl;
        return;
    }
    const std::vector<size_t> &elemNodeTags = _elemNodeTags[0];
    constexpr int n = 10;
    const int nTets = elemNodeTags.size() / n;
    tets.resize(nTets, n);
    // From https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering to https://www.sd.ruhr-uni-bochum.de/downloads/Shape_funct.pdf
    for (int i = 0; i < nTets; ++i) {
         tets.row(i) << elemNodeTags[n*i], elemNodeTags[n*i+1], elemNodeTags[n*i+2], 
                        elemNodeTags[n*i+3], elemNodeTags[n*i+4], elemNodeTags[n*i+5], 
                        elemNodeTags[n*i+6], elemNodeTags[n*i+7], elemNodeTags[n*i+9], elemNodeTags[n*i+8];
    }
    tets.array() -= 1.0;
}

void MeshLib::rectangle(double dx, double dy, Eigen::MatrixX2d &vertices, Eigen::MatrixX3i &faces, bool preview)
{
    gmsh::initialize();
    gmsh::model::add("rectangle");
    namespace factory = gmsh::model::occ;

    factory::addRectangle(0,0,0, dx,dy);

    factory::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeFactor", .2);
    
    gmsh::option::setNumber("General.Verbosity", 2.0); // silent gmsh
    gmsh::model::mesh::setTransfiniteAutomatic({}, 2.35, false); // regular
    gmsh::model::mesh::generate(2); // 2d

    // vertices
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);
    vertices = Eigen::Map<Eigen::MatrixXd>(nodeCoords.data(), 3, nodeCoords.size()/3).transpose().leftCols(2);

    
    // triangles 
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, _elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, _elemNodeTags, 2);
    if (_elemNodeTags.size() < 1) {
        std::cout << "No triangles" << std::endl;
        return;
    }
    const std::vector<size_t> &elemNodeTags = _elemNodeTags[0];
    constexpr int n = 3;
    const int nTris = elemNodeTags.size() / n;
    faces.resize(nTris, n);
    for (int i = 0; i < nTris; ++i) {
        faces.row(i) << elemNodeTags[n*i], elemNodeTags[n*i+1], elemNodeTags[n*i+2];
    }
    faces.array() -= 1.0;

    if (preview) {
        gmsh::fltk::run();
    }
    gmsh::finalize();
}

void MeshLib::structuredSheet(
    int IH, 
    const std::vector<double>& params,
    double inflateWidth,
    const std::array<int, 2>& size,
    double height,
    Eigen::MatrixX3d &vertices, 
    Tet10 &tets,
    bool preview)
{
    const std::array<int, 4> nTiles = {-1, -1, size[0], size[1]};
    std::array<double, 4> trans;
    std::vector<std::vector<double>> tiles = Isohedral::getTiling(IH, params, nTiles, trans);

    std::vector<std::vector<double>> paths = Isohedral::inflate(tiles, 1e2 * inflateWidth); // unit in cm

    // set unit to m
    for (double &t : trans) {
        t *= 1e-2;
    }
    for (auto &p : paths) {
        for (double &v : p) {
            v *= 1e-2;
        }
    }

    gmsh::initialize();
    gmsh::model::add("structured sheet");
    gmsh::option::setNumber("General.Verbosity", 2.0);
    namespace factory = gmsh::model::occ;

    // square sheet
    const int &nx = size[0], &ny = size[1];
    const double dx = nx * std::max(abs(trans[0]), abs(trans[2]));
    const double dy = ny * std::max(abs(trans[1]), abs(trans[3]));
    factory::addRectangle(0.0, 0.0, 0.0, dx, dy);

    // tiling polygons
    std::vector<std::pair<int, int>> tilingDimTags(paths.size());
    for (int i = 0; i < paths.size(); ++i) {
        std::vector<int> pids(paths[i].size()/2);
        for (int j = 0; j < pids.size(); ++j) {
            pids[j] = factory::addPoint(paths[i][2*j], paths[i][2*j+1], 0.0);
        }
        std::vector<int> ls(pids.size());
        for (int j = 0; j < ls.size(); ++j) {
            ls[j] = factory::addLine(pids[j], pids[(j+1)%pids.size()]);
        }
        const int clid = factory::addCurveLoop(ls);
        const int psid = factory::addPlaneSurface({clid});
        tilingDimTags[i] = std::make_pair(2, psid);
    }

    std::vector<std::pair<int, int>> surface;
    std::vector<std::vector<std::pair<int, int>>> ovv1;
    factory::cut({{2, 1}}, tilingDimTags, surface, ovv1);

    std::vector<std::pair<int, int>> volume;
    factory::extrude(surface, 0, 0, height, volume);

    factory::synchronize();

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::setOrder(2);

    getQuadTets(vertices, tets);

    if (preview) {
        gmsh::fltk::run();
    }
    gmsh::finalize();
}

void MeshLib::structuredSheet(
    int IH,
    const std::vector<double>& params,
    double inflateWidth,
    double height,
    Tet10 &tets, 
    Eigen::SparseMatrix<double>& proj, 
    Eigen::VectorXd& reduced,
    std::array<double, 4>& trans,
    bool preview)
{
    const std::array<int, 4> nTiles = {-2, -2, 2, 2};
    std::vector<std::vector<double>> tiles = Isohedral::getTiling(IH, params, nTiles, trans, 0.5);
    // set unit to m
    for (double &t : trans) {
        t *= 1e-2;
    }
    for (auto &p : tiles) {
        for (double &v : p) {
            v *= 1e-2;
        }
    }

    std::vector<std::vector<double>> paths = Isohedral::inflate(tiles, inflateWidth); // unit in cm


    gmsh::initialize();
    gmsh::model::add("structured sheet");
    gmsh::option::setNumber("General.Verbosity", 2.0);
    namespace factory = gmsh::model::occ;

    // periodic region
    const double &dx1 = trans[0], &dy1 = trans[1], &dx2 = trans[2], &dy2 = trans[3];
    factory::addPoint(0.0, 0.0, 0.0, 0.0, 1);
    factory::addPoint(dx1, dy1, 0.0, 0.0, 2);
    factory::addPoint(dx1 + dx2, dy1 + dy2, 0.0, 0.0, 3);
    factory::addPoint(dx2, dy2, 0.0, 0.0, 4);
    factory::addLine(1, 2, 1);
    factory::addLine(2, 3, 2);
    factory::addLine(3, 4, 3);
    factory::addLine(4, 1, 4);
    factory::addCurveLoop({1, 2, 3, 4}, 1);
    factory::addPlaneSurface({1}, 1);

    // tiling polygons
    std::vector<std::pair<int, int>> tilingDimTags(paths.size());
    for (int i = 0; i < paths.size(); ++i) {
        std::vector<int> pids(paths[i].size()/2);
        for (int j = 0; j < pids.size(); ++j) {
            pids[j] = factory::addPoint(paths[i][2*j], paths[i][2*j+1], 0.0);
        }
        std::vector<int> ls(pids.size());
        for (int j = 0; j < ls.size(); ++j) {
            ls[j] = factory::addLine(pids[j], pids[(j+1)%pids.size()]);
        }
        const int clid = factory::addCurveLoop(ls);
        const int psid = factory::addPlaneSurface({clid});
        tilingDimTags[i] = std::make_pair(2, psid);
    }

    std::vector<std::pair<int, int>> surface;
    std::vector<std::vector<std::pair<int, int>>> ovv1;
    factory::cut({{2, 1}}, tilingDimTags, surface, ovv1);

    std::vector<std::pair<int, int>> volume;
    factory::extrude(surface, 0, 0, height, volume);

    factory::synchronize();

    // extrude with subdivision will cause problem if setPeriodic is before generate
    const std::vector<std::pair<int, int>> sDimTags = setPeriodic({{dx1, dy1}, {dx2, dy2}});

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::setOrder(2);

    Eigen::MatrixX3d vertices;
    getQuadTets(vertices, tets);

    getReduced(sDimTags, vertices, proj, reduced);
    
    if (preview) {
        gmsh::fltk::run();
    }
    
    gmsh::finalize();
}

std::vector<std::pair<int,int>> MeshLib::setPeriodic(const std::vector<std::pair<double, double>>& offset)
{
    // get entities on boundary
    std::vector<std::pair<int, int>> meshes, boundaryDimtags;
    gmsh::model::getEntities(meshes, 3);
    gmsh::model::getBoundary(meshes, boundaryDimtags, true, false);

    constexpr double eps = 1e-5;
    std::vector<std::pair<int, int>> sDimTags;
    for (const auto& i : boundaryDimtags) {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        gmsh::model::getBoundingBox(i.first, i.second, xmin, ymin, zmin, xmax, ymax, zmax);

        // find potential pairs
        for (const auto [dx, dy] : offset) {
            std::vector<std::pair<int, int>> js;
            gmsh::model::getEntitiesInBoundingBox(xmin + dx - eps, ymin + dy - eps, zmin - eps, 
                xmax + dx + eps, ymax + dy + eps, zmax + eps, js, i.first);

            for (const auto& j : js) {
                double xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
                gmsh::model::getBoundingBox(j.first, j.second, xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);
                xmin2 -= dx, xmax2 -= dx, ymin2 -= dy, ymax2 -= dy;
                if(std::abs(xmin2 - xmin) < eps && std::abs(xmax2 - xmax) < eps &&
                    std::abs(ymin2 - ymin) < eps && std::abs(ymax2 - ymax) < eps &&
                    std::abs(zmin2 - zmin) < eps && std::abs(zmax2 - zmax) < eps) {

                    // only getPeriodicNodes gives nodes with correspondence
                    sDimTags.push_back(std::make_pair(j.first, j.second));
                    std::vector<double> translation({1, 0, 0, dx, 0, 1, 0, dy, 0, 0, 1, 0, 0, 0, 0, 1});
                    gmsh::model::mesh::setPeriodic(i.first, {j.second}, {i.second}, translation);
                }
            }
        }
    }

    if (sDimTags.empty()) {
        throw std::runtime_error("not found periodic");
    }

    return sDimTags;
}

void MeshLib::getReduced(
    const std::vector<std::pair<int, int>>& sDimTags, 
    const Eigen::MatrixXd& vertices, 
    Eigen::SparseMatrix<double>& proj, 
    Eigen::VectorXd& reduced)
{
    // get correspondence of periodic nodes
    std::vector<std::pair<double, double>> transforms;
    std::vector<size_t> nodeTags, nodeTagsMaster;
    const int nNodes = vertices.rows();
    std::vector<bool> unvisited(nNodes, true);
    for (const auto& dimTag : sDimTags) {
        int tagMaster;
        std::vector<double> affineTransform;
        std::vector<size_t> _nodeTags, _nodeTagsMaster;
        gmsh::model::mesh::getPeriodicNodes(dimTag.first, dimTag.second, tagMaster, _nodeTags, _nodeTagsMaster, affineTransform, true);
        // there could be duplicate nodes if it is not cut
        for (int i = 0; i < _nodeTags.size(); ++i) {
            // index starts from 1 in gmsh
            if (unvisited[_nodeTags[i]-1]) {
                nodeTags.push_back(_nodeTags[i]-1);
                nodeTagsMaster.push_back(_nodeTagsMaster[i]-1);
                transforms.push_back(std::make_pair(affineTransform[3], affineTransform[7]));
            }
        }
    }

    
    // get correspondence between full and reduced indices
    Eigen::VectorXi V = Eigen::VectorXi::LinSpaced(nNodes, 0, nNodes-1);
    for (size_t i = 0; i < nodeTags.size(); ++i) {
        V(nodeTags[i]) = static_cast<int>(nodeTagsMaster[i]);
    }
    // for corner point
    for (int i = 0; i < V.size(); ++i) {
        if (V(i) != V(V(i))) {
            V(i) = V(V(i));
        }
    }
    Eigen::VectorXf SV;
    Eigen::VectorXi SVI, SVJ;
    igl::remove_duplicate_vertices(V.cast<float>(), 0, SV, SVI, SVJ);
    // there seems to be a bug in libigl
    for (int i = 0; i < SVI.size(); ++i) {
        SVI(i) = V(SVI(i));
    }

    // reduced coordinates
    reduced.resize(3 * SV.size() + 3);
    for (int i = 0; i < SV.size(); ++i) {
        reduced.segment(3*i, 3) = vertices.row(SVI(i)).transpose();
    }
    reduced.tail(3) << 1.0, 0.0, 1.0;

    // projection matrix
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(3 * nNodes + 4 * nodeTags.size());
    for (int i = 0; i < nNodes; ++i) {
        for (int d = 0; d < 3; ++d) {
            triplets.emplace_back(3*i+d, 3*SVJ(i)+d, 1.0);
        }
    }
    const int pbc_start = reduced.size() - 3;
    for (size_t i = 0; i < nodeTags.size(); ++i) {
        triplets.emplace_back(3*nodeTags[i], pbc_start, transforms[i].first);
        triplets.emplace_back(3*nodeTags[i], pbc_start+1, transforms[i].second);
        triplets.emplace_back(3*nodeTags[i], pbc_start+2, 0.0);
        triplets.emplace_back(3*nodeTags[i]+1, pbc_start, 0.0);
        triplets.emplace_back(3*nodeTags[i]+1, pbc_start+1, transforms[i].first);
        triplets.emplace_back(3*nodeTags[i]+1, pbc_start+2, transforms[i].second);
    }

    proj.resize(vertices.size(), reduced.size());
    proj.setFromTriplets(triplets.begin(), triplets.end());

    // check
    const Eigen::MatrixXd err = (proj * reduced).reshaped<Eigen::RowMajor>(vertices.rows(), vertices.cols()) - vertices;
    if (err.norm() > 1e-8) {
        throw std::runtime_error("getReduced fails?");
    }
}
