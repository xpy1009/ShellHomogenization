#include <Mesh/MeshLib.h>
#include <gmsh.h>

#include <tiling.hpp>
#include <clipper2/clipper.h>
#include <Eigen/Dense>
#include <iostream>

namespace Eigen { using namespace Eigen::placeholders; }
#include <igl/remove_duplicate_vertices.h>
#include <igl/per_face_normals.h>

constexpr double angle = 0;
constexpr double height = 0.3;

// constexpr double width = 0.1;
// constexpr int IH = 1;
// constexpr bool forceCenter = true;
// constexpr double scale = 1;

constexpr double width = 0.1;
constexpr int IH = 1;
constexpr bool forceCenter = true;
constexpr double scale = 1;

// constexpr double width = 0.1;
// constexpr int IH = 41;
// constexpr bool forceCenter = true;
// constexpr double scale = 1;

// constexpr double width = 0.1;
// constexpr int IH = 21;
// constexpr bool forceCenter = false;
// constexpr double scale = .5;

// constexpr double width = 0.1;
// constexpr int IH = 27;
// constexpr bool forceCenter = false;
// constexpr double scale = .7;

// constexpr double width = 0.1;
// // constexpr double angle = 0;
// // constexpr double height = 0.2;
// constexpr int IH = 42;
// constexpr bool forceCenter = false;
// constexpr double scale = 1;

void MeshLib::getQuadTets(Eigen::MatrixXd &vertices, Eigen::MatrixXi &tets)
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

void MeshLib::generateBox(Eigen::MatrixXd &vertices, Eigen::MatrixXi &tets)
{
    gmsh::initialize();
    gmsh::model::add("tets");
    namespace factory = gmsh::model::occ;

    factory::addBox(0,0,0, .3,.3,.005);

    factory::synchronize();

    // gmsh::option::setNumber("Mesh.MeshSizeFactor", .7);
    // gmsh::model::mesh::setTransfiniteAutomatic({}, 2.35, false);
    // gmsh::option::setNumber("Mesh.Algorithm3D", 4);

    // silent gmsh
    gmsh::option::setNumber("General.Verbosity", 2.0);
    // 3d
    gmsh::model::mesh::generate(3);
    // quadratic
    gmsh::model::mesh::setOrder(2);

    getQuadTets(vertices, tets);

    gmsh::fltk::run();
    gmsh::finalize();
}

void MeshLib::rectangle(Eigen::MatrixX3d &vertices, Eigen::MatrixX3i &faces)
{
    gmsh::initialize();
    gmsh::model::add("tris");
    namespace factory = gmsh::model::occ;

    factory::addRectangle(0,0,0, 1,1);

    factory::synchronize();

    // gmsh::option::setNumber("Mesh.MeshSizeFactor", 10);
    
    gmsh::option::setNumber("General.Verbosity", 2.0); // silent gmsh
    gmsh::model::mesh::setTransfiniteAutomatic({}, 2.35, false); // regular
    gmsh::model::mesh::generate(2); // 2d

    // vertices
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);
    vertices = Eigen::Map<Eigen::MatrixXd>(nodeCoords.data(), 3, nodeCoords.size()/3).transpose();

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

    gmsh::fltk::run();
    gmsh::finalize();
}

void MeshLib::constructTiling(int idx, Eigen::Matrix2d& bbox, std::vector<std::vector<double>>& tiles)
{
    // Construct a tiling of the given type.
    csk::IsohedralTiling t( idx );
    // csk::IsohedralTiling t( csk::tiling_types[ idx ] );
    // Create an array to hold a copy of the tiling vertex parameters.
	double ps[ t.numParameters() ];
	// Now fill the array with the current values of the parameters,
	// which will be set to reasonable defaults when the tiling is
	// created.
	t.getParameters( ps );
    // Perturb the parameters a bit to get a novel tiling.
    // ...
    // ps[0] = 0.25, ps[2] = 0.07; // IH01 negative Poisson ratio
    // ps[1] = ps[2]=0.1;
    // Now send those parameters back to the tiling.
	t.setParameters( ps );

    // Create a vector to hold some edge shapes.  The tiling tells you
	// how many distinct edge shapes you need, but doesn't know anything
	// about how those shapes might be represented.  It simply assumes
	// that each one will be a curve from (0,0) to (1,0).  The tiling
	// provides tools to let you map those curves into position around
	// the outline of a tile.  All the curves below have exactly four
	// control points, so using a vector is overkill; but it offers a 
	// more convenient starting point for experimentation with fancier
	// curves, so I'll keep it.
	std::vector<glm::dvec2> edges[ t.numEdgeShapes() ];

    // Generate some random edge shapes.
	for( csk::U8 idx = 0; idx < t.numEdgeShapes(); ++idx ) {
		std::vector<glm::dvec2> ej;

		// Start by making a random Bezier segment.
        // Note: change to interpolate as linear function
		ej = {glm::dvec2(0, 0), glm::dvec2(1, 0)};

		// Now, depending on the edge shape class, enforce symmetry 
		// constraints on edges.
		switch( t.getEdgeShape( idx ) ) {
		case csk::J: 
			break;
		case csk::U:
			ej[2].x = 1.0 - ej[1].x;
			ej[2].y = ej[1].y;
			break;
		case csk::S:
			ej[2].x = 1.0 - ej[1].x;
			ej[2].y = -ej[1].y;
			break;
		case csk::I:
			ej[1].y = 0.0;
			ej[2].y = 0.0;
			break;
		}
		edges[idx] = ej;
	}

    // Use a vector to hold the control points of the final tile outline.
	std::vector<glm::dvec2> shape;

    // Iterate over the edges of a single tile, asking the tiling to
	// tell you about the geometric information needed to transform 
	// the edge shapes into position.  Note that this iteration is over
	// whole tiling edges.  It's also to iterator over partial edges
	// (i.e., halves of U and S edges) using t.parts() instead of t.shape().
	for( auto i : t.shape() ) {
		// Get the relevant edge shape created above using i->getId().
		const std::vector<glm::dvec2>& ed = edges[ i->getId() ];
		// Also get the transform that maps to the line joining consecutive
		// tiling vertices.
		const glm::dmat3& T = i->getTransform();

		// If i->isReversed() is true, we need to run the parameterization
		// of the path backwards.
		if( i->isReversed() ) {
			for( size_t idx = 1; idx < ed.size(); ++idx ) {
				shape.push_back( T * glm::dvec3( ed[ed.size()-1-idx], 1.0 ) );
			}
		} else {
			for( size_t idx = 1; idx < ed.size(); ++idx ) {
				shape.push_back( T * glm::dvec3( ed[idx], 1.0 ) );
			}
		}
	}

    double x1 = bbox(0,0), y1 = bbox(0,1), x2 = bbox(1,0), y2 = bbox(1,1);
    // repeating pattern
    glm::dvec2 t1 = t.getT1();
    glm::dvec2 t2 = t.getT2();
    bbox << t1.x,t2.x,t1.y,t2.y; 

    // normalize with square root of area
    const double n = std::sqrt(abs(bbox.determinant())) * scale;
    bbox /= n;

    // shift with center at origin
    glm::dvec2 center(0,0);
    if (forceCenter) {
        for (const glm::dvec2& p : shape) {
            center += p;
        }
        center /= shape.size();
        // for (glm::dvec2& p : shape) {
        //     p -= center;
        // }
    }
    Eigen::Matrix2d R;
    R << cos(angle), sin(angle), -sin(angle), cos(angle);
    bbox = R * bbox;

    // Note: not sure how to use fillRegion. Use a large enough region for now. There should be an elegant way.
    
    // Ask the tiling to generate (approximately) enough tiles to
	// fill the bounding box below.  The bounding box is a bit bigger
	// than the box we actually want to display in the document, to
	// hopefully ensure that it completely covers that box.
	for( auto i : t.fillRegion( x1, y1, x2, y2 ) ) {
	// for( auto i : t.fillRegion( -5, -15, 35, 15 ) ) {
	// for( auto i : t.fillRegion( -15, -2, 15, 35 ) ) {
		// The region filling algorithm will give us a transform matrix
		// that takes a tile in default position to its location in the
		// tiling.
		glm::dmat3 T = i->getTransform();

        std::vector<double> list(2 * shape.size());
        for(size_t j = 0; j < shape.size(); ++j) {
            glm::dvec2 p = T * glm::dvec3(shape[j], 1.0);
            p -= center;
            list[2*j] = (cos(angle) * p.x + sin(angle) * p.y) / n;
            list[2*j+1] = (-sin(angle) * p.x + cos(angle) * p.y) / n;
        }
        tiles.push_back(list);
    }
    
}

std::vector<std::vector<double>> MeshLib::inflatePaths(const std::vector<std::vector<double>>& lists)
{
    // mult is used to increase precision of Clipper2, necessary for gmsh setPeriodic
    constexpr float mult = 1e8;
    Clipper2Lib::PathsD polygon(lists.size());
    for (size_t i = 0; i < lists.size(); ++i)
    {
        Clipper2Lib::PathD pi = Clipper2Lib::MakePathD(lists[i]);
        for (auto& vi : pi) {
            vi.x *= mult;
            vi.y *= mult;
        }
        // Might need to check orientation
        if (Clipper2Lib::IsPositive(pi)) {
            std::reverse(pi.begin(), pi.end());
        }
        polygon[i] = pi;
    }
    constexpr float delta = -width * mult;
    Clipper2Lib::PathsD solution = Clipper2Lib::InflatePaths(
        polygon, delta, Clipper2Lib::JoinType::Bevel, Clipper2Lib::EndType::Polygon);
    
    std::vector<std::vector<double>> paths(solution.size());
    for (size_t i = 0; i < solution.size(); ++i)
    {
        paths[i].resize(2 * solution[i].size());
        for (size_t j = 0; j < solution[i].size(); ++j)
        {
            paths[i][2*j] = solution[i][j].x / mult;
            paths[i][2*j+1] = solution[i][j].y / mult;
        }
    }

    return paths;
}

void MeshLib::structuredSheet(Eigen::MatrixXd &vertices, Eigen::MatrixXi &tets)
{
    Eigen::Matrix2d bbox; 
    std::vector<double> tile;
    std::vector<std::vector<double>> tiles;
    bbox << -5,-5,15, 15;
    constructTiling(IH, bbox, tiles);
    std::cout << "bbox: " << std::endl << bbox << std::endl;

    std::vector<std::vector<double>> paths = inflatePaths(tiles);

    gmsh::initialize();
    gmsh::model::add("test");
    namespace factory = gmsh::model::occ;

    // square sheet
    int nx = 6;
    int ny = 5;
    double dx = abs(nx * bbox(0,1));
    double dy = abs(ny * bbox(1,0));
    std::cout << nx << " " << ny << std::endl;
    std::cout << dx << " " << dy << std::endl;
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

    gmsh::option::setNumber("General.Verbosity", 2.0);

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::setOrder(2);

    getQuadTets(vertices, tets);
    vertices.col(2).array() -= 0.5 * height;
    vertices *= 1e-2;

    gmsh::fltk::run();
    gmsh::finalize();
    // std::exit(0);
}

void MeshLib::structuredSheetPeriodic(
    Eigen::MatrixXi &tets, Eigen::SparseMatrix<double>& paramsToCoords, Eigen::VectorXd& params, Eigen::Matrix2d& bbox, 
    std::vector<int>* bdyIdx1, std::vector<int>* bdyIdx2)
{
    std::vector<std::vector<double>> tiles;
    bbox << -5,-5,5,5;
    constructTiling(IH, bbox, tiles);
    std::cout << "bbox: " << std::endl << bbox << std::endl;

    std::vector<std::vector<double>> paths = inflatePaths(tiles);

    gmsh::initialize();
    gmsh::model::add("test");
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
    // gmsh::option::setNumber("Mesh.MeshSizeFactor", 2);

    // get surface
    std::vector<std::pair<int, int>> meshes, boundaryDimtags;
    gmsh::model::getEntities(meshes, 3);
    gmsh::model::getBoundary(meshes, boundaryDimtags, true, false);
    const std::vector<std::pair<int, int>> bdyDimTags1 = setPeriodic(boundaryDimtags, dx1, dy1);
    const std::vector<std::pair<int, int>> bdyDimTags2 = setPeriodic(boundaryDimtags, dx2, dy2);

    gmsh::option::setNumber("General.Verbosity", 2.0);

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::setOrder(2);

    Eigen::MatrixXd vertices;
    getQuadTets(vertices, tets);
    // move midsurface to z=0 for moment computation
    vertices.col(2).array() -= 0.5 * height;

    // build reduced coordinates
    buildReduced(vertices, bdyDimTags1, bdyDimTags2, paramsToCoords, params, bdyIdx1, bdyIdx2);

    gmsh::fltk::run();
    gmsh::finalize();
    // std::exit(0);
}

void MeshLib::buildReduced(
    const Eigen::MatrixXd& vertices, const std::vector<std::pair<int, int>>& dimTags1, const std::vector<std::pair<int, int>>& dimTags2,
    Eigen::SparseMatrix<double>& reducedToFull, Eigen::VectorXd& reduced, 
    std::vector<int>* bdyIdx1, std::vector<int>* bdyIdx2)
{
    // get correspondence of periodic nodes
    std::vector<std::pair<double, double>> transforms;
    std::vector<size_t> nodeTags, nodeTagsMaster;
    auto getCorrespondene = [&](const std::vector<std::pair<int, int>>& dimTags, std::vector<int>* bdyIdx) {
        const int offset = nodeTags.size();
        for (const auto& dimTag : dimTags) {
            int tagMaster;
            std::vector<double> affineTransform;
            std::vector<size_t> _nodeTags, _nodeTagsMaster;
            gmsh::model::mesh::getPeriodicNodes(dimTag.first, dimTag.second, tagMaster, _nodeTags, _nodeTagsMaster, affineTransform, true);
            for (int i = 0; i < _nodeTags.size(); ++i) {
                // duplicate nodes if it is not cut
                if (std::find(std::next(nodeTags.begin(), offset), nodeTags.end(), _nodeTags[i]) == nodeTags.end()) {
                    nodeTags.push_back(_nodeTags[i]);
                    nodeTagsMaster.push_back(_nodeTagsMaster[i]);
                    transforms.push_back(std::make_pair(affineTransform[3], affineTransform[7]));
                    if (bdyIdx != nullptr) {
                        bdyIdx->push_back(static_cast<int>(_nodeTags[i] - 1));
                    }
                }
            }
        }
    };
    getCorrespondene(dimTags1, bdyIdx1);
    getCorrespondene(dimTags2, bdyIdx2);

    // index starts from 1 in gmsh
    for (size_t i = 0; i < nodeTags.size(); ++i) {
        nodeTags[i] -= 1;
        nodeTagsMaster[i] -= 1;
    }
    
    // get correspondence between full and reduced indices
    const int nNodes = vertices.rows();
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
    Eigen::VectorXd SV;
    Eigen::VectorXi SVI, SVJ;
    igl::remove_duplicate_vertices(V.cast<double>(), 0, SV, SVI, SVJ);

    // reduced coordinates
    Eigen::VectorXd full = vertices.reshaped<Eigen::RowMajor>();
    for (size_t i = 0; i < nodeTags.size(); i++) {
        full.segment(3*nodeTags[i], 3) = full.segment(3*nodeTagsMaster[i], 3);
    }
    reduced.resize(3 * SV.size() + 3);
    for (int i = 0; i < SV.size(); i++) {
        reduced.segment(3*i, 3) = full.segment(3*SVI(i), 3);
    }
    reduced.tail(3) << 1.0, 0.0, 1.0;

    // map reduced to full
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
        triplets.emplace_back(3*nodeTags[i]+1, pbc_start+1, transforms[i].first);
        triplets.emplace_back(3*nodeTags[i]+1, pbc_start+2, transforms[i].second);
    }

    reducedToFull.resize(full.size(), reduced.size());
    reducedToFull.setFromTriplets(triplets.begin(), triplets.end());

    // check
    const Eigen::MatrixXd err = (reducedToFull * reduced).reshaped<Eigen::RowMajor>(vertices.rows(), vertices.cols()) - vertices;
    if (err.norm() > 1e-8) {
        throw std::runtime_error("Build reduced fails?");
    }
}

std::vector<std::pair<int,int>> MeshLib::setPeriodic(
    const std::vector<std::pair<int, int>>& dimTags, double dx, double dy)
{
    std::vector<std::pair<int,int>> bdyDimTags;
    constexpr double eps = 1e-3;
    std::vector<double> translation({1, 0, 0, dx, 0, 1, 0, dy, 0, 0, 1, 0, 0, 0, 0, 1});
    for (const auto& i : dimTags) {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        gmsh::model::getBoundingBox(i.first, i.second, xmin, ymin, zmin, xmax, ymax, zmax);

        // find potential
        std::vector<std::pair<int, int> > js;
        gmsh::model::getEntitiesInBoundingBox(xmin + dx - eps, ymin + dy - eps, zmin - eps, 
            xmax + dx + eps, ymax + dy + eps, zmax + eps, js, i.first);

        int cnt = 0;
        for (const auto& j : js) {
            double xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
            gmsh::model::getBoundingBox(j.first, j.second, xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);
            xmin2 -= dx, xmax2 -= dx, ymin2 -= dy, ymax2 -= dy;
            if(std::abs(xmin2 - xmin) < eps && std::abs(xmax2 - xmax) < eps &&
                std::abs(ymin2 - ymin) < eps && std::abs(ymax2 - ymax) < eps &&
                std::abs(zmin2 - zmin) < eps && std::abs(zmax2 - zmax) < eps) {

                if (cnt > 0) {
                    throw std::runtime_error("Multiple periodic?");
                }
                // only getPeriodicNodes gives nodes with correspondence
                bdyDimTags.push_back(std::make_pair(j.first, j.second));
                gmsh::model::mesh::setPeriodic(i.first, {j.second}, {i.second}, translation);
                cnt++;
            }

        }

    }
    if (bdyDimTags.empty()) {
        throw std::runtime_error("No periodic?");
    }
    return bdyDimTags;
}

void MeshLib::structuredSheetMidsurface(Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces)
{
    Eigen::Matrix2d bbox; 
    std::vector<double> tile;
    std::vector<std::vector<double>> tiles;
    bbox << -5,-5,15, 15;
    constructTiling(IH, bbox, tiles);
    std::cout << "bbox: " << std::endl << bbox << std::endl;

    std::vector<std::vector<double>> paths = inflatePaths(tiles);

    gmsh::initialize();
    gmsh::model::add("test");
    namespace factory = gmsh::model::occ;

    // square sheet
    int nx = abs(floor(4 / bbox(0,1)));
    int ny = abs(floor(3 / bbox(1,0)));
    nx = 6, ny = 6;
    double dx = abs(nx * bbox(0,1));
    double dy = abs(ny * bbox(1,0));
    std::cout << dx << " " << dy << std::endl;
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

    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::pair<int, int>>> ovv1;
    factory::fragment({{2, 1}}, tilingDimTags, out, ovv1);

    factory::synchronize();

    // We then retrieve all the volumes in the bounding box of the original cube,
    // and delete all the parts outside it:
    double eps = 1e-3;
    std::vector<std::pair<int, int> > in;
    gmsh::model::getEntitiesInBoundingBox(-eps, -eps, -eps, dx + eps, dy + eps, eps, in, 2);
    for(auto i : in) {
        auto it = std::find(out.begin(), out.end(), i);
        if(it != out.end()) out.erase(it);
    }
    gmsh::model::removeEntities(out, true); // Delete outside parts recursively

    gmsh::option::setNumber("General.Verbosity", 2.0);

    gmsh::model::mesh::generate(2);

    // vertices
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);
    vertices = Eigen::Map<Eigen::MatrixXd>(nodeCoords.data(), 3, nodeCoords.size()/3).transpose();

    // triangles 
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, _elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, _elemNodeTags, 2);
    if (_elemNodeTags.size() < 1) {
        std::cout << "Not triangles" << std::endl;
        return;
    }
    const std::vector<size_t> &elemNodeTags = _elemNodeTags[0];
    constexpr int n = 3;
    const int nF = elemNodeTags.size() / n;
    faces.resize(nF, n);
    for (int i = 0; i < nF; ++i) {
         faces.row(i) << elemNodeTags[n*i], elemNodeTags[n*i+1], elemNodeTags[n*i+2];
    }
    faces.array() -= 1.0;

    Eigen::MatrixXd N;
    igl::per_face_normals(vertices, faces, N);
    for (int i = 0; i < N.rows(); ++i) {
        if (N(i, 2) < 0) {
            std::swap(faces(i, 0), faces(i, 1));
        }
    }
    vertices *= 1e-2;

    gmsh::fltk::run();
    gmsh::finalize();
    // std::exit(0);

}



