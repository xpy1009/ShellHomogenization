#include <Mesh/Isohedral.h>

#include <tiling.hpp>
#include <glm/glm.hpp>

#include <clipper2/clipper.h>

std::vector<std::vector<double>> Isohedral::getTiling(
    int IH, 
    const std::vector<double>& params, 
    const std::array<int, 4>& nTiles, 
    std::array<double, 4>& trans,
	double normalizePos,
	bool normalizeArea)
{
    // Construct a tiling of the given type.
    csk::IsohedralTiling t( IH );
    // Create an array to hold a copy of the tiling vertex parameters.
	double ps[ t.numParameters() ];
	// Now fill the array with the current values of the parameters,
	// which will be set to reasonable defaults when the tiling is
	// created.
	t.getParameters( ps );
    // Perturb the parameters a bit to get a novel tiling.
	if (!params.empty() && params.size() != t.numParameters()) {
		throw std::runtime_error("isohedral parameters not match");
	}
	for (size_t i = 0; i < params.size(); ++i) {
		ps[i] = params[i];
	}
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
	// Note: simplify every edge to a single line segment
	const std::vector<std::vector<glm::dvec2>> edges(static_cast<size_t>(t.numEdgeShapes()), 
													 std::vector<glm::dvec2>({glm::dvec2(0, 0), glm::dvec2(1, 0)}));


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

    // periodic translation
    const glm::dvec2& t1 = t.getT1();
    const glm::dvec2& t2 = t.getT2();
    trans = { t1.x, t1.y, t2.x, t2.y }; 

	// normalize to unit area
	const double normalizer = normalizeArea ? sqrt(abs(t1.x * t2.y - t1.y * t2.x)) : 1.0;
	for (double &t : trans) {
		t /= normalizer;
	}

	// shift to shape center
	glm::dvec2 shapeCenter(0.0, 0.0);
	for (const glm::dvec2& p : shape) {
		shapeCenter += p;
	}
	shapeCenter *= normalizePos / shape.size();

	
    const double dx = std::max(abs(t1.x), abs(t2.x));
    const double dy = std::max(abs(t1.y), abs(t2.y));

    std::vector<std::vector<double>> tiles;
    // Ask the tiling to generate (approximately) enough tiles to
	// fill the bounding box below.  The bounding box is a bit bigger
	// than the box we actually want to display in the document, to
	// hopefully ensure that it completely covers that box.
	for( auto i : t.fillRegion( nTiles[0]*dx, nTiles[1]*dy, nTiles[2]*dx, nTiles[3]*dy ) ) {
		// The region filling algorithm will give us a transform matrix
		// that takes a tile in default position to its location in the
		// tiling.
		glm::dmat3 T = i->getTransform();

        std::vector<double> list(2 * shape.size());
        for(size_t j = 0; j < shape.size(); ++j) {
            const glm::dvec2 p = T * glm::dvec3(shape[j], 1.0);
            list[2*j] = (p.x - shapeCenter.x) / normalizer;
            list[2*j+1] = (p.y - shapeCenter.y) / normalizer;
        }
        tiles.push_back(list);
    }
    return tiles;
}


std::vector<std::vector<double>> Isohedral::inflate(
	const std::vector<std::vector<double>>& tiles,
	double width)
{
    // mult is used to increase precision of Clipper2, necessary for gmsh setPeriodic
    constexpr double mult = 1e8;
    Clipper2Lib::PathsD polygon(tiles.size());
    for (size_t i = 0; i < tiles.size(); ++i) {
        Clipper2Lib::PathD pi = Clipper2Lib::MakePathD(tiles[i]);
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
    const double delta = -0.5 * width * mult;
    Clipper2Lib::PathsD solution = Clipper2Lib::InflatePaths(
        polygon, delta, Clipper2Lib::JoinType::Bevel, Clipper2Lib::EndType::Polygon);
    
    std::vector<std::vector<double>> paths(solution.size());
    for (size_t i = 0; i < solution.size(); ++i) {
        paths[i].resize(2 * solution[i].size());
        for (size_t j = 0; j < solution[i].size(); ++j) {
            paths[i][2*j] = solution[i][j].x / mult;
            paths[i][2*j+1] = solution[i][j].y / mult;
        }
    }

    return paths;
}