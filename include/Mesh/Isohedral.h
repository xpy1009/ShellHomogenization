#ifndef MESH_ISOHEDRAL_H
#define MESH_ISOHEDRAL_H

#include <vector>
#include <array>

class Isohedral
{    
public:
    Isohedral() = delete;

    // from tactile/demo/psdemo.cpp
    static std::vector<std::vector<double>> getTiling(
        int IH, 
        const std::vector<double>& params,
        const std::array<int, 4>& nTiles,
        std::array<double, 4>& trans,
        double normalizePos=1.0,
        bool normalizeArea=true);

    // inflate via Clipper2
    static std::vector<std::vector<double>> inflate(
        const std::vector<std::vector<double>>& tiles,
        double width);
};

#endif