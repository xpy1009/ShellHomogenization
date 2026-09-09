#include <Mesh/MeshLib.h>
#include <Element/Shell.h>
#include <Compression.h>

#include <torch/script.h>
#include <igl/writeOBJ.h>

int main()
{
    const auto start = std::chrono::high_resolution_clock::now();

    const torch::jit::Module stretchModel(torch::jit::load("../data/stretch.pt"));
    const torch::jit::Module bendModel(torch::jit::load("../data/bend.pt"));

    constexpr double dx = 0.0930605;
    constexpr std::array<double, 3> dys = {0.0322371, 0.0537285, 0.0752199};
    constexpr bool verbose = false;
    for (int i = 0; i < 3; ++i) {
        const int n = 3 + 2 * i;
        const double& dy = dys[i];
        Eigen::MatrixX2d restPos;
        Eigen::MatrixX3i F;
        MeshLib::rectangle(dx, dy, restPos, F);
        Shell mesh(F);

        std::ofstream file("../data/H" + std::to_string(n) + ".txt");
        const Eigen::MatrixX3d pos = Compression::run(mesh, restPos, stretchModel, bendModel, file, verbose);
        igl::writeOBJ("../data/H" + std::to_string(n) + ".obj", pos, F);
    }

    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> diff = end - start;
    std::cout << diff << std::endl;

    return 0;
}