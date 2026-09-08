# Taking a Moment to Characterize the Bending Response of Thin Sheet Materials

A reference implementation for the paper [Taking a Moment to Characterize the Bending Response of Thin Sheet Materials](https://onlinelibrary.wiley.com/doi/10.1111/cgf.70513).

## Dependencies

- gmsh
- Eigen
- libigl
- tactile
- Clipper2
- LibTorch

## Building (Tested on MacOS)


Download the code with:
```
git clone --recurse-submodules https://git.ista.ac.at/wojtan-group/peiyuan-xie/boundaryblaze.git
```
Repos of tactile and libigl are automatically downloaded.

Install the dependencies using some package manager, such as Homebrew:
```
brew install gmsh clipper2 pytorch 
```
Eigen5 should be automatically installed.


Python will be needed for training, but the version should not be important. Only `pytorch` and `matplotlib` needs to be installed, which can be done via `pip`:
```
pip install torch matplotlib
```


Build using `CMake` and your favorite build system:
```
cd boundaryblaze
mkdir build 
cd build
cmake ..
make -j8
```


### Replicating Results(Timed on MacBook Air M2)
To reproduce **Fig. 8** in the paper that shows the Gaussian curvature of thick structured sheets and homogenized thin sheets with varying width, run the following commands:

```
./run_native
python ../train.py
./run_homogenized
python ../plot.py
```

- `./run_native` simulates the native-scale model and generates data for homogenization and the compression experiments described in **Section 4.2**. This step takes about 25 minutes.

- `python ../train.py` trains the neural stretching and bending models from the generated data. This step takes about 10 minutes.

- `./run_homogenized` uses the trained models to simulate a thin sheet. This step takes about 1 minute.

- `python ../plot.py` plots the evolution of Gaussian curvatures from different models under compression.

All data including the corresponding meshes are saved under the folder `data/`. 
