## Dependencies
### yaml-cpp
```bash
https://github.com/jbeder/yaml-cpp
```
### Eigen 3.4.0
```bash
https://eigen.tuxfamily.org/index.php?title=Main_Page
```
See the ```CMakeLists.txt``` and modify it if necessary. 

## Running
Build the executables:
```bash
cd <path_to_MultiDomain-ROM>
mkdir build
cmake ..
make
```

Create a data folder in teh root directory to store some simulation data:
```bash
cd ..
mkdir data
```

Change values in ```config/config.yaml``` (2D case) or ```config/config_3D.yaml``` (3D case) to play with the controller.

To run the the simulation:
```bash
cd build
./a    # 2D case
./s    # 3D case
```

Plot the result in MATLAB by running the ```plot_data.m``` or ```plot_data_3D.m``` MATLAB script.

