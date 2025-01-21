// standard includes
#include <iostream>
#include <chrono>
#include <fstream>

// package includes
#include "Eigen/Dense"
#include "yaml-cpp/yaml.h"

// custom includes
#include "../inc/types.h"
#include "../inc/dynamics.h"
#include "../inc/control.h"

int main()
{
    // load parameters from yaml file
    std::string config_file_path = "../config/config.yaml";   
    YAML::Node config_file = YAML::LoadFile(config_file_path);
    
    // create dynamics object
    Dynamics dynamics(config_file);

    // initial conditions
    Vector_8d x0;
    Vector_2d p0_left, p0_right;
    Vector_2d_List p0_feet(2);
    Domain d0(2, Contact::SWING);

    std::vector<double> x0_temp = config_file["STATE"]["x0"].as<std::vector<double>>();
    x0 << x0_temp[0],    // px_com
          x0_temp[1],    // pz_com
          x0_temp[2],    // vx_com
          x0_temp[3],    // vz_com
          x0_temp[4],    // l0_command (left)
          x0_temp[5],    // theta_command (left)
          x0_temp[6],    // l0_command (right)
          x0_temp[7];    // theta_command (right)
    p0_left << -0.1, 0;
    p0_right << 0.1, 0;
    p0_feet[0] = p0_left;
    p0_feet[1] = p0_right;
    Vector_2d_List u(2);
    u[0] << 0.0, 0.0;
    u[1] << 0.0, 0.0;

    // query the dynamics
    Vector_8d xdot = dynamics.dynamics(x0, u, p0_feet, d0);

    std::cout << "xdot: " << xdot.transpose() << std::endl;

    


    return 0;
}
