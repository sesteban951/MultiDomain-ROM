#ifndef DYNAMICS_H
#define DYNAMICS_H

// standard includes
#include <iostream>

// package includes
#include "yaml-cpp/yaml.h"

// custom includes
#include "types.h"

// class for system dynamics
class Dynamics
{
    public:
        
        // Constructor and Destructor
        Dynamics(YAML::Node config_file);  
        ~Dynamics(){};

        // NonLinear System dynamics, xdot = f(x, u, d)
        Vector_8d dynamics(Vector_8d x, 
                           Vector_2d_List u, 
                           Vector_2d_List p_feet,
                           Domain d);

        // System parameters
        SystemParams params;

        // fixed parameters
        const double nx_sys = 8;  // vector size of system state variables
        const double nx_leg = 4;  // vector size of control inputs
        const double nu = 4;      // vector size of control inputs
        const int n_leg = 2;      // vector size of legs
};

#endif