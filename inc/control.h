#ifndef CONTROL_H
#define CONTROL_H

// standard includes
#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>

// package includes
#include "yaml-cpp/yaml.h"

// custom includes
#include "types.h"
#include "dynamics.h"

// class for my controller
class Controller
{
    public:
        
        // Constructor and Destructor
        Controller(YAML::Node config_file);
        ~Controller(){};

        // to initialize the initial distribution
        void initialize_distribution(YAML::Node config_file);

        // sample the input trajectories from the distribution
        Vector_2d_Traj_Bundle sample_input_trajectory(int K);

        // update the distribution parameters from a bundle of control inputs
        void update_distribution_params(Vector_2d_Traj_Bundle U_bundle);

        // internal dynamics object
        Dynamics dynamics;

        // Control parameters
        ControlParams params;

        // distribution parameters
        GaussianDistribution dist;

        // random number generator
        std::mt19937 rand_generator;
        std::normal_distribution<double> normal_dist;

        // desired reference trajectory
        double pz_des;
        double vx_des;
        double r_des;
};

#endif