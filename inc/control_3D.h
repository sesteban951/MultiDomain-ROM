#pragma once

// standard includes
#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>

// package includes
#include "yaml-cpp/yaml.h"

// custom includes
#include "types_3D.h"
#include "dynamics_3D.h"

// class for my controller
class Controller
{
    public:
        
        // Constructor and Destructor
        Controller(YAML::Node config_file);
        ~Controller(){};

        // initialize the cost matrices
        void initialize_costs(YAML::Node config_file);

        // initialize the global time reference trajecotry
        void initialize_reference_trajectories(YAML::Node config_file);

        // to initialize the parametric distribution
        void initialize_distribution(YAML::Node config_file);

        // // sample the input trajectories from the distribution
        // Vector_4d_Traj_Bundle sample_input_trajectory(int K);

        // // update the distribution parameters from a bundle of control inputs
        // void update_distribution_params(const Vector_4d_Traj_Bundle& U_bundle);

        // // generate a reference trajectory for the predictive control to track
        // Reference generate_reference_trajectory(double t_sim, Vector_4d x0_com);

        // // evaluate log barrier function cost on legs
        // double cost_log_barrier(Vector_8d x_leg);

        // // evaluate the cost function given a solution
        // double cost_function(const Reference& ref, const Solution& Sol, const Vector_4d_Traj& U);

        // // perform open loop rollouts
        // MC_Result monte_carlo(double t_sim, Vector_8d x0_sys, Vector_4d p0_feet, Domain d0, int K);

        // // select solutions based on cost
        // void sort_trajectories(const Solution_Bundle&  S,       const Vector_4d_Traj_Bundle& U,        const Vector_1d_List& J,
        //                        Solution_Bundle& S_elite, Vector_4d_Traj_Bundle& U_elite, Vector_1d_List& J_elite);

        // // perform sampling predictive control
        // RHC_Result sampling_predictive_control(double t, Vector_8d x0_sys, Vector_4d p0_feet, Domain d0);

        // internal dynamics object
        Dynamics dynamics;

        // Control parameters
        ControlParams params;

        // distribution parameters
        GaussianDistribution dist;

        // reference
        Reference ref;

        // random number generator
        std::mt19937 rand_generator;
        std::normal_distribution<double> normal_dist;

        // minimum covariance norm (theoretical)
        double min_cov_norm;

        // if multi threading is enabled
        bool threading_enabled;

        // info 
        bool verbose;
};
