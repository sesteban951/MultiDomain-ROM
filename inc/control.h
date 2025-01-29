#pragma once

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

        // construct the long horizon reference trajecotry
        void initialize_reference_trajectories(YAML::Node config_file);

        // to initialize the initial distribution
        void initialize_distribution(YAML::Node config_file);

        // sample the input trajectories from the distribution
        Vector_4d_Traj_Bundle sample_input_trajectory(int K);

        // update the distribution parameters from a bundle of control inputs
        void update_distribution_params(Vector_4d_Traj_Bundle U_bundle);

        // generate a reference trajectory for the predictive control to track
        Reference generate_reference_trajectory(double t_sim, Vector_4d x0_com);

        // evaluate log barrier function cost on legs
        double cost_log_barrier(Vector_8d x_leg);

        // evaluate the cost function given a solution
        double cost_function(Reference ref, Solution Sol, Vector_4d_Traj U);

        // perform open loop rollouts
        MC_Result monte_carlo(double t_sim, Vector_8d x0_sys, Vector_4d p0_feet, Domain d0, int K);

        // select solutions based on cost
        void sort_trajectories(const Solution_Bundle&  S,       const Vector_4d_Traj_Bundle& U,        const Vector_1d_List& J,
                               Solution_Bundle& S_elite, Vector_4d_Traj_Bundle& U_elite, Vector_1d_List& J_elite);

        // perform sampling predictive control
        RHC_Result sampling_predictive_control(double t, Vector_8d x0_sys, Vector_4d p0_feet, Domain d0);

        // internal dynamics object
        Dynamics dynamics;

        // Control parameters
        ControlParams params;

        // distribution parameters
        GaussianDistribution dist;

        // random number generator
        std::mt19937 rand_generator;
        std::normal_distribution<double> normal_dist;

        // if multi threading is enabled
        bool threading_enabled;

        // desired reference trajectory
        double pz_des;
        double vx_des;
        double r_des;
        double theta_des;
        double T_cycle;
        double T_SSP;
        Vector_1d_Traj t_ref;
        Vector_2d_Traj p_com_ref;

        // minimum covariance norm (theoretical)
        double min_cov_norm;

        // info 
        bool verbose;

};
