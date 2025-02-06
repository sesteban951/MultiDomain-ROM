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

        // perform open loop rollouts
        MC_Result monte_carlo(double t_sim, Vector_12d x0_sys, Vector_6d p0_feet, Domain d0);

        // perform sampling predictive control
        RHC_Result sampling_predictive_control(double t, Vector_12d x0_sys, Vector_6d p0_feet, Domain d0);

        // internal dynamics object
        Dynamics dynamics;

        // Control parameters
        ControlParams params;

        // distribution parameters
        GaussianDistribution dist;

    private:

        // initialize the cost matrices
        void initialize_costs(YAML::Node config_file);

        // initialize the global time reference trajecotry
        void initialize_reference_trajectories(YAML::Node config_file);

        // to initialize the parametric distribution
        void initialize_distribution(YAML::Node config_file);

        // initialize the class variables
        void initialize_variables();

        // sample the input trajectories from the distribution
        void sample_input_trajectory();

        // update the distribution parameters from a bundle of control inputs
        void update_distribution_params(const Vector_6d_Traj_Bundle& U_elite);

        // generate a reference trajectory for the predictive control to track
        void generate_reference_trajectory(double t_sim, const Vector_6d& x0_com);

        // evaluate log barrier function cost on legs
        double cost_limits(const Vector_12d& x_leg);

        // evaluate the cost function given a solution
        double cost_function(const ReferenceLocal& ref, const Solution& Sol, const Vector_6d_Traj& U);

        // select solutions based on cost
        void sort_trajectories(const Solution_Bundle&  S,      const Vector_6d_Traj_Bundle& U,      const Vector_1d_List& J,
                                     Solution_Bundle&  S_elite,      Vector_6d_Traj_Bundle& U_elite,      Vector_1d_List& J_elite);

        // reference types
        ReferenceGlobal ref_sim;    // global reference (throughout the entire simulation)
        ReferenceLocal ref_horizon; // local reference (for the horizon)

        // random number generator
        std::mt19937 rand_generator;
        std::normal_distribution<double> normal_dist;

        // minimum covariance norm (theoretical)
        double min_cov_norm;

        // if multi threading is enabled
        bool threading_enabled;

        // info 
        bool verbose;

        // sampling sectors and matrices
        Matrix_d L;              // cholesky factor
        Vector_d Z_vec_sample;   // standard normal vector
        Vector_d U_vec_sample;   // multi-varaite gaussian vector
        Vector_6d_Traj U_traj_sample;          // sampled input trajectory
        Vector_6d_Traj_Bundle U_traj_samples;  // sampled input trajectories

        // for updating the distribution parameters
        Vector_d mu;                 // mean
        Matrix_d Sigma;              // covariance
        Vector_6d_Traj U_elite_traj; // elite input trajectory
        Vector_d U_elite_vec;        // elite input vector
        Matrix_d U_elite_matrix;     // elite input matrix
        Vector_d eigval;           // eigenvalues of covariance matrix
        Matrix_d eigvec;           // eigenvectors of covariance matrix
        Matrix_d eigvec_inv;       // inverse of eigenvectors
        Matrix_d I;  // identity matrix

};
