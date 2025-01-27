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

    Controller controller(config_file);

    // initial conditions
    Vector_8d x0_sys;
    Vector_2d p0_left, p0_right;
    Vector_4d p0_feet;

    std::vector<double> x0_temp = config_file["STATE"]["x0"].as<std::vector<double>>();
    x0_sys << x0_temp[0],    // px_com
              x0_temp[1],    // pz_com
              x0_temp[2],    // vx_com
              x0_temp[3],    // vz_com
              x0_temp[4],    // l0_command (left)
              x0_temp[5],    // theta_command (left)
              x0_temp[6],    // l0dot_command (left)
              x0_temp[7];    // thetadot_command (left)
    p0_left << -0.1, 0.0;
    p0_right << 0.1, 0.0;
    p0_feet << p0_left, p0_right;
    Domain d0(2, Contact::SWING);

    ////////////////////////////////// Function testing //////////////////////////////////

    Vector_4d u_const;
    Vector_2d u_L, u_R;
    u_L << 0.0, 0.0;
    u_R << 0.0, -0.0;
    u_const << u_L, u_R;

    // example rollout of the dynamics
    int N = 140;
    Vector_1d_Traj T_x(N);
    for (int i = 0; i < N; i++) {
        T_x[i] = i * 0.01;
    }

    int Nu = 70;
    Vector_1d_Traj T_u(Nu);
    Vector_4d_Traj U(Nu);
    for (int i = 0; i < Nu; i++) {
        T_u[i] = i * 0.02;
        U[i] = u_const;
    }

    // query the dynamics
    // DynamicsResult res = dynamics.dynamics(x0_sys, u_const, p0_feet, d0);
    // Vector_8d xdot = res.xdot;
    // std::cout << "xdot: " << xdot.transpose() << std::endl;

    // // leg state
    // Vector_8d x_leg = dynamics.compute_leg_state(x0_sys, u_const, p0_feet, d0);
    // std::cout << "x_leg: " << x_leg.transpose() << std::endl;

    // // foot state
    // Vector_8d x_foot = dynamics.compute_foot_state(x0_sys, x_leg, p0_feet, d0);
    // std::cout << "x_foot: " << x_foot.transpose() << std::endl;
    
    // Solution sol = dynamics.RK3_rollout(T_x, T_u, x0_sys, p0_feet, d0, U);
    // // std::cout << "Rollout compl
    // generate inputs
    // Vector_4d_Traj_Bundle U_bundle = controller.sample_input_trajectory(1000);

    // // update the distribution parameters
    // controller.update_distribution_params(U_bundle);
    // std::cout << "mean: \n" << controller.dist.mean.transpose() << std::endl;
    // std::cout << "cov: \n" << controller.dist.cov << std::endl;

    // Do a rollout of the dynamics
    // U = U_bundle[0];
    // Solution sol = dynamics.RK3_rollout(T_x, T_u, x0_sys, p0_feet, d0, U);
    // // std::cout << "Rollout complete." << std::endl;

    // // generate a reference trajectory
    // Vector_12d_Traj X_ref = controller.generate_reference_trajectory(x0_sys.head<4>());

    // // test cost function
    // double J = controller.cost_function(X_ref, sol, U);
    // std::cout << "cost: " << J << std::endl;

    // test the monte carlo simulation
    // MC_Result mc = controller.monte_carlo(x0_sys, p0_feet, d0, 10);
    // std::cout << "Monte Carlo complete." << std::endl;

    ////////////////////////////////// Nominal testing //////////////////////////////////

    // compute time stuff
    double duration = config_file["SIM"]["duration"].as<double>();
    int N_sim = std::floor(duration / controller.params.dt);
    double dt = controller.params.dt;

    // run the simulation
    RHC_Result rhc_res;
    Solution sol;

    // run the simulation
    // Vector_8d xk_sys;
    // Vector_2d pk_feet;
    // Domain dk = d0;
    // for (int k = 0; k < N_sim; k++) {
    //     std::cout << "Simulation time: " << k * dt << " sec" << std::endl;
    // }
    
    // do a simulation
    rhc_res = controller.sampling_predictive_control(x0_sys, p0_feet, d0);
    sol = rhc_res.S;

    ////////////////////////////////// Logging //////////////////////////////////

    // unpack the solution
    Vector_1d_Traj t = sol.t;
    Vector_8d_Traj x_sys_t = sol.x_sys_t;
    Vector_8d_Traj x_leg_t = sol.x_leg_t;
    Vector_8d_Traj x_foot_t = sol.x_foot_t;
    Vector_4d_Traj u_t = sol.u_t;
    Vector_4d_Traj lambda_t = sol.lambda_t;
    Vector_2d_Traj tau_t = sol.tau_t;
    Domain_Traj domain_t = sol.domain_t;
    bool viability = sol.viability;

    // save the solution to a file
    std::string time_file = "../data/time.csv";
    std::string x_sys_file = "../data/state_sys.csv";
    std::string x_leg_file = "../data/state_leg.csv";
    std::string x_foot_file = "../data/state_foot.csv";
    std::string u_file = "../data/input.csv";
    std::string lambda_file = "../data/lambda.csv";
    std::string tau_file = "../data/tau.csv";
    std::string domain_file = "../data/domain.csv";

    int N_ = t.size();

    // save the solution to a file
    std::ofstream file;

    file.open(time_file);
    for (int i = 0; i < N_; i++) {
        file << t[i] << std::endl;
    }
    file.close();
    std::cout << "Saved time trajectory." << std::endl;

    file.open(x_sys_file);
    for (int i = 0; i < N_; i++) {
        file << x_sys_t[i].transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved system state trajectory." << std::endl;

    file.open(x_leg_file);
    for (int i = 0; i < N_; i++) {
        file << x_leg_t[i].transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved leg state trajectory." << std::endl;

    file.open(x_foot_file);
    for (int i = 0; i < N_; i++) {
        file << x_foot_t[i].transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved foot state trajectory." << std::endl;

    file.open(u_file);
    for (int i = 0; i < N_; i++) {
        file << u_t[i].transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved control input trajectory." << std::endl;

    file.open(lambda_file);
    for (int i = 0; i < N_; i++) {
        file << lambda_t[i].transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved leg force trajectory." << std::endl;

    file.open(tau_file);
    for (int i = 0; i < N_; i++) {
        file << tau_t[i].transpose() << std::endl;
    }
    file.close();

    Domain domain_t_(2);
    Vector_2i domain_;
    file.open(domain_file);
    for (int i = 0; i < N_; i++) {
        
        domain_t_ = domain_t[i];

        if (domain_t_[0] == Contact::STANCE) {
            domain_[0] = 1;
        }
        else {
            domain_[0] = 0;
        }

        if (domain_t_[1] == Contact::STANCE) {
            domain_[1] = 1;
        }
        else {
            domain_[1] = 0;
        }

        file << domain_.transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved domain trajectory." << std::endl;

    std::cout << "Final boss defeated." << std::endl;

    return 0;
}
