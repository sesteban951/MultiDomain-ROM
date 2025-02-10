// standard includes
#include <iostream>
#include <chrono>
#include <fstream>

// package includes
#include "Eigen/Dense"
#include "yaml-cpp/yaml.h"

// custom includes
#include "../inc/types/types.h"
#include "../inc/dynamics/dynamics.h"
#include "../inc/control/control.h"

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

    // // Vector_4d u_const;
    // Vector_4d u_const;
    // Vector_2d u_L, u_R;
    // u_L << 0.0, 0.0;
    // u_R << 0.0, -0.0;
    // u_const << u_L, u_R;

    // // example rollout of the dynamics
    // int N = 140;
    // Vector_1d_Traj T_x(N);
    // for (int i = 0; i < N; i++) {
    //     T_x[i] = i * 0.01;
    // }

    // int Nu = 70;
    // Vector_1d_Traj T_u(Nu);
    // Vector_4d_Traj U(Nu);
    // for (int i = 0; i < Nu; i++) {
    //     T_u[i] = i * 0.02;
    //     U[i] = u_const;
    // }

    // // query the dynamics
    // Dynamics_Result res = dynamics.dynamics(x0_sys, u_const, p0_feet, d0);
    // Vector_8d xdot = res.xdot;
    // std::cout << "xdot: " << xdot.transpose() << std::endl;

    // // leg state
    // Vector_8d x_leg = dynamics.compute_leg_state(x0_sys, u_const, p0_feet, d0);
    // std::cout << "x_leg: " << x_leg.transpose() << std::endl;

    // // foot state
    // Vector_8d x_foot = dynamics.compute_foot_state(x0_sys, x_leg, p0_feet, d0);
    // std::cout << "x_foot: " << x_foot.transpose() << std::endl;
    
    // // allocate solution conatiners
    // Solution sol;
    // dynamics.resizeSolution(sol, T_x);

    // // do the rollouts
    // dynamics.RK3_rollout(T_x, T_u, x0_sys, p0_feet, d0, U, sol);
    // std::cout << "Rollout complete." << std::endl;
    
    // // generate inputs
    // Vector_4d_Traj_Bundle U_bundle = controller.sample_input_trajectory(1000);

    // // update the distribution parameters
    // controller.update_distribution_params(U_bundle);
    // controller.update_distribution_params(U_bundle);
    // std::cout << "updated distribution." << std::endl;

    // // Do a rollout of the dynamics
    // U = U_bundle[0];
    // dynamics.RK3_rollout(T_x, T_u, x0_sys, p0_feet, d0, U, sol);
    // std::cout << "Rollout complete." << std::endl;

    // generate a reference trajectory
    // Reference ref = controller.generate_reference_trajectory(x0_sys.head<4>());

    // // test cost function
    // double J = controller.cost_function(X_ref, sol, U);
    // std::cout << "cost: " << J << std::endl;

    // test the monte carlo simulation
    // Solution sol;
    // dynamics.resizeSolution(sol, controller.params.T_x);
    // MC_Result mc = controller.monte_carlo(0.0, x0_sys, p0_feet, d0, controller.params.K);
    // std::cout << "Monte Carlo complete." << std::endl;

    // do a simulation
    // RHC_Result rhc_res = controller.sampling_predictive_control(0.0, x0_sys, p0_feet, d0);
    // Solution sol = rhc_res.S;

    ////////////////////////////////// Nominal testing //////////////////////////////////

    // compute time stuff
    double duration = config_file["SIM"]["duration"].as<double>();
    int N_sim = std::floor(duration / controller.params.dt_x);
    double dt = controller.params.dt_x;
    int Nu = controller.params.N_u;

    // create the one step integration time vector
    Vector_1d_Traj T_x_onestep(2);
    T_x_onestep[0] = 0.0;
    T_x_onestep[1] = dt;
    
    // solution container to use for logging
    Solution sol;
    sol.t.resize(N_sim);
    sol.x_sys_t.resize(N_sim);
    sol.x_leg_t.resize(N_sim);
    sol.x_foot_t.resize(N_sim);
    sol.u_t.resize(N_sim);
    sol.lambda_t.resize(N_sim);
    sol.tau_t.resize(N_sim);
    sol.domain_t.resize(N_sim);
    
    Solution sol_rhc;
    controller.dynamics.resizeSolution(sol_rhc, controller.params.T_x);

    // to time time stuff
    double T_rhc = 0.0;
    double T_total = 0.0;

    // run the simulation
    RHC_Result rhc_res;
    Vector_4d_Traj U_opt(Nu), U_opt_(Nu);
    Vector_d U_opt_vec(2 * Nu * controller.dynamics.n_leg);
    Vector_8d xk_sys = x0_sys;
    Vector_4d pk_feet = p0_feet;
    Vector_8d xk_feet;
    Domain dk = d0;
    double t_sim;
    auto t0_tot = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < N_sim; k++) {
        t_sim = k * dt;
        std::cout << "Sim time: " << t_sim << " sec" << std::endl;

        // do predictive control 
        auto t0 = std::chrono::high_resolution_clock::now();
        rhc_res = controller.sampling_predictive_control(t_sim, xk_sys, pk_feet, dk);
        auto tf = std::chrono::high_resolution_clock::now();
        T_rhc += std::chrono::duration<double>(tf - t0).count();

        // extract the optimal input sequence
        U_opt = rhc_res.U;

        for (int i = 0; i < controller.params.N_u; i++) {
            U_opt_[i] = U_opt[i];
            U_opt_vec.segment<4>(4*i) = U_opt_[i];
        }

        // integrate the dynamics
        sol_rhc = dynamics.RK3_rollout(T_x_onestep, 
                             controller.params.T_u, 
                             xk_sys, pk_feet, dk, U_opt_);
        
        // save into solution bundle
        sol.t[k] = t_sim;
        sol.x_sys_t[k] = sol_rhc.x_sys_t[1];
        sol.x_leg_t[k] = sol_rhc.x_leg_t[1];
        sol.x_foot_t[k] = sol_rhc.x_foot_t[1];
        sol.u_t[k] = sol_rhc.u_t[1];
        sol.lambda_t[k] = sol_rhc.lambda_t[1];
        sol.tau_t[k] = sol_rhc.tau_t[1];
        sol.domain_t[k] = sol_rhc.domain_t[1];

        // update the state
        xk_sys = sol_rhc.x_sys_t[1];
        dk = sol_rhc.domain_t[1];        
        xk_feet = sol_rhc.x_foot_t[1];
        pk_feet(0) = xk_feet(0);
        pk_feet(1) = xk_feet(1);
        pk_feet(2) = xk_feet(4);
        pk_feet(3) = xk_feet(5);
    }
    auto tf_tot = std::chrono::high_resolution_clock::now();

    // print some info
    double T_tot = std::chrono::duration<double>(tf_tot - t0_tot).count();
    std::cout << "Total time: " << T_tot << " [sec]" << std::endl;

    double T_rhc_avg = (T_rhc / N_sim) * 1000.0;
    std::cout << "Average time for RHC: " << T_rhc_avg << " [ms], " << 1000.0 / T_rhc_avg << " [Hz]" << std::endl;

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
    std::string time_file = "../data/2D/time.csv";
    std::string x_sys_file = "../data/2D/state_sys.csv";
    std::string x_leg_file = "../data/2D/state_leg.csv";
    std::string x_foot_file = "../data/2D/state_foot.csv";
    std::string u_file = "../data/2D/input.csv";
    std::string lambda_file = "../data/2D/lambda.csv";
    std::string tau_file = "../data/2D/tau.csv";
    std::string domain_file = "../data/2D/domain.csv";

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
