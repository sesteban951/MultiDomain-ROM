// standard includes
#include <iostream>
#include <chrono>
#include <fstream>

// package includes
#include "Eigen/Dense"
#include "yaml-cpp/yaml.h"

// custom includes
#include "../inc/types_3D.h"
#include "../inc/dynamics_3D.h"
#include "../inc/control_3D.h"

int main()
{
    // load parameters from yaml file
    std::string config_file_path = "../config/config_3D.yaml";
    YAML::Node config_file = YAML::LoadFile(config_file_path);
    
    // create dynamics object
    Dynamics dynamics(config_file);
    Controller controller(config_file);

    // initial conditions
    Vector_12d x0_sys;
    Vector_3d p0_left, p0_right;
    Vector_6d p0_feet;

    std::vector<double> x0_temp = config_file["STATE"]["x0"].as<std::vector<double>>();
    x0_sys << x0_temp[0],    // px_com
              x0_temp[1],    // py_com
              x0_temp[2],    // pz_com
              x0_temp[3],    // vx_com
              x0_temp[4],    // vy_com
              x0_temp[5],    // vz_com
              x0_temp[6],    // l0_command (left)
              x0_temp[7],    // thetadot_x_command (left)
              x0_temp[8],    // thetadot_y_command (left)
              x0_temp[9],    // l0_command (right)
              x0_temp[10],   // thetadot_x_command (right)
              x0_temp[11];   // thetadot_y_command (right)
    p0_left <<  0.0,  0.1, 0.0;
    p0_right << 0.0, -0.1, 0.0;
    p0_feet << p0_left, p0_right;
    Domain d0(2, Contact::SWING);

    ////////////////////////////////// Simulation testing //////////////////////////////////

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

    // for timing stuff
    double T_rhc = 0.0;
    double T_tot = 0.0;

    // run the simulation
    RHC_Result rhc_res;
    Vector_6d_Traj U_opt(Nu), U_opt_(Nu);
    Vector_d U_opt_vec(3 * Nu * controller.dynamics.n_leg);
    Vector_12d xk_sys = x0_sys;
    Vector_6d pk_feet = p0_feet;
    Vector_12d xk_feet;
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
            U_opt_vec.segment<6>(6*i) = U_opt_[i];
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
        pk_feet(2) = xk_feet(2);
        pk_feet(3) = xk_feet(6);
        pk_feet(4) = xk_feet(7);
        pk_feet(5) = xk_feet(8);
    }
    auto tf_tot = std::chrono::high_resolution_clock::now();

    // print some info
    T_tot = std::chrono::duration<double>(tf_tot - t0_tot).count();
    std::cout << "\nTotal time: " << T_tot << " sec" << std::endl;

    double rt_ratio = (dt * N_sim) / T_tot;
    std::cout << "Real-time ratio: " << rt_ratio *100 << "%" << std::endl;

    double T_rhc_avg = (T_rhc / N_sim) * 1000.0;
    std::cout << "Average time for RHC: " << T_rhc_avg << " [ms], " << 1000.0 / T_rhc_avg << " [Hz]" << std::endl;

    ////////////////////////////////// Logging //////////////////////////////////

    // unpack the solution
    Vector_1d_Traj t = sol.t;
    Vector_12d_Traj x_sys_t = sol.x_sys_t;
    Vector_12d_Traj x_leg_t = sol.x_leg_t;
    Vector_12d_Traj x_foot_t = sol.x_foot_t;
    Vector_6d_Traj u_t = sol.u_t;
    Vector_6d_Traj lambda_t = sol.lambda_t;
    Vector_6d_Traj tau_t = sol.tau_t;
    Domain_Traj domain_t = sol.domain_t;
    bool viability = sol.viability;

    // save the solution to a file
    std::string time_file = "../data/3D/time.csv";
    std::string x_sys_file = "../data/3D/state_sys.csv";
    std::string x_leg_file = "../data/3D/state_leg.csv";
    std::string x_foot_file = "../data/3D/state_foot.csv";
    std::string u_file = "../data/3D/input.csv";
    std::string lambda_file = "../data/3D/lambda.csv";
    std::string tau_file = "../data/3D/tau.csv";
    std::string domain_file = "../data/3D/domain.csv";

    int N_ = t.size();

    // save the solution to a file
    std::ofstream file;

    file.open(time_file);
    for (int i = 0; i < N_; i++) {
        file << t[i] << std::endl;
    }
    file.close();
    std::cout << "\nSaved time trajectory." << std::endl;

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

        domain_[0] = (domain_t_[0] == Contact::STANCE) ? 1 : 0;
        domain_[1] = (domain_t_[1] == Contact::STANCE) ? 1 : 0;

        file << domain_.transpose() << std::endl;
    }
    file.close();
    std::cout << "Saved domain trajectory." << std::endl;

    std::cout << "Final boss defeated." << std::endl;

    return 0;
}
