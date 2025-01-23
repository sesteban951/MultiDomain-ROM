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
    Vector_8d x0;
    Vector_2d p0_left, p0_right;
    Vector_2d_List p0_feet(2);

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
    Domain d0(2, Contact::SWING);

    // Vector_2d_List u(2);
    // u[0] << 0.04, 0.5;
    // u[1] << -0.04, -0.5;

    // // example rollout of the dynamics
    // int N = 300;
    // Vector_1d_Traj T_x(N);
    // for (int i = 0; i < N; i++) {
    //     T_x[i] = i * 0.01;
    // }

    // int Nu = 150;
    // Vector_1d_Traj T_u(Nu);
    // Vector_2d_Traj U(Nu);
    // Vector_2d U1, U2;
    // Vector_2d_List Ui(2);
    // U1 << 0., -0.0;
    // U2 << 0., 0.0;
    // for (int i = 0; i < Nu; i++) {
    //     T_u[i] = i * 0.02;
    //     Ui[0] = U1;
    //     Ui[1] = U2;
    //     U[i] = Ui;
    // }

    // // Do a rollout of the dynamics
    // Solution sol = dynamics.RK3_rollout(T_x, T_u, x0, p0_feet, d0, U);

    // // generate a reference trajectory
    // Vector_12d_List X_ref = controller.generate_reference_trajectory(x0.head<4>());

    // double J = controller.cost_function(X_ref, sol, U);

    // std::cout << "cost: " << J << std::endl;

    Solution sol = controller.sampling_predictive_control(x0, p0_feet, d0);

    // unpack the solution
    Vector_1d_List t = sol.t;
    Vector_8d_List x_sys_t = sol.x_sys_t;
    Vector_4d_Traj x_leg_t = sol.x_leg_t;
    Vector_4d_Traj x_foot_t = sol.x_foot_t;
    Vector_2d_Traj u_t = sol.u_t;
    Domain_List domain_t = sol.domain_t;
    bool viability = sol.viability;

    // save the solution to a file
    std::string time_file = "../data/time.csv";
    std::string x_sys_file = "../data/state_sys.csv";
    std::string x_leg_file = "../data/state_leg.csv";
    std::string x_foot_file = "../data/state_foot.csv";
    std::string u_file = "../data/input.csv";
    std::string domain_file = "../data/domain.csv";

    // save the solution to a file
    std::ofstream file;

    file.open(time_file);
    for (int i = 0; i < controller.params.N; i++) {
        file << t[i] << std::endl;
    }
    file.close();

    std::cout << "\nSaved time trajectory." << std::endl;

    file.open(x_sys_file);
    for (int i = 0; i < controller.params.N; i++) {
        file << x_sys_t[i].transpose() << std::endl;
    }
    file.close();

    std::cout << "Saved system state trajectory." << std::endl;

    Vector_8d x_leg_;
    Vector_4d_List x_leg_k;
    file.open(x_leg_file);
    for (int i = 0; i < controller.params.N; i++) {
        x_leg_k = x_leg_t[i];
        x_leg_ << x_leg_k[0], x_leg_k[1];
        file << x_leg_.transpose() << std::endl;
    }
    file.close();

    std::cout << "Saved leg state trajectory." << std::endl;

    Vector_8d x_foot_;
    Vector_4d_List x_foot_k;
    file.open(x_foot_file);
    for (int i = 0; i < controller.params.N; i++) {
        x_foot_k = x_foot_t[i];
        x_foot_ << x_foot_k[0], x_foot_k[1];
        file << x_foot_.transpose() << std::endl;
    }
    file.close();

    std::cout << "Saved foot state trajectory." << std::endl;

    Vector_4d u_;
    Vector_2d_List u_k;
    file.open(u_file);
    for (int i = 0; i < controller.params.N; i++) {
        u_k = u_t[i];
        u_ << u_k[0], u_k[1];
        file << u_.transpose() << std::endl;
    }
    file.close();

    std::cout << "Saved input trajectory." << std::endl;

    Domain domain_t_(2);
    Vector_2i domain_;
    file.open(domain_file);
    for (int i = 0; i < controller.params.N; i++) {
        
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

    return 0;
}
