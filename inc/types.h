#pragma once

// standard libraries
#include <Eigen/Dense>
#include <vector>

// ***********************************************************************************
// Domain Class
// ***********************************************************************************

enum class Contact {SWING, STANCE};   // STANCE: in contanct, SWING: not in contact
using Domain = std::vector<Contact>;  // For multiple legs

// ***********************************************************************************
// Leg Identification Class
// ***********************************************************************************

using Leg_Idx = int;

// ***********************************************************************************
// ARRAYS
// ***********************************************************************************

// Fixed size arrays
using Vector_2i = Eigen::Matrix<int, 2, 1>;
using Vector_2d = Eigen::Matrix<double, 2, 1>;
using Vector_4d = Eigen::Matrix<double, 4, 1>;
using Vector_6d = Eigen::Matrix<double, 6, 1>;
using Vector_8d = Eigen::Matrix<double, 8, 1>;
using Vector_12d = Eigen::Matrix<double, 12, 1>;
// using Vector_2i = Eigen::Vector<int, 2>;
// using Vector_2d = Eigen::Vector<double, 2>;
// using Vector_4d = Eigen::Vector<double, 4>;
// using Vector_6d = Eigen::Vector<double, 6>;
// using Vector_8d = Eigen::Vector<double, 8>;
// using Vector_12d = Eigen::Vector<double, 12>;

using Matrix_2d = Eigen::Matrix<double, 2, 2>;
using Matrix_4d = Eigen::Matrix<double, 4, 4>;
using Matrix_8d = Eigen::Matrix<double, 8, 8>;
using Matrix_12d = Eigen::Matrix<double, 12, 12>;

// Dynamic arrays
// using Vector_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix_d = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

using Vector_d = Eigen::Vector<double, Eigen::Dynamic>;

// Trajectory types
using Vector_1i_Traj = std::vector<int>;
using Vector_1d_Traj = std::vector<double>;
using Vector_2i_Traj = std::vector<Vector_2i>;
using Vector_2d_Traj = std::vector<Vector_2d>;
using Vector_4d_Traj = std::vector<Vector_4d>;
using Vector_8d_Traj = std::vector<Vector_8d>;
using Vector_12d_Traj = std::vector<Vector_12d>;

using Domain_Traj = std::vector<Domain>;

// List types
using Vector_1i_List = std::vector<int>;
using Vector_1d_List = std::vector<double>;

// ***********************************************************************************
// STRUCTS
// ***********************************************************************************

// struct to store system parameters
struct SystemParams
{
    double m;                 // mass [kg]
    double g;                 // gravity [m/s^2]
    double k;                 // spring constant [N/m]
    double b;                 // damping constant [Ns/m]
    double l0;                // nominal rest length [m]
    double r_min;             // minimum rest length [m]
    double r_max;             // maximum rest length [m]
    double theta_min;         // minimum leg angle from vertical [rad]
    double theta_max;         // maximum leg angle from vertical [rad]
    double rdot_lim;          // maximum leg extension velocity [m/s]
    double thetadot_lim;      // maximum leg angle velocity [rad/s]
    bool torque_ankle;        // enable ankle torque 
    double torque_ankle_lim;  // enable ankle torque 
    double torque_ankle_kp;   // proportional gain for ankle torque
    double torque_ankle_kd;   // derivative gain for ankle torque
    char interp;              // control interpolation type
};

// struct to hold control parameters
struct ControlParams
{
    int N_x;                  // number of system dynamics integration steps
    int N_u;                  // number of control points
    double dt_x;              // time step [sec]
    double dt_u;              // time step [sec], for the control inputs
    Vector_1d_Traj T_x;       // time array for system dynamics
    Vector_1d_Traj T_u;       // time array for control inputs
    int K;                    // number of parallel 
    int N_elite;              // number of elite control sequences
    int CEM_iters;            // number of CEM iterations
    Matrix_4d Q_com;          // diagonal elements of com Q matrix
    Matrix_4d Qf_com;         // diagonal elements of com Qf matrix
    Matrix_8d Q_leg;          // diagonal elements of leg Q matrix
    Matrix_8d Qf_leg;         // diagonal elements of leg Qf matrix
    Matrix_4d R;              // diagonal elements of R matrix
    Matrix_4d R_rate;         // diagonal elements of R rate matrix
    double gait_cycle_weight; // weight for gait cycle cost
    bool log_barrier_enabled; // enable log barrier function
    double r_weight;          // weight for log barrier function
    double theta_weight;      // weight for log barrier function
    double rdot_weight;       // weight for log barrier function
    double thetadot_weight;   // weight for log barrier function
};

// Distribution struct
struct GaussianDistribution
{
    Vector_d mean;      // mean of the distribution
    Matrix_d cov;       // covariance of the distribution
    double epsilon;     // small value to add to the diagonal of the covariance
    bool diag_cov;      // choose to enfoce diagonal covariance
    bool seed_enabled;  // enable random number generator seed
    int seed;           // random number generator seed
};

// dynamics solution struct (flow result)
struct Solution
{
    Vector_1d_Traj t;        // time trajectory
    Vector_8d_Traj x_sys_t;  // system state trajectory
    Vector_8d_Traj x_leg_t;  // leg state trajectory
    Vector_8d_Traj x_foot_t; // foot state trajectory
    Vector_4d_Traj u_t;      // interpolated control input trajectory
    Vector_4d_Traj lambda_t; // leg force trajectory
    Vector_2d_Traj tau_t;         // ankle torque trajectory
    Domain_Traj domain_t;    // domain trajectory
    bool viability;          // viability of the trajectory
};

// Reference type
struct Reference
{
    Vector_4d_Traj X_com_ref; // reference trajectory
    Vector_8d_Traj X_leg_ref; // reference trajectory
    Vector_2i_Traj D_ref;     // reference trajectory
};

// ***********************************************************************************
// Bundles
// ***********************************************************************************

// Bundle of Trajectories
using Vector_2d_Traj_Bundle = std::vector<Vector_2d_Traj>;
using Vector_4d_Traj_Bundle = std::vector<Vector_4d_Traj>;
using Vector_8d_Traj_Bundle = std::vector<Vector_8d_Traj>;
using Vector_12d_Traj_Bundle = std::vector<Vector_12d_Traj>;

// Bundle of Solutions
using Solution_Bundle = std::vector<Solution>;

// ***********************************************************************************
// Result Types
// ***********************************************************************************

// dynamics integration step result
struct Dynamics_Result
{
    Vector_8d xdot;    // state derivative
    Vector_4d lambdas; // leg force
    Vector_2d taus;    // ankle torque
};

struct MC_Result
{
    Solution_Bundle S;       // Solutions
    Vector_4d_Traj_Bundle U; // Control Inputs  
    Vector_1d_List J;        // Costs
};

struct RHC_Result
{
    Solution S;       // Solution
    Vector_4d_Traj U; // Optimal Control Inputs  
};
