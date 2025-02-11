#pragma once

// standard libraries
#include <Eigen/Dense>
#include <vector>
#include <cmath>

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
using Vector_3d = Eigen::Matrix<double, 3, 1>;
using Vector_6d = Eigen::Matrix<double, 6, 1>;
using Vector_8d = Eigen::Matrix<double, 8, 1>;
using Vector_12d = Eigen::Matrix<double, 12, 1>;

using Matrix_3d = Eigen::Matrix<double, 3, 3>;
using Matrix_6d = Eigen::Matrix<double, 6, 6>;
using Matrix_8d = Eigen::Matrix<double, 8, 8>;
using Matrix_12d = Eigen::Matrix<double, 12, 12>;

// Dynamic arrays
using Vector_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix_d = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

// Trajectory types
using Vector_1i_Traj = std::vector<int>;
using Vector_1d_Traj = std::vector<double>;
using Vector_2i_Traj = std::vector<Vector_2i>;
using Vector_3d_Traj = std::vector<Vector_3d>;
using Vector_6d_Traj = std::vector<Vector_6d>;
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
    double pz_offset;        // z offset of the feet [m]
    double r_min;             // minimum rest length [m]
    double r_max;             // maximum rest length [m]
    double theta_x_min;       // minimum leg angle from vertical [rad]
    double theta_x_max;       // maximum leg angle from vertical [rad]
    double theta_y_min;       // minimum leg angle from vertical [rad]
    double theta_y_max;       // maximum leg angle from vertical [rad]
    double rdot_lim;          // maximum leg extension velocity [m/s]
    double thetadot_lim;      // maximum leg angle velocity [rad/s]
    bool torque_ankle;        // enable ankle torque 
    double torque_ankle_lim;  // enable ankle torque 
    double torque_ankle_kp;   // proportional gain for ankle torque
    double torque_ankle_kd;   // derivative gain for ankle torque
    double friction_coeff;    // friction coefficient
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
    Matrix_6d Q_com;          // diagonal elements of com Q matrix
    Matrix_6d Qf_com;         // diagonal elements of com Qf matrix
    Matrix_12d Q_leg;         // diagonal elements of leg Q matrix
    Matrix_12d Qf_leg;        // diagonal elements of leg Qf matrix
    Matrix_6d R;              // diagonal elements of R matrix
    Matrix_6d R_rate;         // diagonal elements of R rate matrix
    bool limits_enabled;      // enable kineamtic limits cost
    bool pos_limits;          // enable position limits cost
    bool vel_limits;          // enable velocity limits cost
    double w_pos_limits;      // kinematic limits cost
    double w_vel_limits;      // velocity limits cost
    bool friction_enabled;    // enable friction cone cost
    double w_friction;        // friction cone cost
    bool gait_enabled;        // enable gait cycle cost
    double w_gait;            // gait cycle cost
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
    bool saturate;      // saturate the control inputs
};

// dynamics solution struct (flow result)
struct Solution
{
    Vector_1d_Traj t;         // time trajectory
    Vector_12d_Traj x_sys_t;  // system state trajectory
    Vector_12d_Traj x_leg_t;  // leg state trajectory
    Vector_12d_Traj x_foot_t; // foot state trajectory
    Vector_6d_Traj u_t;       // interpolated control input trajectory
    Vector_6d_Traj lambda_t;  // leg force trajectory
    Vector_6d_Traj tau_t;     // ankle torque trajectory
    Domain_Traj domain_t;     // domain trajectory
    bool viability;           // viability of the trajectory
};

// Reference type
struct ReferenceGlobal  // throughout the entire simulation
{
    Vector_1d_Traj t_ref;     // time reference
    Vector_3d_Traj p_com_ref; // com position reference
    Vector_12d X_leg_ref;     // leg state reference
    int N_ref;                // number of reference points
    bool dynamic_cycle;       // dynamic gait cycle
    double T_cycle_static;    // gait cycle total time (in P2 sense)
    double contact_static;    // contact ratio (% of the gait cycle where the leg is in contact)
    double phase_static;      // phase offset (static gait)
    double v_min;             // minimum com velocity for spline (dynamic gait)
    double v_max;             // maximum com velocity for spline (dynamic gait)
    Vector_1d_List T_coeffs;  // quadratic coefficients for gait period
    Vector_1d_List c_coeffs;  // cubic coefficients for contact ratio
};

struct ReferenceLocal  // throughtout a horizon
{
    Vector_6d_Traj  X_com_ref;    // com position reference
    Vector_12d_Traj X_leg_ref;    // leg state reference
    Vector_2i_Traj  domain_ref;   // domain reference
    double T_cycle;               // gait cycle total time
    double contact_ratio;         // contact ratio
    double phase;                 // phase offset
};

// ***********************************************************************************
// Bundles
// ***********************************************************************************

// Bundle of Trajectories
using Vector_3d_Traj_Bundle = std::vector<Vector_3d_Traj>;
using Vector_6d_Traj_Bundle = std::vector<Vector_6d_Traj>;
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
    Vector_12d xdot;    // state derivative
    Vector_6d lambdas; // leg force
    Vector_6d taus;    // ankle torque
};

// monte carlo simulation result
struct MC_Result
{
    Solution_Bundle S;       // Solutions
    Vector_6d_Traj_Bundle U; // Control Inputs  
    Vector_1d_List J;        // Costs
};

// result of the sampling predictive control
struct RHC_Result
{
    Solution S;       // Solution
    Vector_6d_Traj U; // Optimal Control Inputs  
};
