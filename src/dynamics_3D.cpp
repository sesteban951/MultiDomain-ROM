#include "../inc/dynamics_3D.h"


// populate the system parameters
Dynamics::Dynamics(YAML::Node config_file)
{
    // set the system parameters
    this->params.m = config_file["SYS_PARAMS"]["m"].as<double>();
    this->params.g = config_file["SYS_PARAMS"]["g"].as<double>();
    this->params.k = config_file["SYS_PARAMS"]["k"].as<double>();
    this->params.b = config_file["SYS_PARAMS"]["b"].as<double>();
    this->params.l0 = config_file["SYS_PARAMS"]["l0"].as<double>();
    this->params.r_min = config_file["SYS_PARAMS"]["r_min"].as<double>();
    this->params.r_max = config_file["SYS_PARAMS"]["r_max"].as<double>();
    this->params.theta_min = config_file["SYS_PARAMS"]["theta_min"].as<double>();
    this->params.theta_max = config_file["SYS_PARAMS"]["theta_max"].as<double>();
    this->params.rdot_lim = config_file["SYS_PARAMS"]["rdot_lim"].as<double>();
    this->params.thetadot_lim = config_file["SYS_PARAMS"]["thetadot_lim"].as<double>();
    this->params.torque_ankle = config_file["SYS_PARAMS"]["torque_ankle"].as<bool>();
    this->params.torque_ankle_lim = config_file["SYS_PARAMS"]["torque_ankle_lim"].as<double>();
    this->params.torque_ankle_kp = config_file["SYS_PARAMS"]["torque_ankle_kp"].as<double>();
    this->params.torque_ankle_kd = config_file["SYS_PARAMS"]["torque_ankle_kd"].as<double>();
    this->params.interp = config_file["SYS_PARAMS"]["interp"].as<char>();
}


// NonLinear Dynamics function, xdot = f(x, u, d)
Dynamics_Result_3D Dynamics::dynamics(const Vector_12d& x_sys, const Vector_6d& u, const Vector_6d& p_feet, const Domain& d) 
{
    // want to compute the state derivative
    Vector_12d xdot;     // state derivative, [xdot_com, xdot_command]
    Vector_6d lambdas;   // leg force vectors
    Vector_6d taus;      // ankle torques

    // access some system parameters
    double m = this->params.m;
    double g = this->params.g;
    double k = this->params.k;
    double b = this->params.b;

    // unpack the state vector
    Vector_6d x_com;
    x_com = x_sys.segment<6>(0);

    // vectors to use in the calculations
    Vector_3d v_com, a_com, g_vec, f_com;
    Vector_6d v_commands;
    v_com = x_com.segment<3>(3);
    f_com.setZero();
    v_commands = u;

    // temporary wrench variables
    Vector_3d leg_force;      //  one leg force vector applied to the COM
    Vector_3d ankle_torque;   //  one ankle torque vector applied to the COM
    double tau_x_ankle, tau_y_ankle;
    leg_force.setZero();
    ankle_torque.setZero();
    tau_x_ankle = 0.0;
    tau_y_ankle = 0.0;

    // contact state
    Contact c;

    // loop through the legs to compute the forces and torques on COM
    for (int i = 0; i < this->n_leg; i++)
    {
        // contact state
        c = d[i];

        // The i'th leg is in NOT in contact, no forces or torques
        if (c == Contact::SWING) {
            // in swing, no forces or torques on the COM
            leg_force.setZero();
            ankle_torque.setZero();

            // save the leg force and ankle torque
            lambdas.segment<3>(3*i) = leg_force;
            taus.segment<3>(3*i) = ankle_torque;
        }

        // The i'th leg is in contact, there are forces and torques
        else if (c == Contact::STANCE) {

            // in stance, compute the forces and torques on the COM
            Vector_3d f_com_ankle;

            // COM state
            Vector_3d p_com;
            p_com = x_com.segment<3>(0);

            // LEG Command State
            Vector_3d x_leg_command;
            x_leg_command = x_sys.segment<3>(6 + 3*i);

            // LEG input
            Vector_3d u_leg;
            u_leg = u.segment<3>(3*i);

            // foot position
            Vector_3d p_foot;
            p_foot = p_feet.segment<3>(3*i);
            
            // useful variables
            Vector_3d r_vec, rdot_vec, r_hat;
            double r_norm, rdot_norm;
            double l0_command, l0dot_command;

            // leg r vector
            r_vec = p_foot - p_com;
            r_norm = r_vec.norm();
            r_hat = r_vec / r_norm;
            
            // leg rdot vector
            rdot_vec = -v_com;     // assumes no sliding feet
            rdot_norm = -v_com.dot(r_hat);

            // compute the force along the leg
            l0_command = x_leg_command(0);
            l0dot_command = u_leg(0);
            leg_force = -r_hat * (k * (l0_command - r_norm) + b * (l0dot_command - rdot_norm)); 

            // TODO: think about how having torque on both legs when both are in stance -- overactuated
            // compute the ankle torque
            if (this->params.torque_ankle == true) {

                // useful leg pos/vel variables
                double r_x, r_y, r_z, rdot_x, rdot_y, rdot_z;
                r_x = r_vec(0);
                r_y = r_vec(1);
                r_z = r_vec(2);
                rdot_x = rdot_vec(0);
                rdot_y = rdot_vec(1);
                rdot_z = rdot_vec(2);

                // actual angle state and commanded angle states
                double theta_x, theta_x_command, thetadot_x, thetadot_x_command;
                double theta_y, theta_y_command, thetadot_y, thetadot_y_command;
                theta_x =  std::atan2(r_y, -r_z);
                theta_y = -std::atan2(r_x, -r_z);
                thetadot_x = (r_y * rdot_z - r_z * rdot_y) / (r_z*r_z + r_y*r_y);
                thetadot_y = (r_z * rdot_x - r_x * rdot_z) / (r_z*r_z + r_x*r_x);
                theta_x_command = x_leg_command(1);
                theta_y_command = x_leg_command(2);
                thetadot_x_command = u_leg(1);
                thetadot_y_command = u_leg(2);

                // compute the ankle torque
                double kp = this->params.torque_ankle_kp;
                double kd = this->params.torque_ankle_kd;
                tau_x_ankle = kp * (theta_x_command - theta_x) + kd * (thetadot_x_command - thetadot_x);
                tau_y_ankle = kp * (theta_y_command - theta_y) + kd * (thetadot_y_command - thetadot_y);

                // saturate the torque
                tau_x_ankle = std::max(-this->params.torque_ankle_lim, std::min(this->params.torque_ankle_lim, tau_x_ankle));
                tau_y_ankle = std::max(-this->params.torque_ankle_lim, std::min(this->params.torque_ankle_lim, tau_y_ankle));

                // compute the required force at the COM to generate the torque (note that a tau_z torque is induced)
                Matrix_3d A_inv;
                Vector_3d b;
                double sigma1 = r_norm * r_norm * r_z;
                double sigma2 = r_norm * r_norm;
                A_inv << -(r_x * r_y)/sigma1,         -(r_y*r_y + r_z*r_z)/sigma1, -r_x/sigma2,
                          (r_x*r_x + r_z*r_z)/sigma1,  (r_x * r_y)/sigma1,         -r_y/sigma2,
                          -r_y/sigma2,                  r_x/sigma2,                -r_z/sigma2;
                b << tau_x_ankle, tau_y_ankle, 0.0;
                f_com_ankle = A_inv * b;

                // NOTE: this method induces a torque about the z-axis, which is not actuated, we aim to minimize tau_z induced
                double tau_z_ankle = r_y * f_com_ankle(0) - r_x * f_com_ankle(1);
                ankle_torque << tau_x_ankle, tau_y_ankle, tau_z_ankle;
            }
            // no ankle torque
            else {
                ankle_torque.setZero();
                f_com_ankle.setZero();
            }

            // save the leg force and ankle torque
            lambdas.segment<3>(3*i) = leg_force;
            taus.segment<3>(3*i) = ankle_torque;

            // add the leg force and torque to the COM dynamics
            f_com += leg_force + f_com_ankle;
        }
    }

    // build the gravity vector
    g_vec << 0.0, 0.0, -g;

    // compute acceleration of the center of mass
    a_com << (1/m) * f_com + g_vec;

    // final form of the state derivative
    xdot << v_com, a_com, v_commands;

    // pack into the struct
    Dynamics_Result_3D res;
    res.xdot = xdot;
    res.lambdas = lambdas;
    res.taus = taus;

    // return the state derivative
    return res;
}


// compute the leg state
Vector_12d Dynamics::compute_leg_state(const Vector_12d& x_sys, const Vector_6d& u, const Vector_6d& p_feet, const Domain& d)
{
    // want to compute the leg states
    Vector_12d x_legs;
    Vector_6d x_leg;

    // indivudual states
    double r, theta_x, theta_y, rdot, thetadot_x, thetadot_y;

    // contact state
    Contact c;

    // loop though legs and compute the leg states
    for (int i = 0; i < this->n_leg; i++)
    {
        // contact state
        c = d[i];

        // in swing (single integrator dynamics)
        if (c == Contact::SWING) {
        
            // leg command vector
            Vector_3d x_leg_command = x_sys.segment<3>(6 + 3*i);

            // leg input vector
            Vector_3d u_leg = u.segment<3>(3*i);

            // single integrator state
            r = x_leg_command(0);
            theta_x = x_leg_command(1);
            theta_y = x_leg_command(2);
            rdot = u_leg(0);
            thetadot_x = u_leg(1);
            thetadot_y = u_leg(2);

            // populate the polar state vector
            x_leg << r, theta_x, theta_y, rdot, thetadot_x, thetadot_y;
        }

        // in stance (leg state goverened by the COM dynamics)
        else if (c == Contact::STANCE) {

            // COM state
            Vector_3d p_com, v_com;
            p_com = x_sys.segment<3>(0);
            v_com = x_sys.segment<3>(3);

            // get the foot position
            Vector_3d p_foot = p_feet.segment<3>(3*i);

            // leg vectors
            Vector_3d r_vec, r_hat, rdot_vec;
            double r_norm, r_x, r_y, r_z, rdot_x, rdot_y, rdot_z;
            r_vec = p_foot - p_com;
            r_norm = r_vec.norm();
            r_hat = r_vec / r_norm;
            r_x = r_vec(0);
            r_y = r_vec(1);
            r_z = r_vec(2);

            rdot_vec = -v_com; // assumes no sliding feet
            rdot_x = rdot_vec(0);
            rdot_y = rdot_vec(1);
            rdot_z = rdot_vec(2);

            // compute the polar leg state 
            r = r_norm;
            rdot = -v_com.dot(r_hat);
            theta_x =  std::atan2(r_y, -r_z);
            theta_y = -std::atan2(r_x, -r_z);
            thetadot_x = (r_y * rdot_z - r_z * rdot_y) / (r_z*r_z + r_y*r_y);
            thetadot_y = (r_z * rdot_x - r_x * rdot_z) / (r_z*r_z + r_x*r_x);

            // populate the polar state vector
            x_leg << r, theta_x, theta_y, rdot, thetadot_x, thetadot_y;
        }

        // populate the leg state vector
        x_legs.segment<6>(6*i) = x_leg;
    }

    return x_legs;
}


// compute foot state in world frame
Vector_12d Dynamics::compute_foot_state(const Vector_12d& x_sys, const Vector_12d& x_legs, const Vector_6d& p_feet, const Domain& d)
{
    // want to compute the foot states
    Vector_12d x_feet;
    Vector_6d x_foot;

    // contact state
    Contact c;

    // loop though legs and compute the foot states
    for (int i = 0; i < this->n_leg; i++)
    {
        // contact state
        c = d[i];

        // Foot is in SWING
        if (c == Contact::SWING) {

            // com varaibles
            double px_com, py_com, pz_com, vx_com, vy_com, vz_com;
            px_com = x_sys(0);
            py_com = x_sys(1);
            pz_com = x_sys(2);
            vx_com = x_sys(3);
            vy_com = x_sys(4);
            vz_com = x_sys(5);

            // leg states
            Vector_6d x_leg = x_legs.segment<6>(6*i);
            double r, theta_x, theta_y, rdot, thetadot_x ,thetadot_y;
            r = x_leg(0);
            theta_x = x_leg(1);
            theta_y = x_leg(2);
            rdot = x_leg(3);
            thetadot_x = x_leg(4);
            thetadot_y = x_leg(5);

            // WARNING: you get singularities when |theta_x|, |theta_y| >= pi/2 
            // compute the foot state via fwd kinematics
            double n, d;
            d = std::tan(theta_x)*std::tan(theta_x) + std::tan(theta_y)*std::tan(theta_y) + 1;
            n = std::tan(theta_x)/(std::cos(theta_x)*std::cos(theta_x))*thetadot_x + std::tan(theta_y)/(std::cos(theta_y)*std::cos(theta_y))*thetadot_y;
            
            // foot positions and velocities (w.r.t. the COM frame)
            double px_foot, py_foot, pz_foot, vx_foot, vy_foot, vz_foot;
            pz_foot = -std::sqrt(r*r/d);             
            px_foot =  pz_foot * std::tan(theta_y);  
            py_foot = -pz_foot * std::tan(theta_x);  
 
            // foot velocities (w.r.t. the COM frame)
            vz_foot = -rdot/std::sqrt(d) + (r*n)/std::sqrt(d*d*d);
            vx_foot =  vz_foot * std::tan(theta_y) + pz_foot * thetadot_y / (std::cos(theta_y)*std::cos(theta_y));
            vy_foot = -vz_foot * std::tan(theta_x) - pz_foot * thetadot_x / (std::cos(theta_x)*std::cos(theta_x));

            // now in world frame
            px_foot += px_com;
            py_foot += py_com;
            pz_foot += pz_com;
            vx_foot += vx_com;
            vy_foot += vy_com;
            vz_foot += vz_com;

            // populate the foot state vector
            x_foot << px_foot, py_foot, pz_foot, vx_foot, vy_foot, vz_foot;
        }

        // Foot is in STANCE
        else if (c == Contact::STANCE) {
            // foot state is the same as the foot position w/ zero velocity
            Vector_3d p_foot = p_feet.segment<3>(3*i);
            x_foot << p_foot(0), p_foot(1), p_foot(2), 0.0, 0.0, 0.0;
        }

        // insert into the foot state vector
        x_feet.segment<6>(6*i) = x_foot;
    }

    return x_feet;
}


// Touch-Down (TD) Switching Surface -- checks individual legs
bool Dynamics::S_TD(Vector_12d x_feet, Leg_Idx leg_idx) 
{   
    // foot variable
    Vector_6d x_foot = x_feet.segment<6>(6*leg_idx);

    // unpack the state variables
    double pz_foot, vz_foot;
    pz_foot = x_foot(2);
    vz_foot = x_foot(5);

    // check the switching surface conditions
    bool gnd_pos, neg_vel, touchdown;
    gnd_pos = pz_foot <= 0.0;       // foot penetrated the ground
    neg_vel = vz_foot <= 0.0;       // foot is moving downward
    touchdown = gnd_pos && neg_vel; // if true, the foot touched down 

    return touchdown;
}


// Take-Off (TO) Switching Surface -- checks individual legs
bool Dynamics::S_TO(Vector_12d x_sys, Vector_12d x_legs, Vector_6d u, Leg_Idx leg_idx) 
{    
    // leg state
    Vector_6d x_leg = x_legs.segment<6>(6*leg_idx);

    // commanded leg state
    Vector_3d x_leg_command = x_sys.segment<3>(6 + 3*leg_idx);

    // leg input
    Vector_3d u_leg = u.segment<3>(3*leg_idx);

    // actual leg state variables
    double r, rdot;
    r = x_leg(0);
    rdot = x_leg(3);

    // commanded leg state variables
    double l0_command, l0dot_command;
    l0_command = x_leg_command(0);
    l0dot_command = u_leg(0);

    // check the switching surface condition (zero force in leg) 
    bool nom_length, pos_vel, takeoff;
    nom_length = r >= l0_command;     // leg is at or above the commanded length
    pos_vel = rdot >= l0dot_command;  // leg rate is exapnding the leg length
    takeoff = nom_length && pos_vel;  // if true, the leg took off

    return takeoff;
}


// Check if a switching event has occurred
Domain Dynamics::check_switching_event(const Vector_12d& x_sys, const Vector_12d& x_legs, const Vector_12d& x_feet, const Vector_6d& u, const Domain& d_current)
{
    // return the next domain after checking the switching surfaces
    Domain d_next(this->n_leg);

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {
        // contact state
        Contact c = d_current[i];

        // the i-th leg is in CONTACT
        if (c == Contact::SWING) {
            // check for a touchdown event
            bool touchdown = this->S_TD(x_feet, i);

            // if a touchdown event is detected, switch the leg to SWING
            d_next[i] = touchdown ? Contact::STANCE : Contact::SWING;
        }

        // the i-th leg is in STANCE
        else if (c == Contact::STANCE) {
            // check for a takeoff event
            bool takeoff = this->S_TO(x_sys, x_legs, u, i);

            // if a takeoff event is detected, switch the leg to STANCE
            d_next[i] = takeoff ? Contact::SWING : Contact::STANCE;
        }
    }

    return d_next;
}


// apply the reset map
void Dynamics::reset_map(Vector_12d& x_sys, Vector_12d& x_legs, Vector_12d& x_feet, Vector_6d& u, Domain d_prev, Domain d_next)
{
    // states to apply reset map to 
    Vector_12d x_sys_post;
    Vector_12d x_legs_post;
    Vector_12d x_feet_post;
    x_sys_post = x_sys;    // we will intialize as prev and modify for post as needed
    x_legs_post = x_legs;  // we will intialize as prev and modify for post as needed
    x_feet_post = x_feet;  // we will intialize as prev and modify for post as needed

    // foot positions post
    Vector_6d p_feet_post;
    for (int i = 0; i < this->n_leg; i++) {
        p_feet_post.segment<3>(3*i) = x_feet.segment<3>(6*i);
    }

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {
        // determine if there was a switch
        bool switched = d_prev[i] != d_next[i];

        // this leg did not switch
        if (switched == false) {
            // skip this leg, nothing to modify
            continue;
        }

        // this leg went through a switch
        else if (switched == true) {

            // temporary variables
            Vector_6d x_leg_post_;
            Vector_6d x_leg_i_post;
            Vector_6d x_foot_i_post;

            // the i-th leg is now in CONTACT
            if (d_prev[i] == Contact::SWING && d_next[i] == Contact::STANCE) {

                // unpack the state variables
                Vector_3d p_foot = p_feet_post.segment<3>(3*i);
                Vector_3d p_foot_i_post;

                // update the foot location (based on hueristic)
                p_foot_i_post << p_foot(0), p_foot(1), 0.0;
                p_feet_post.segment<3>(3*i) = p_foot_i_post;

                // update the system state
                x_leg_post_ = x_legs.segment<6>(6*i);
                x_sys_post(6 + 3*i) = x_leg_post_(0); // reset leg length command to the actual leg length at TD
                x_sys_post(7 + 3*i) = x_leg_post_(1); // reset leg x angle command to the actual leg angle at TD
                x_sys_post(8 + 3*i) = x_leg_post_(2); // reset leg y angle command to the actual leg angle at TD

                // update the leg state
                Vector_12d x_legs_ = this->compute_leg_state(x_sys_post, u, p_feet_post, d_next);
                x_leg_i_post = x_legs_.segment<6>(6*i);
                x_legs_post.segment<6>(6*i) = x_leg_i_post;

                // update the foot state
                Vector_12d x_feet_ = this->compute_foot_state(x_sys_post, x_legs_post, p_feet_post, d_next);
                x_foot_i_post = x_feet_.segment<6>(6*i);
                x_feet_post.segment<6>(6*i) = x_foot_i_post;
            }

            // the i-th leg is now in STANCE
            else if (d_prev[i] == Contact::STANCE && d_next[i] == Contact::SWING) {

                // update the system state
                x_leg_post_ = x_legs.segment<6>(6*i);
                x_sys_post(6 + 3*i) = x_leg_post_(0); // reset leg length command to the actual leg length at TD
                x_sys_post(7 + 3*i) = x_leg_post_(1); // reset leg angle command to the actual leg angle at TD
                x_sys_post(8 + 3*i) = x_leg_post_(2); // reset leg angle command to the actual leg angle at TD

                // leg input
                Vector_3d u_leg = u.segment<3>(3*i);

                // update the leg state
                x_leg_i_post = x_legs_post.segment<6>(6*i);
                x_leg_i_post(3) = u(0); // leg velocity is commanded velocity
                x_leg_i_post(4) = u(1); // leg angular x velocity is commanded angular velocity
                x_leg_i_post(5) = u(2); // leg angular y velocity is commanded angular velocity

                // update the foot state
                Vector_12d x_feet_ = this->compute_foot_state(x_sys_post, x_legs_post, p_feet_post, d_next);
                x_foot_i_post = x_feet_.segment<6>(6*i);
                x_feet_post.segment<6>(6*i) = x_foot_i_post;
            }
        }
    }

    // update the states
    x_sys = x_sys_post;
    x_legs = x_legs_post;
    x_feet = x_feet_post;
}


// interpolate an input signal
Vector_6d Dynamics::interpolate_control_input(double t, const Vector_1d_Traj& T_u, const Vector_6d_Traj& U) 
{
    // we want to find the interpolated control input
    Vector_6d u;

    // find the first element in the time vector that is greater than the current time
    auto it = std::upper_bound(T_u.begin(), T_u.end(), t);
    int idx = std::distance(T_u.begin(), it) - 1; // returns -1 if before the first element,
                                                  // returns N - 1 if after the last element

    // Zero-order hold (Z)
    if (this->params.interp == 'Z') {
        // contant control input
        u = U[idx];
    }

    // Linear Interpolation (L)
    else if (this->params.interp == 'L') {
        // beyond the last element
        if (idx == T_u.size() - 1) {
            u = U[idx];
        }
        // before the last element (shouldn't happen)
        else if (idx == -1) {
            u = U[0];
        }
        // within a time interval
        else {
            // find the time interval
            double t0 = T_u[idx];
            double tf = T_u[idx + 1];
            Vector_6d u0 = U[idx];
            Vector_6d uf = U[idx + 1];

            // linear interpolation
            u = u0 + (uf - u0) * (t - t0) / (tf - t0);
        }
    }

    return u;
}


// resize the solution bundle to the same as the time vector
void Dynamics::resizeSolution(Solution_3D& sol, const Vector_1d_Traj& T_x) {
    const int N = T_x.size();
    sol.x_sys_t.resize(N);
    sol.x_leg_t.resize(N);
    sol.x_foot_t.resize(N);
    sol.u_t.resize(N);
    sol.lambda_t.resize(N);
    sol.tau_t.resize(N);
    sol.domain_t.resize(N);
}


// RK3 Integration of the system dynamics
Solution_3D Dynamics::RK3_rollout(const Vector_1d_Traj& T_x, const Vector_1d_Traj& T_u, 
                                     const Vector_12d& x0_sys, const Vector_6d& p0_feet, const Domain& d0, 
                                     const Vector_6d_Traj& U) 
{
    // time integration parameters
    double dt = T_x[1] - T_x[0]; // WARNING: assumes uniform time steps
    int N = T_x.size();

    Solution_3D sol;
    this->resizeSolution(sol, T_x);

    // initial condition
    Vector_12d x0_legs;
    Vector_12d x0_feet;
    x0_legs = this->compute_leg_state(x0_sys, U[0], p0_feet, d0);
    x0_feet = this->compute_foot_state(x0_sys, x0_legs, p0_feet, d0);

    // first iteration just to get the initial leg forces and torques
    Vector_6d lambda0;
    Vector_6d tau0;
    Dynamics_Result_3D res = this->dynamics(x0_sys, U[0], p0_feet, d0);
    lambda0 = res.lambdas;
    tau0 = res.taus;

    // populate the initial conditions
    sol.x_sys_t[0] = x0_sys;
    sol.x_leg_t[0] = x0_legs;
    sol.x_foot_t[0] = x0_feet;
    sol.u_t[0] = U[0];
    sol.lambda_t[0] = lambda0;
    sol.tau_t[0] = tau0;
    sol.domain_t[0] = d0;

    // current state variables
    Vector_12d xk_sys = x0_sys;
    Vector_12d xk_legs = x0_legs;
    Vector_12d xk_feet = x0_feet;
    Vector_6d p_feet = p0_feet;
    Vector_6d lambdak = lambda0;
    Vector_6d tauk = tau0;
    Domain dk = d0;
    Domain dk_next;

    // ************************************* RK Integration *************************************
    // viability variable (for viability kernel)
    bool viability = true;

    // intermmediate times, inputs, and vector fields
    double tk, t1, t2, t3;
    Vector_6d u1, u2, u3;
    Vector_12d f1, f2, f3;
    Dynamics_Result_3D res1, res2, res3;

    // forward propagate the system dynamics
    for (int k = 1; k < N; k++) {

        // interpolation times
        tk = k * dt;
        t1 = tk;
        t2 = tk + 0.5 * dt;
        t3 = tk + dt;

        // interpolate the control input
        u1 = this->interpolate_control_input(t1, T_u, U);
        u2 = this->interpolate_control_input(t2, T_u, U);
        u3 = this->interpolate_control_input(t3, T_u, U);

        // vector fields for RK3 integration
        res1 = this->dynamics(xk_sys, 
                              u1, p_feet, dk);
        f1 = res1.xdot;
        res2 = this->dynamics(xk_sys + 0.5 * dt * f1,
                              u2, p_feet, dk);
        f2 = res2.xdot;
        res3 = this->dynamics(xk_sys - dt * f1 + 2.0 * dt * f2,
                              u3, p_feet, dk);
        f3 = res3.xdot;

        // take the RK3 step
        xk_sys = xk_sys + (dt / 6.0) * (f1 + 4.0 * f2 + f3);
        xk_legs = this->compute_leg_state(xk_sys, u3, p_feet, dk);
        xk_feet = this->compute_foot_state(xk_sys, xk_legs, p_feet, dk);

        // check for switching events
        dk_next = this->check_switching_event(xk_sys, xk_legs, xk_feet, u3, dk);

        // if there was a switching event, apply the reset map
        if (dk_next != dk) {

            // update all the states
            this->reset_map(xk_sys, xk_legs, xk_feet, u3, dk, dk_next);

            // update the foot positions
            for (int i = 0; i < this->n_leg; i++) {
                p_feet.segment<3>(3*i) = xk_feet.segment<3>(6*i);
            }

            // update the domain
            dk = dk_next;
        }

        // do a dynamics query to get the leg forces and torques
        res = this->dynamics(xk_sys, u3, p_feet, dk);
        lambdak = res.lambdas;
        tauk = res.taus;

        // store the states
        sol.x_sys_t[k] = xk_sys;
        sol.x_leg_t[k] = xk_legs;
        sol.x_foot_t[k] = xk_feet;
        sol.u_t[k] = u3;
        sol.lambda_t[k] = lambdak;
        sol.tau_t[k] = tauk;
        sol.domain_t[k] = dk_next;
    }

    // pack the solution into the solution struct
    sol.t = T_x;
    sol.viability = viability;

    return sol;
}