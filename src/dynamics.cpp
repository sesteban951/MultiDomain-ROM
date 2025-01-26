#include "../inc/dynamics.h"


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

    // build the linear dynamics matrices for leg command dynamics
    Matrix_8d A;
    Matrix_8x4d B;
    A.setZero();
    B.setZero();
    A.block<2,2>(0,2) = Eigen::Matrix2d::Identity();
    A.block<2,2>(4,6) = Eigen::Matrix2d::Identity();
    B.block<2,2>(2,0) = Eigen::Matrix2d::Identity();
    B.block<2,2>(6,2) = Eigen::Matrix2d::Identity();
    this->A = A;
    this->B = B;
}


// NonLinear Dynamics function, xdot = f(x, u, d)
Vector_12d Dynamics::dynamics(Vector_12d x_sys, Vector_4d u, Vector_4d p_feet, Domain d) 
{
    // want to compute the state derivative
    Vector_12d xdot;       // state derivative, [xdot_com, xdot_legs]

    // access some system parameters
    double m = this->params.m;
    double g = this->params.g;
    double k = this->params.k;
    double b = this->params.b;

    // linear dynamics matrices for the leg commands
    Matrix_8d A = this->A;
    Matrix_8x4d B = this->B;

    // unpack the state vector
    Vector_4d x_com;
    x_com = x_sys.segment<4>(0);

    // vectors to use in the calculations
    Vector_2d v_com, a_com, g_vec, f_com; // vectors that affect COM dynamics 
    Vector_8d x_legs_command, xdot_legs_command;   // state, state derivative leg
    double tau_ankle;              // ankle torque  
    f_com.setZero();
    tau_ankle = 0.0;
    v_com = x_com.segment<2>(2);

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
        }

        // The i'th leg is in contact, there are forces and torques
        else if (c == Contact::STANCE) {

            // in stance, compute the forces and torques on the COM
            Vector_2d lambd_leg, f_com_ankle;

            // COM state
            Vector_2d p_com;
            p_com = x_com.segment<2>(0);

            // LEG Command State
            Vector_4d x_leg_command;
            x_leg_command = x_sys.segment<4>(4 + 4*i);

            // foot position
            Vector_2d p_foot;
            p_foot = p_feet.segment<2>(2*i);
            
            // useful variables
            Vector_2d r_vec, rdot_vec, r_hat;
            double r_norm, rdot_norm, r_x, r_z, rdot_x, rdot_z;
            double l0_command, l0dot_command;

            // leg r vector
            r_vec = p_foot - p_com;
            r_norm = r_vec.norm();
            r_hat = r_vec / r_norm;
            r_x = r_vec(0);
            r_z = r_vec(1);

            // leg rdot_vec
            rdot_vec = -v_com;    // assumes no sliding feet
            rdot_x = rdot_vec(0);
            rdot_z = rdot_vec(1);
            rdot_norm = -v_com.dot(r_hat);

            // compute the force along the leg
            l0_command = x_leg_command(0);
            l0dot_command = x_leg_command(2);
            lambd_leg = -r_hat * (k * (l0_command - r_norm) - b * (l0dot_command - rdot_norm)); // TODO: add acceleration feed forward?
                                                   

            // TODO: think about how having torque on both legs when both are in stance -- overactuated
            // compute the ankle torque
            if (this->params.torque_ankle == true) {

                // actual angle state and commanded angle states
                double theta, theta_command, thetadot, thetadot_command;
                theta = -std::atan2(r_x, -r_z);
                thetadot = (r_z * rdot_x - r_x * rdot_z) / (r_norm * r_norm);
                theta_command = x_leg_command(1);
                thetadot_command = x_leg_command(3);

                // compute the ankle torque
                double kp = this->params.torque_ankle_kp;
                double kd = this->params.torque_ankle_kd;
                tau_ankle = kp * (theta_command - theta) + kd * (thetadot_command - thetadot); // TODO: add acceleration feed forward?

                // convert leg torque into COM force perpendicular to the leg
                Vector_2d f_unit;
                double f_mag = tau_ankle / r_norm;
                f_unit << std::cos(theta), -std::sin(theta);
                f_com_ankle = f_mag * f_unit;
            }
            // no ankle torque
            else {
                tau_ankle = 0.0;
                f_com_ankle << 0.0, 0.0;
            }

            // add the leg force and torque to the COM dynamics
            f_com += lambd_leg + f_com_ankle;
        }
    }

    // build the gravity vector
    g_vec << 0.0, -g;

    // compute acceleration of the center of mass
    a_com << (1/m) * f_com + g_vec;

    // compute the leg vlocities
    x_legs_command = x_sys.segment<8>(4);
    xdot_legs_command = A * x_legs_command + B * u;

    // final form of the state derivative
    xdot << v_com, a_com, xdot_legs_command;

    // return the state derivative
    return xdot;
}


// compute the leg state
Vector_8d Dynamics::compute_leg_state(Vector_12d x_sys, Vector_4d p_feet, Domain d)
{
    // want to compute the leg states
    Vector_8d x_legs;
    Vector_4d x_leg;

    // indivudual states
    double r, theta, rdot, thetadot;

    // contact state
    Contact c;

    // loop though legs and compute the leg states
    for (int i = 0; i < this->n_leg; i++)
    {
        // contact state
        c = d[i];

        // in swing (double integrator dynamics)
        if (c == Contact::SWING) {
        
            // leg command vector
            Vector_4d x_leg_command = x_sys.segment<4>(4 + 4*i);

            // double integrator state
            r = x_leg_command(0);
            theta = x_leg_command(1);
            rdot = x_leg_command(2);
            thetadot = x_leg_command(3);

            // populate the polar state vector
            x_leg << r, theta, rdot, thetadot;
        }

        // in stance (leg state goverened by the COM dynamics)
        else if (c == Contact::STANCE) {

            // COM state
            Vector_2d p_com, v_com;
            p_com = x_sys.segment<2>(0);
            v_com = x_sys.segment<2>(2);

            // get the foot position
            Vector_2d p_foot = p_feet.segment<2>(2*i);

            // leg vectors
            Vector_2d r_vec, r_hat, rdot_vec;
            double r_norm, r_x, r_z, rdot_x, rdot_z;
            r_vec = p_foot - p_com;
            r_norm = r_vec.norm();
            r_hat = r_vec / r_norm;
            r_x = r_vec(0);
            r_z = r_vec(1);

            rdot_vec = -v_com; // assumes no sliding feet
            rdot_x = rdot_vec(0);
            rdot_z = rdot_vec(1);

            // compute the polar leg state 
            r = r_norm;
            rdot = -v_com.dot(r_hat);
            theta = -std::atan2(r_x, -r_z);
            thetadot = (r_z * rdot_x - r_x * rdot_z) / (r * r);

            // populate the polar state vector
            x_leg << r, theta, rdot, thetadot;
        }

        // populate the leg state vector
        x_legs.segment<4>(4*i) = x_leg;
    }

    return x_legs;
}


// compute foot state in world frame
Vector_8d Dynamics::compute_foot_state(Vector_12d x_sys, Vector_8d x_legs, Vector_4d p_feet, Domain d)
{
    // want to compute the foot states
    Vector_8d x_feet;
    Vector_4d x_foot;

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
            double px_com, pz_com, vx_com, vz_com;
            px_com = x_sys(0);
            pz_com = x_sys(1);
            vx_com = x_sys(2);
            vz_com = x_sys(3);

            // leg states
            Vector_4d x_leg = x_legs.segment<4>(4*i);
            double r, theta, rdot, thetadot;
            r = x_leg(0);
            theta = x_leg(1);
            rdot = x_leg(2);
            thetadot = x_leg(3);

            // compute the foot state via fwd kinematics
            double px_foot, pz_foot, vx_foot, vz_foot;
            px_foot = px_com - r * std::sin(theta);
            pz_foot = pz_com - r * std::cos(theta);
            vx_foot = vx_com - rdot * std::sin(theta) - r * thetadot * std::cos(theta);
            vz_foot = vz_com - rdot * std::cos(theta) + r * thetadot * std::sin(theta);

            // populate the foot state vector
            x_foot << px_foot, pz_foot, vx_foot, vz_foot;
        }

        // Foot is in STANCE
        else if (c == Contact::STANCE) {
            // foot state is the same as the foot position w/ zero velocity
            Vector_2d p_foot = p_feet.segment<2>(2*i);
            x_foot << p_foot(0), p_foot(1), 0.0, 0.0;
        }

        // insert into the foot state vector
        x_feet.segment<4>(4*i) = x_foot;
    }

    return x_feet;
}


// compute the leg force
Vector_4d Dynamics::compute_leg_force(Vector_12d x_sys, Vector_8d x_legs, Vector_4d p_feet, Vector_4d u, Domain d)
{
    // want to compute the leg force
    Vector_4d lambdas;
    Vector_2d lambd;

    // contact state
    Contact c;

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {
        // contact state
        c = d[i];

        // In SWING, no force
        if (c == Contact::SWING) {
            lambd << 0.0, 0.0;
        }

        // In STANCE, compute the leg force
        else if (c == Contact::STANCE) {
            
            // com states
            Vector_2d p_com, v_com;
            p_com << x_sys.segment<2>(0);
            v_com << x_sys.segment<2>(2);

            // leg sates 
            Vector_4d x_leg = x_legs.segment<4>(4*i);
            double r, rdot;
            r = x_leg(0);
            rdot = x_leg(2);

            // foot sates
            Vector_2d p_foot = p_feet.segment<2>(2*i);

            // compute the relevant vectors
            Vector_2d r_vec, r_hat, rdot_vec;
            r_vec = p_foot - p_com;
            r_hat = r_vec / r;

            // leg commands 
            Vector_4d x_leg_command = x_sys.segment<4>(4 + 4*i);
            double l0_command, l0dot_command;
            l0_command = x_leg_command(0);
            l0dot_command = x_leg_command(2);
            
            // compute the leg force
            double k = this->params.k;
            double b = this->params.b;
            lambd = -r_hat * (k * (l0_command - r) - b * (l0dot_command - rdot)); 
        }

        // insert into the leg force vector
        lambdas.segment<2>(2*i) = lambd;
    }

    return lambdas;
}


// Touch-Down (TD) Switching Surface -- checks individual legs
bool Dynamics::S_TD(Vector_8d x_feet, Leg_Idx leg_idx) 
{   
    // foot variable
    Vector_4d x_foot = x_feet.segment<4>(4*leg_idx);

    // unpack the state variables
    double pz_foot, vz_foot;
    pz_foot = x_foot(1);
    vz_foot = x_foot(3);

    // check the switching surface conditions
    bool gnd_pos, neg_vel, touchdown;
    gnd_pos = pz_foot <= 0.0;       // foot penetrated the ground
    neg_vel = vz_foot <= 0.0;       // foot is moving downward
    touchdown = gnd_pos && neg_vel; // if true, the foot touched down 

    return touchdown;
}


// Take-Off (TO) Switching Surface -- checks individual legs
bool Dynamics::S_TO(Vector_12d x_sys, Vector_8d x_legs, Leg_Idx leg_idx) 
{    
    // leg state
    Vector_4d x_leg = x_legs.segment<4>(4*leg_idx);

    // commanded leg state
    Vector_4d x_leg_command = x_sys.segment<4>(4 + 4*leg_idx);

    // actual leg state variables
    double r, rdot;
    r = x_leg(0);
    rdot = x_leg(2);

    // commanded leg state variables
    double l0_command, l0dot_command;
    l0_command = x_leg_command(0);
    l0dot_command = x_leg_command(2);

    // check the switching surface condition (zero force in leg) TODO: do I want to add acceleration feed forward? Would need the input
    bool nom_length, pos_vel, takeoff;
    nom_length = r >= l0_command;     // leg is at or above the commanded length
    pos_vel = rdot >= l0dot_command;  // leg is moving upward 
    takeoff = nom_length && pos_vel;  // if true, the leg took off

    return takeoff;
}


// Check if a switching event has occurred
Domain Dynamics::check_switching_event(Vector_12d x_sys, Vector_8d x_legs, Vector_8d x_feet, Domain d_current)
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
            bool takeoff = this->S_TO(x_sys, x_legs, i);

            // if a takeoff event is detected, switch the leg to STANCE
            d_next[i] = takeoff ? Contact::SWING : Contact::STANCE;
        }
    }

    return d_next;
}


// apply the reset map
void Dynamics::reset_map(Vector_12d& x_sys, Vector_8d& x_legs, Vector_8d& x_feet, Domain d_prev, Domain d_next)
{
    // states to apply reset map to 
    Vector_12d x_sys_post;
    Vector_8d x_legs_post;
    Vector_8d x_feet_post;
    x_sys_post = x_sys;    // we will intialize as prev and modify for post as needed
    x_legs_post = x_legs;  // we will intialize as prev and modify for post as needed
    x_feet_post = x_feet;  // we will intialize as prev and modify for post as needed

    // foot positions post
    Vector_4d p_feet_post;
    for (int i = 0; i < this->n_leg; i++) {
        p_feet_post.segment<2>(2*i) = x_feet.segment<2>(4*i);
    }

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {
        // determine if there was a switch
        bool switched = d_prev[i] != d_next[i];

        // this leg did not switch
        if (switched == false) {
            // skip this tieration of the for loop
            continue;
        }

        // this leg went through a switch
        else if (switched == true) {
            // unpack the state variables
            Vector_2d p_foot = p_feet_post.segment<2>(2*i);

            // temporary variables
            Vector_4d x_leg_i_post;
            Vector_4d x_foot_i_post;
            Vector_2d p_foot_i_post;

            // the i-th leg is now in CONTACT
            if (d_prev[i] == Contact::SWING && d_next[i] == Contact::STANCE) {
                
                // update the foot location (based on hueristic)
                p_foot_i_post << p_foot(0), 0.0;
                p_feet_post.segment<2>(2*i) = p_foot_i_post;

                // update the leg state
                Vector_8d x_legs_ = this->compute_leg_state(x_sys_post, p_feet_post, d_next);
                x_leg_i_post = x_legs_.segment<4>(4*i);
                x_legs_post.segment<4>(4*i) = x_leg_i_post;

                // update the system state
                x_sys_post.segment<4>(4 + 4*i) = x_leg_i_post; // TODO: decide if I want to reset the whole command
                                                               // set veloicty to zero and keep position?

                // update the foot state
                Vector_8d x_feet_ = this->compute_foot_state(x_sys_post, x_legs_post, p_feet_post, d_next);
                x_foot_i_post = x_feet_.segment<4>(4*i);
                x_feet_post.segment<4>(4*i) = x_foot_i_post;
            }

            // the i-th leg is now in STANCE
            else if (d_prev[i] == Contact::STANCE && d_next[i] == Contact::SWING) {

                // update the leg state
                x_leg_i_post = x_legs_post.segment<4>(4*i);
                
                // update the system state
                x_sys_post.segment<4>(4 + 4*i) = x_leg_i_post; // TODO: decide if I want to reset the whole command
                                                               // set veloicty to zero and keep position?

                // update the foot state
                Vector_8d x_feet_ = this->compute_foot_state(x_sys_post, x_legs_post, p_feet_post, d_next);
                x_foot_i_post = x_feet_.segment<4>(4*i);
            }
        }
    }

    // update the states
    x_sys = x_sys_post;
    x_legs = x_legs_post;
    x_feet = x_feet_post;
}


// // interpolate an input signal
// Vector_2d_List Dynamics::interpolate_control_input(double t, Vector_1d_Traj T_u, Vector_2d_Traj U) 
// {
//     // we want to find the interpolated control input
//     Vector_2d_List u(this->n_leg);

//     // find the first element in the time vector that is greater than the current time
//     auto it = std::upper_bound(T_u.begin(), T_u.end(), t);
//     int idx = std::distance(T_u.begin(), it) - 1; // returns -1 if before the first element,
//                                                   // returns N - 1 if after the last element

//     // Zero-order hold (Z)
//     if (this->params.interp == 'Z') {
//         // contant control input
//         u = U[idx];
//     }

//     // Linear Interpolation (L)
//     else if (this->params.interp == 'L') {
//         // beyond the last element
//         if (idx == T_u.size() - 1) {
//             u = U[idx];
//         }
//         // before the last element (shouldn't happen)
//         else if (idx == -1) {
//             u = U[0];
//         }
//         // within a time interval
//         else {
//             // find the time interval
//             double t0 = T_u[idx];
//             double tf = T_u[idx + 1];
//             Vector_2d_List u0 = U[idx];
//             Vector_2d_List uf = U[idx + 1];

//             // linear interpolation
//             for (int i = 0; i < this->n_leg; i++) {
//                 u[i] = u0[i] + (uf[i] - u0[i]) * (t - t0) / (tf - t0);
//             }
//         }
//     }

//     return u;
// }


// RK3 Integration of the system dynamics
Solution Dynamics::RK3_rollout(Vector_1d_Traj T_x, Vector_1d_Traj T_u, 
                               Vector_12d x0_sys, Vector_4d p0_feet, Domain d0, 
                               Vector_4d_Traj U) 
{
    // time integration parameters
    double dt = T_x[1] - T_x[0]; // CAUTION: assumes uniform time steps
    int N = T_x.size();

    // make the solutiuon trajectory containers
    Vector_12d_Traj x_sys_t(N);  // system state trajectory
    Vector_8d_Traj x_leg_t(N);   // leg state trajectory
    Vector_8d_Traj x_foot_t(N);  // foot state trajectory
    Vector_4d_Traj u_t(N);       // interpolated control input trajectory
    Vector_4d_Traj lambda_t(N);  // leg force trajectory
    Domain_Traj domain_t(N);     // domain trajectory

    // initial condition
    Vector_8d x0_legs;
    Vector_8d x0_feet;
    Vector_4d lambda0;
    x0_legs = this->compute_leg_state(x0_sys, p0_feet, d0);
    x0_feet = this->compute_foot_state(x0_sys, x0_legs, p0_feet, d0);
    lambda0 = this->compute_leg_force(x0_sys, x0_legs, p0_feet, U[0], d0);

    // populate the initial conditions
    x_sys_t[0] = x0_sys;
    x_leg_t[0] = x0_legs;
    x_foot_t[0] = x0_feet;
    u_t[0] = U[0];
    lambda_t[0] = lambda0;
    domain_t[0] = d0;

    // current state variables
    Vector_12d xk_sys = x0_sys;
    Vector_8d xk_legs = x0_legs;
    Vector_8d xk_feet = x0_feet;
    Vector_4d p_feet = p0_feet;
    Vector_4d lambdak = lambda0;
    Domain dk = d0;
    Domain dk_next;

    // ************************************* RK Integration *************************************
    // viability variable (for viability kernel)
    bool viability = true;

    // intermmediate times, inputs, and vector fields
    double tk, t1, t2, t3;
    Vector_4d u1, u2, u3;
    Vector_12d f1, f2, f3;

    // forward propagate the system dynamics
    for (int k = 1; k < N; k++) {

        // interpolation times
        tk = k * dt;
        t1 = tk;
        t2 = tk + 0.5 * dt;
        t3 = tk + dt;

        // interpolate the control input
        // u1 = this->interpolate_control_input(t1, T_u, U);
        // u2 = this->interpolate_control_input(t2, T_u, U);
        // u3 = this->interpolate_control_input(t3, T_u, U);
        u1 = U[0];
        u2 = U[0];
        u3 = U[0];

        // vector fields for RK3 integration
        f1 = this->dynamics(xk_sys, 
                            u1, p_feet, dk);
        f2 = this->dynamics(xk_sys + 0.5 * dt * f1,
                            u2, p_feet, dk);
        f3 = this->dynamics(xk_sys - dt * f1 + 2 * dt * f2,
                            u3, p_feet, dk);

        // take the RK3 step
        xk_sys = xk_sys + (dt / 6) * (f1 + 4 * f2 + f3);
        xk_legs = this->compute_leg_state(xk_sys, p_feet, dk);
        xk_feet = this->compute_foot_state(xk_sys, xk_legs, p_feet, dk);

        // check for switching events
        dk_next = this->check_switching_event(xk_sys, xk_legs, xk_feet, dk);

        // if there was a switching event, apply the reset map
        if (dk_next != dk) {
            
            // update all the states
            this->reset_map(xk_sys, xk_legs, xk_feet, dk, dk_next);

            // update the foot positions
            for (int i = 0; i < this->n_leg; i++) {
                p_feet.segment<2>(2*i) = xk_feet.segment<2>(4*i);
            }

            // update the domain
            dk = dk_next;
        }

        // compute the leg forces
        lambdak = this->compute_leg_force(xk_sys, xk_legs, p_feet, u3, dk);

        // store the states
        x_sys_t[k] = xk_sys;
        x_leg_t[k] = xk_legs;
        x_foot_t[k] = xk_feet;
        u_t[k] = u3;
        lambda_t[k] = lambdak;
        domain_t[k] = dk_next;
    }

    // pack the solution into the solution struct
    Solution sol;
    sol.t = T_x;
    sol.x_sys_t = x_sys_t;
    sol.x_leg_t = x_leg_t;
    sol.x_foot_t = x_foot_t;
    sol.u_t = u_t;
    sol.lambda_t = lambda_t;
    sol.domain_t = domain_t;
    sol.viability = viability;

    return sol;
}