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
}


// NonLinear Dynamics function, xdot = f(x, u, d)
Vector_8d Dynamics::dynamics(Vector_8d x, Vector_2d_List u, Vector_2d_List p_feet, Domain d) 
{
    // access some system parameters
    double m = this->params.m;
    double g = this->params.g;
    double k = this->params.k;
    double b = this->params.b;

    // unpack the state vector
    Vector_2d p_com, v_com;
    p_com << x(0), x(1);
    v_com << x(2), x(3);

    // vectors to use in the calculations
    Vector_2d a_com, g_vec, f_com; // vectors that affect COM dynamics 
    Vector_4d v_leg;      // velocity state of the legs
    double tau_ankle;     // ankle torque  
    Vector_8d xdot;       // state derivative
    f_com.setZero();
    tau_ankle = 0.0;

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {

        // The i'th leg is in NOT in contact, no forces or torques
        if (d[i] == Contact::SWING) {
            // in swing, no forces or torques on the COM
        }

        // The i'th leg is in contact, there are forces and torques
        else if (d[i] == Contact::STANCE) {

            // in stance, compute the forces and torques on the COM
            Vector_2d lambd_leg, f_com_ankle;
            
            // useful variables
            Vector_2d r_vec, rdot_vec, r_hat, p_foot;
            double r_norm, rdot_norm, r_x, r_z, rdot_x, rdot_z;
            double theta, thetadot;
            double l0_command, l0dot_command;

            // leg r vector
            r_vec = p_feet[i] - p_com;
            r_norm = r_vec.norm();
            r_hat = r_vec / r_norm;
            r_x = r_vec(0);
            r_z = r_vec(1);

            // leg rdot_vec
            rdot_vec = -v_com;
            rdot_x = rdot_vec(0);
            rdot_z = rdot_vec(1);
            rdot_norm = -v_com.dot(r_hat);

            // compute the force along the leg
            l0_command = x(4 + 2*i);
            l0dot_command = u[i][0]; // <-- TODO: add this to leg dynamics (next line)
            lambd_leg = -r_hat * (k * (l0_command - r_norm) - b * rdot_norm); // TODO: don't I need to add the input velocity to track?
                                                                              // But then you need to modify the take off switching manifold

            // TODO: think about how having torque on both legs when both are in stance. Doesn't make sense...
            // compute the ankle torque
            if (this->params.torque_ankle == true) {
                // actual angle state and commanded angle states
                double theta, theta_command, thetadot, thetadot_command;
                theta = -std::atan2(r_x, -r_z);
                thetadot = (r_z * rdot_x - r_x * rdot_z) / (r_norm * r_norm);
                theta_command = x(5 + 2*i);
                thetadot_command = u[i][1];

                // compute the ankle torque
                double kp = this->params.torque_ankle_kp;
                double kd = this->params.torque_ankle_kd;
                tau_ankle = kp * (theta_command - theta) + kd * (thetadot_command - thetadot);

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
            f_com += f_com_ankle + lambd_leg;
        }
    }

    // build the gravity vector
    g_vec << 0.0, -g;

    // compute the leg vlocities
    v_leg << u[0], u[1]; 

    // compute acceleration of the center of mass
    a_com << (1/m) * f_com + g_vec;

    // final form of the state derivative
    xdot << v_com, a_com, v_leg;

    // return the state derivative
    return xdot;
}


// compute the leg state
Vector_4d Dynamics::compute_leg_state(Vector_8d x_sys, Vector_2d_List p_feet, Vector_2d_List u_, Domain d, Leg_Idx leg_idx)
{
    // leg state variable
    Vector_4d x_leg;

    // unpack variables for leg index
    Contact c = d[leg_idx];
    Vector_2d u = u_[leg_idx];

    // indivudual states
    double r, theta, rdot, thetadot;

    // in flight (single integrator dynamics)
    if (c == Contact::SWING) {
       
        // leg positions are the integrated velocity commands
        r = x_sys(4 + 2*leg_idx);
        theta = x_sys(5 + 2*leg_idx);
        rdot = u(0);
        thetadot = u(1);

        // populate the polar state vector
        x_leg << r, theta, rdot, thetadot;
    }

    // on ground (leg state goverened by the COM dynamics)
    else if (c == Contact::STANCE) {
       
        // unpack COM state
        Vector_2d p_com, v_com;
        p_com << x_sys(0), x_sys(1);
        v_com << x_sys(2), x_sys(3);

        // get the foot position
        Vector_2d p_foot = p_feet[leg_idx];

        // leg vectors
        Vector_2d r_vec, r_hat, rdot_vec;
        double r_norm, r_x, r_z, rdot_x, rdot_z;
        r_vec = p_foot - p_com;
        r_norm = r_vec.norm();
        r_hat = r_vec / r_norm;
        rdot_vec = -v_com;

        r_x = r_vec(0);
        r_z = r_vec(1);
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

    return x_leg;
}


// compute foot state in world frame
Vector_4d Dynamics::compute_foot_state(Vector_8d x_sys, Vector_4d_List x_leg_, Vector_2d_List p_feet, Domain d, Leg_Idx leg_idx)
{
    // foot state variable
    Vector_4d x_foot;

    // unpack variables for leg index
    Contact c = d[leg_idx];

    // Foot is in SWING
    if (c == Contact::SWING) {

        // com varaibles
        double px_com, pz_com, vx_com, vz_com;
        px_com = x_sys(0);
        pz_com = x_sys(1);
        vx_com = x_sys(2);
        vz_com = x_sys(3);

        // leg states
        Vector_4d x_leg = x_leg_[leg_idx];
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
        Vector_2d p_foot = p_feet[leg_idx];
        x_foot << p_foot(0), p_foot(1), 0.0, 0.0;
    }

    return x_foot;
}


// compute the leg force
Vector_2d Dynamics::compute_leg_force(Vector_8d x_sys, Vector_4d_List x_leg_, Vector_2d_List p_feet, Vector_2d_List u_, Domain d, Leg_Idx leg_idx)
{
    // leg force to compute
    Vector_2d lambd;

    // contact
    Contact c = d[leg_idx];

    // In SWING, no force
    if (c == Contact::SWING) {
        lambd << 0.0, 0.0;
    }

    // In STANCE, compute the leg force
    else if (c == Contact::STANCE) {
        
        // com states
        Vector_2d p_com, v_com;
        p_com << x_sys(0), x_sys(1);
        v_com << x_sys(2), x_sys(3);

        // leg sates 
        Vector_4d x_leg = x_leg_[leg_idx];
        double r, rdot;
        r = x_leg(0);
        rdot = x_leg(2);

        // foot sates
        Vector_2d p_foot = p_feet[leg_idx];

        // compute the relevant vectors
        Vector_2d r_vec, r_hat, rdot_vec;
        r_vec = p_foot - p_com;
        r_hat = r_vec / r;

        // leg commands 
        double l0_command, l0dot_command;
        l0_command = x_sys(4 + 2*leg_idx);
        l0dot_command = u_[leg_idx](0); // TODO: add this in damping term

        // compute the leg force
        double k = this->params.k;
        double b = this->params.b;
        lambd = -r_hat * (k * (l0_command - r) - b * rdot); // TODO: add this in damping term
    }

    return lambd;
}


// Touch-Down (TD) Switching Surface -- checks individual legs
bool Dynamics::S_TD(Vector_4d x_foot) 
{   
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
bool Dynamics::S_TO(Vector_8d x_sys, Vector_4d x_leg, Leg_Idx leg_idx) 
{    
    // unpack the state variables
    double r, rdot, l0_command, l0dot_command;
    r = x_leg(0);
    rdot = x_leg(2);
    l0_command = x_sys(4 + 2*leg_idx);
    l0dot_command = x_sys(5 + 2*leg_idx); // TODO: would put this in the pos_vel boolean

    // check the switching surface condition (zero force in leg)
    bool nom_length, pos_vel, takeoff;
    nom_length = r >= l0_command;  // leg is at or above the commanded length
    pos_vel = rdot >= 0.0;         // leg is moving upward   // TODO: if I end up tracking the commanded velocity, need to modify this to (l0dot - rdot) <= 0.0
    takeoff = nom_length && pos_vel; // if true, the leg took off

    return takeoff;
}


// Check if a switching event has occurred
Domain Dynamics::check_switching_event(Vector_8d x_sys, Vector_4d_List x_leg, Vector_4d_List x_foot, Domain d_current)
{
    // return the next domain after checking the switching surfaces
    Domain d_next(this->n_leg);

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {
        // unpack the state variables
        Vector_4d x_foot_i = x_foot[i];
        Vector_4d x_leg_i = x_leg[i];
        Contact d_i = d_current[i];

        // the i-th leg is in CONTACT
        if (d_i == Contact::SWING) {
            // check for a touchdown event
            bool touchdown = this->S_TD(x_foot_i);

            // if a touchdown event is detected, switch the leg to SWING
            d_next[i] = touchdown ? Contact::STANCE : Contact::SWING;
        }

        // the i-th leg is in STANCE
        else if (d_i == Contact::STANCE) {
            // check for a takeoff event
            bool takeoff = this->S_TO(x_sys, x_leg_i, i);

            // if a takeoff event is detected, switch the leg to STANCE
            d_next[i] = takeoff ? Contact::SWING : Contact::STANCE;
        }
    }

    return d_next;
}


// apply the reset map
void Dynamics::reset_map(Vector_8d& x_sys, Vector_4d_List& x_leg, Vector_4d_List& x_foot, Vector_2d_List u, Domain d_prev, Domain d_next)
{
    // states to apply reset map to 
    Vector_8d x_sys_post;
    Vector_4d_List x_leg_post(this->n_leg);
    Vector_4d_List x_foot_post(this->n_leg);
    x_sys_post = x_sys;    // we will intialize as prev and modify for post as needed
    x_leg_post = x_leg;    // we will intialize as prev and modify for post as needed
    x_foot_post = x_foot;  // we will intialize as prev and modify for post as needed

    // useful intermmediate variables
    Vector_2d_List p_feet_post(this->n_leg);
    for (int i = 0; i < this->n_leg; ++i) {
        p_feet_post[i] << x_foot[i](0), 
                          x_foot[i](1); // we will initialize as prev and modify for post as needed
    }
    Vector_2d p_foot_i_post;

    // loop through the legs
    for (int i = 0; i < this->n_leg; i++)
    {
        // determine if there was a switch
        bool switched = d_prev[i] != d_next[i];

        // this leg did not switch
        if (switched == false) {
            // skip this tieration of the for lop
            continue;
        }

        // this leg went through a switch
        else if (switched == true) {
            // unpack the state variables
            Vector_4d x_foot_i = x_foot[i];
            Vector_4d x_leg_i = x_leg[i];
            Vector_2d p_foot_i = p_feet_post[i];
            Contact d_i_prev = d_prev[i];
            Contact d_i_next = d_next[i];

            // the i-th leg is now in CONTACT
            if (d_i_prev == Contact::SWING && d_i_next == Contact::STANCE) {
                
                // update the foot location (based on hueristic)
                p_foot_i_post << p_foot_i(0), 0.0;
                p_feet_post[i] = p_foot_i_post;

                // update the leg state
                x_leg_post[i] = this->compute_leg_state(x_sys_post, p_feet_post, u, d_next, i);

                // update the system state
                x_sys_post(4 + 2*i) = x_leg_post[i](0); // reset leg length command to the actual leg length at TD
                x_sys_post(5 + 2*i) = x_leg_post[i](1); // reset leg angle command to the actual leg angle at TD

                // update the foot state
                x_foot_post[i] << this->compute_foot_state(x_sys_post, x_leg_post, p_feet_post, d_next, i);
            }

            // the i-th leg is now in STANCE
            else if (d_i_prev == Contact::STANCE && d_i_next == Contact::SWING) {
                
                // update the foot location (based on hueristic) FIXME: probably dont need to update pos foot here
                // p_foot_i_post << p_foot_i(0), p_foot_i(1);
                // p_feet_post[i] = p_foot_i_post;

                // update the leg state
                x_leg_post[i](2) = u[i](0); // leg velocity is the commanded velocity
                x_leg_post[i](3) = u[i](1); // leg angular velocity is the commanded angular velocity

                // update the system state
                x_sys_post(4 + 2*i) = x_leg_post[i](0); // reset leg length command to the actual leg length at TD
                x_sys_post(5 + 2*i) = x_leg_post[i](1); // reset leg angle command to the actual leg angle at TD

                // update the foot state
                x_foot_post[i] = this->compute_foot_state(x_sys_post, x_leg_post, p_feet_post, d_next, i);
            }
        }
    }

    // update the states
    x_sys = x_sys_post;
    x_leg = x_leg_post;
    x_foot = x_foot_post;
}


// interpolate an input signal
Vector_2d_List Dynamics::interpolate_control_input(double t, Vector_1d_Traj T_u, Vector_2d_Traj U) 
{
    // we want to find the interpolated control input
    Vector_2d_List u(this->n_leg);

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
            Vector_2d_List u0 = U[idx];
            Vector_2d_List uf = U[idx + 1];

            // linear interpolation
            for (int i = 0; i < this->n_leg; i++) {
                u[i] = u0[i] + (uf[i] - u0[i]) * (t - t0) / (tf - t0);
            }
        }
    }

    return u;
}


// RK3 Integration of the system dynamics
Solution Dynamics::RK3_rollout(Vector_1d_Traj T_x, Vector_1d_Traj T_u, 
                              Vector_8d x0_sys, Vector_2d_List p0_feet, Domain d0, 
                              Vector_2d_Traj U)
{
    // time integration parameters
    double dt = T_x[1] - T_x[0];
    int N = T_x.size();

    // make the solutiuon trajectory containers
    Vector_8d_List x_sys_t(N);  // system state trajectory
    Vector_4d_Traj x_leg_t(N);  // leg state trajectory
    Vector_4d_Traj x_foot_t(N); // foot state trajectory
    Vector_2d_Traj u_t(N);      // interpolated control input trajectory
    Vector_2d_Traj lambd_t(N);   // leg force trajectory
    Domain_List domain_t(N);    // domain trajectory

    // initial condition
    Vector_4d_List x0_leg(this->n_leg);
    Vector_4d_List x0_foot(this->n_leg);
    Vector_2d_List lambd0(this->n_leg);
    for (int i = 0; i < this->n_leg; i++) {
        x0_leg[i] = this->compute_leg_state(x0_sys, p0_feet, U[0], d0, i);
        x0_foot[i] = this->compute_foot_state(x0_sys, x0_leg, p0_feet, d0, i);
        lambd0[i] = this->compute_leg_force(x0_sys, x0_leg, p0_feet, U[0], d0, i);
    }

    x_sys_t[0] = x0_sys;
    x_leg_t[0] = x0_leg;
    x_foot_t[0] = x0_foot;
    u_t[0] = U[0];
    lambd_t[0] = lambd0;
    domain_t[0] = d0;

    // current state variables
    Vector_8d xk_sys = x0_sys;
    Vector_4d_List xk_leg = x0_leg;
    Vector_4d_List xk_foot = x0_foot;
    Vector_2d_List p_feet = p0_feet;
    Vector_2d_List uk = U[0];
    Vector_2d_List lambdk = lambd0;
    Domain dk = d0;
    Domain dk_next;

    // ************************************* RK Integration *************************************
    // viability variable (for viability kernel)
    bool viability = true;

    // intermmediate times, inputs, and vector fields
    double tk, t1, t2, t3;
    Vector_2d_List u1, u2, u3;
    Vector_8d f1, f2, f3;

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
        f1 = this->dynamics(xk_sys, 
                            u1, p_feet, dk);
        f2 = this->dynamics(xk_sys + 0.5 * dt * f1,
                            u2, p_feet, dk);
        f3 = this->dynamics(xk_sys - dt * f1 + 2 * dt * f2,
                            u3, p_feet, dk);

        // take the RK3 step
        xk_sys = xk_sys + (dt / 6) * (f1 + 4 * f2 + f3);
        for (int i = 0; i < this->n_leg; i++) {
            xk_leg[i] = this->compute_leg_state(xk_sys, p_feet, u3, dk, i);
            xk_foot[i] = this->compute_foot_state(xk_sys, xk_leg, p_feet, dk, i);
        }

        // check for switching events
        dk_next = this->check_switching_event(xk_sys, xk_leg, xk_foot, dk);

        // if there was a switching event, apply the reset map
        if (dk_next != dk) {
            
            // update all the states
            this->reset_map(xk_sys, xk_leg, xk_foot, u3, dk, dk_next);

            // update the foot positions
            for (int i = 0; i < this->n_leg; i++) {
                p_feet[i] << xk_foot[i](0), 
                             xk_foot[i](1);
            }

            // update the domain
            dk = dk_next;
        }

        // compute the leg forces
        for (int i = 0; i < this->n_leg; i++) {
            lambdk[i] = this->compute_leg_force(xk_sys, xk_leg, p_feet, u3, dk, i);
        }

        // store the states
        x_sys_t[k] = xk_sys;
        x_leg_t[k] = xk_leg;
        x_foot_t[k] = xk_foot;
        u_t[k] = u3;
        lambd_t[k] = lambdk;
        domain_t[k] = dk_next;
    }

    // pack the solution into the solution struct
    Solution sol;
    sol.t = T_x;
    sol.x_sys_t = x_sys_t;
    sol.x_leg_t = x_leg_t;
    sol.x_foot_t = x_foot_t;
    sol.u_t = u_t;
    sol.lambd_t = lambd_t;
    sol.domain_t = domain_t;
    sol.viability = viability;

    return sol;
}