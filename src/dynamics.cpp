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
    else if (c == Contact::SWING) {
       
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
    Vector_4d x_leg = x_leg_[leg_idx];
    Vector_2d p_foot = p_feet[leg_idx];

    // Foot is in SWING
    if (c == Contact::SWING) {
        // com varaibles
        double px_com, pz_com, vx_com, vz_com;
        px_com = x_sys(0);
        pz_com = x_sys(1);
        vx_com = x_sys(2);
        vz_com = x_sys(3);

        // leg states
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
        x_foot << p_foot(0), p_foot(1), 0.0, 0.0;
    }

    return x_foot;
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
    double r, rdot, l0_command, thetadot_command;
    r = x_leg(0);
    rdot = x_leg(2);
    l0_command = x_sys(4 + 2*leg_idx);
    thetadot_command = x_sys(5 + 2*leg_idx); // TODO: would put this in the pos_vel boolean

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
    Domain d_next(2);

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
