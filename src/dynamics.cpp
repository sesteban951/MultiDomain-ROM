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