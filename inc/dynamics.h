#ifndef DYNAMICS_H
#define DYNAMICS_H

// standard includes
#include <iostream>

// package includes
#include "yaml-cpp/yaml.h"

// custom includes
#include "types.h"

// class for system dynamics
class Dynamics
{
    public:
        
        // Constructor and Destructor
        Dynamics(YAML::Node config_file);  
        ~Dynamics(){};

        // NonLinear System dynamics, xdot = f(x, u, d)
        Vector_8d dynamics(Vector_8d x, 
                           Vector_2d_List u, 
                           Vector_2d_List p_feet,
                           Domain d);

        // compute the leg state
        Vector_4d compute_leg_state(Vector_8d x_sys, 
                                    Vector_2d_List p_feet, 
                                    Vector_2d_List u, 
                                    Domain d, 
                                    Leg_Idx leg_idx);

        // compute foot state in world frame
        Vector_4d compute_foot_state(Vector_8d x_sys, 
                                     Vector_4d_List x_leg_, 
                                     Vector_2d_List p_feet, 
                                     Domain d, 
                                     Leg_Idx leg_idx);
        
        // Switching Surfaces
        bool S_TD(Vector_4d x_foot);
        bool S_TO(Vector_8d x_sys, Vector_4d x_leg, Leg_Idx leg_idx);

        // check if a switching event has occurred
        Domain check_switching_event(Vector_8d x_sys, 
                                     Vector_4d_List x_leg, 
                                     Vector_4d_List x_foot, 
                                     Domain d);

        // apply the reset map
        void reset_map(Vector_8d& x_sys, 
                      Vector_4d_List& x_leg, 
                      Vector_4d_List& x_foot, 
                      Vector_2d_List u, 
                      Domain d_prev, 
                      Domain d_next);

        // interpolate the input signal
        Vector_2d_List interpolate_control_input(double t, 
                                                Vector_1d_Traj T_u, 
                                                Vector_2d_Traj U);
        
        // RK forwaqrd propagation
        Solution RK3_rollout(Vector_1d_Traj T_x, 
                            Vector_1d_Traj T_u, 
                            Vector_8d x0_sys, 
                            Vector_2d_List p0_feet, 
                            Domain d0, 
                            Vector_2d_Traj U);

        // System parameters
        SystemParams params;

        // fixed parameters
        const int n_leg = 2;      // vector size of legs
};

#endif
