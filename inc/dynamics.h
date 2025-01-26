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
        Vector_12d dynamics(Vector_12d x, 
                           Vector_4d u, 
                           Vector_4d p_feet,
                           Domain d);

        // compute the leg state
        Vector_8d compute_leg_state(Vector_12d x_sys, 
                                    Vector_4d p_feet, 
                                    Domain d);

        // compute foot state in world frame
        Vector_8d compute_foot_state(Vector_12d x_sys, 
                                     Vector_8d x_legs,
                                     Vector_4d p_feet,
                                     Domain d);
        
        // compute the leg force
        Vector_4d compute_leg_force(Vector_12d x_sys, 
                                   Vector_8d x_legs, 
                                   Vector_4d p_feet, 
                                   Vector_4d u, 
                                   Domain d);
        
        // Switching Surfaces
        bool S_TD(Vector_8d x_feet, Leg_Idx leg_idx);
        bool S_TO(Vector_12d x_sys, Vector_8d x_legs, Leg_Idx leg_idx);

        // check if a switching event has occurred
        Domain check_switching_event(Vector_12d x_sys, 
                                     Vector_8d x_legs, 
                                     Vector_8d x_feet, 
                                     Domain d);

        // apply the reset map
        void reset_map(Vector_12d& x_sys, 
                       Vector_8d& x_legs, 
                       Vector_8d& x_feet, 
                       Domain d_prev, 
                       Domain d_next);

        // interpolate the input signal
        Vector_4d interpolate_control_input(double t, 
                                            Vector_1d_Traj T_u, 
                                            Vector_4d_Traj U);
        
        // RK forwaqrd propagation
        Solution RK3_rollout(Vector_1d_Traj T_x, 
                             Vector_1d_Traj T_u, 
                             Vector_12d x0_sys,
                             Vector_4d p0_feet,
                             Domain d0,
                             Vector_4d_Traj U);

        // System parameters
        SystemParams params;

        // fixed parameters
        const int n_leg = 2;      // vector size of legs

        // linear dynamics matrices
        Matrix_8d A;
        Matrix_8x4d B;
    
};

#endif
