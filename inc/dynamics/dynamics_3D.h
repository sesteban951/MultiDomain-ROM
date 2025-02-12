#pragma once

// standard includes
#include <iostream>

// package includes
#include "yaml-cpp/yaml.h"

// custom includes
#include "../types/types_3D.h"

// class for system dynamics
class Dynamics
{
    public:
        
        // Constructor and Destructor
        Dynamics(YAML::Node config_file);  
        ~Dynamics(){};

        // NonLinear System dynamics, xdot = f(x, u, d)
        Dynamics_Result dynamics(const Vector_12d& x_sys, 
                                       const Vector_6d& u, 
                                       const Vector_6d& p_feet,
                                       const Domain& d);

        // compute the leg state
        Vector_12d compute_leg_state(const Vector_12d& x_sys, 
                                        const Vector_6d& u,
                                        const Vector_6d& p_feet, 
                                        const Domain& d);

        // compute foot state in world frame
        Vector_12d compute_foot_state(const Vector_12d& x_sys, 
                                         const Vector_12d& x_legs,
                                         const Vector_6d& p_feet,
                                         const Domain& d);

        // Switching Surfaces
        bool S_TD(Vector_12d x_feet, Leg_Idx leg_idx);
        bool S_TO(Vector_12d x_sys, Vector_12d x_legs, Vector_6d u, Leg_Idx leg_idx);

        // check if a switching event has occurred
        Domain check_switching_event(const Vector_12d& x_sys, 
                                        const Vector_12d& x_legs, 
                                        const Vector_12d& x_feet, 
                                        const Vector_6d& u,
                                        const Domain& d_current);

        // apply the reset map
        void reset_map(Vector_12d& x_sys, 
                          Vector_12d& x_legs, 
                          Vector_12d& x_feet, 
                          Vector_6d& u,
                          Vector_6d& p_feet,
                          Domain& d_prev, 
                          Domain d_next);

        // interpolate the input signal
        Vector_6d interpolate_control_input(double t, 
                                               const Vector_1d_Traj& T_u, 
                                               const Vector_6d_Traj& U);
        
        // resize the solution bundle to the same as the time vector
        void resizeSolution(Solution& sol, const Vector_1d_Traj& T_x);

        // RK forwaqrd propagation
        Solution RK3_rollout(const Vector_1d_Traj& T_x, 
                                const Vector_1d_Traj& T_u, 
                                const Vector_12d& x0_sys,
                                const Vector_6d& p0_feet,
                                const Domain& d0,
                                const Vector_6d_Traj& U);

        // System parameters
        SystemParams params;

        // fixed parameters
        const int n_leg = 2;   // numberof legs

};
