#pragma once

// standard includes
#include <iostream>
#include <cmath>

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
        Dynamics_Result dynamics(const Vector_8d& x_sys, 
                                 const Vector_4d& u, 
                                 const Vector_4d& p_feet,
                                 const Domain& d);
        Dynamics_Result_3D dynamics_3D(const Vector_12d& x_sys, 
                                       const Vector_6d& u, 
                                       const Vector_6d& p_feet,
                                       const Domain& d);

        // compute the leg state
        Vector_8d compute_leg_state(const Vector_8d& x_sys, 
                                    const Vector_4d& u,
                                    const Vector_4d& p_feet, 
                                    const Domain& d);
        Vector_12d compute_leg_state_3D(const Vector_12d& x_sys, 
                                        const Vector_6d& u,
                                        const Vector_6d& p_feet, 
                                        const Domain& d);

        // compute foot state in world frame
        Vector_8d compute_foot_state(const Vector_8d& x_sys, 
                                     const Vector_8d& x_legs,
                                     const Vector_4d& p_feet,
                                     const Domain& d);
        Vector_12d compute_foot_state_3D(const Vector_12d& x_sys, 
                                         const Vector_12d& x_legs,
                                         const Vector_6d& p_feet,
                                         const Domain& d);

        // Switching Surfaces
        bool S_TD(Vector_8d x_feet, Leg_Idx leg_idx);
        bool S_TO(Vector_8d x_sys, Vector_8d x_legs, Vector_4d u, Leg_Idx leg_idx);
        bool S_TD_3D(Vector_12d x_feet, Leg_Idx leg_idx);
        bool S_TO_3D(Vector_12d x_sys, Vector_12d x_legs, Vector_6d u, Leg_Idx leg_idx);

        // check if a switching event has occurred
        Domain check_switching_event(const Vector_8d& x_sys, 
                                     const Vector_8d& x_legs, 
                                     const Vector_8d& x_feet, 
                                     const Vector_4d& u,
                                     const Domain& d);
        Domain check_switching_event_3D(const Vector_12d& x_sys, 
                                        const Vector_12d& x_legs, 
                                        const Vector_12d& x_feet, 
                                        const Vector_6d& u,
                                        const Domain& d_current);

        // apply the reset map
        void reset_map(Vector_8d& x_sys, 
                       Vector_8d& x_legs, 
                       Vector_8d& x_feet, 
                       Vector_4d& u,
                       Domain d_prev, 
                       Domain d_next);
        void reset_map_3D(Vector_12d& x_sys, 
                          Vector_12d& x_legs, 
                          Vector_12d& x_feet, 
                          Vector_6d& u,
                          Domain d_prev, 
                          Domain d_next);

        // interpolate the input signal
        Vector_4d interpolate_control_input(double t, 
                                            const Vector_1d_Traj& T_u, 
                                            const Vector_4d_Traj& U);
        Vector_6d interpolate_control_input_3D(double t, 
                                               const Vector_1d_Traj& T_u, 
                                               const Vector_6d_Traj& U);
        
        // resize the solution bundle to the same as the time vector
        void resizeSolution(Solution& sol, const Vector_1d_Traj& T_x);
        void resizeSolution_3D(Solution_3D& sol, const Vector_1d_Traj& T_x);

        // RK forwaqrd propagation
        Solution RK3_rollout(const Vector_1d_Traj& T_x, 
                             const Vector_1d_Traj& T_u, 
                             const Vector_8d& x0_sys,
                             const Vector_4d& p0_feet,
                             const Domain& d0,
                             const Vector_4d_Traj& U);
        Solution_3D RK3_rollout_3D(const Vector_1d_Traj& T_x, 
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
