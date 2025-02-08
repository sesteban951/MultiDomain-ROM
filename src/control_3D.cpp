#include "../inc/control_3D.h"


// Constructor
Controller::Controller(YAML::Node config_file) : dynamics(config_file)
{
    // set the control parameters
    this->params.N_x = config_file["CTRL_PARAMS"]["N_x"].as<int>();
    this->params.N_u = config_file["CTRL_PARAMS"]["N_u"].as<int>();
    this->params.dt_x = config_file["CTRL_PARAMS"]["dt_x"].as<double>();
    this->params.K = config_file["CTRL_PARAMS"]["K"].as<int>();
    this->params.N_elite = config_file["CTRL_PARAMS"]["N_elite"].as<int>();
    this->params.CEM_iters = config_file["CTRL_PARAMS"]["CEM_iters"].as<int>();
    
    // ensure the parameters make sense
    if (this->params.N_x < this->params.N_u) {
        throw std::runtime_error("N_x must be greater than N_u!");
    }
    if (this->params.N_elite > this->params.K) {
        throw std::runtime_error("N_elite must be less than or equal to K!");
    }

    // compute u(t) dt (N_x of the integration is not necessarily equal to the number of control points, N_u)
    double T = (this->params.N_x-1) * this->params.dt_x;
    this->params.dt_u = T / (this->params.N_u-1);

    // generate the time arrays
    Vector_1d_Traj T_x;
    Vector_1d_Traj T_u;
    T_x.resize(this->params.N_x);
    T_u.resize(this->params.N_u);
    for (int i = 0; i < this->params.N_x; i++) {
        T_x[i] = i * this->params.dt_x;
    }
    for (int i = 0; i < this->params.N_u; i++) {
        T_u[i] = i * this->params.dt_u;
    }
    this->params.T_x = T_x;
    this->params.T_u = T_u;

    // initialize the class variables
    this->initialize_variables();

    // initialize the cost matrices
    this->initialize_costs(config_file);

    // construct the reference trajectory
    this->initialize_reference_trajectories(config_file);

    // construct the initial distribution
    this->initialize_distribution(config_file);

    // set the number of parallel threads to use
    this->threading_enabled = config_file["THREADING"]["enabled"].as<bool>();
    if (threading_enabled == true) {

        // set the number of threads
        int num_threads = config_file["THREADING"]["num_threads"].as<int>();
        omp_set_num_threads(num_threads); // Set the number of threads for parallel regions
        
        // enable nested parallelism
        bool nested = config_file["THREADING"]["nested"].as<bool>();
        if (nested == true) { omp_set_nested(1); }
    }
    else {
        omp_set_num_threads(1);
    }

    // logging prinouts
    this->verbose = config_file["INFO"]["verbose"].as<bool>();
}


// initialize the class variables
void Controller::initialize_variables()
{
    int n_leg = this->dynamics.n_leg;
    int N_u = this->params.N_u;

    // vectors and matrices that we use to sample from the distribution
    this->L.resize(3 * N_u * n_leg, 3 * N_u * n_leg);
    this->Z_vec_sample.resize(3 * N_u * n_leg);
    this->U_vec_sample.resize(3 * N_u * n_leg);
    this->U_traj_sample.resize(N_u);

    // vectors and matrices for the distribution updates
    this->mu.resize(3 * N_u * n_leg);
    this->Sigma.resize(3 * N_u * n_leg, 3 * N_u * n_leg);
    this->U_elite_matrix.resize(3 * N_u * n_leg, this->params.N_elite);
    this->U_elite_traj.resize(N_u);
    this->U_elite_vec.resize(3 * N_u * n_leg);
    this->eigval.resize(3 * N_u * n_leg);
    this->eigvec.resize(3 * N_u * n_leg, 3 * N_u * n_leg);
    this->eigvec_inv.resize(3 * N_u * n_leg, 3 * N_u * n_leg);
    this->I = Matrix_d::Identity(3 * N_u * n_leg, 3 * N_u * n_leg);

    // Monte Carlo results
    this->S_mc.resize(this->params.K);
    this->U_mc.resize(this->params.K);
    this->J_mc.resize(this->params.K);

    // Receding Horizon results
    this->S_elite.resize(this->params.N_elite);
    this->U_elite.resize(this->params.N_elite);
    this->J_elite.resize(this->params.N_elite);
}

void Controller::initialize_costs(YAML::Node config_file)
{
    // build the cost matrices from the diagonal elements
    Vector_1d_List Q_com_diags_temp = config_file["COST"]["Q_com_diags"].as<std::vector<double>>();
    Vector_1d_List Q_leg_diags_temp = config_file["COST"]["Q_leg_diags"].as<std::vector<double>>();
    Vector_1d_List Qf_com_diags_temp = config_file["COST"]["Qf_com_diags"].as<std::vector<double>>();
    Vector_1d_List Qf_leg_diags_temp = config_file["COST"]["Qf_leg_diags"].as<std::vector<double>>();
    Vector_1d_List R_diags_temp = config_file["COST"]["R_diags"].as<std::vector<double>>();
    Vector_1d_List R_rate_diags_temp = config_file["COST"]["R_rate_diags"].as<std::vector<double>>();

    Vector_6d Q_com_diags, Qf_com_diags;
    Vector_6d Q_leg_diag, Qf_leg_diag;
    for (int i = 0; i < Q_com_diags_temp.size(); i++) {
        Q_com_diags[i] = Q_com_diags_temp[i];
        Qf_com_diags[i] = Qf_com_diags_temp[i];
        Q_leg_diag[i] = Q_leg_diags_temp[i];
        Qf_leg_diag[i] = Qf_leg_diags_temp[i];
    }
    Vector_12d Q_leg_diags, Qf_leg_diags;
    Q_leg_diags << Q_leg_diag, Q_leg_diag;
    Qf_leg_diags << Qf_leg_diag, Qf_leg_diag;
    
    Vector_3d R_leg_diags, R_rate_leg_diags;
    for (int i = 0; i < R_diags_temp.size(); i++) {
        R_leg_diags[i] = R_diags_temp[i];
        R_rate_leg_diags[i] = R_rate_diags_temp[i];
    }
    Vector_6d R_diags, R_rate_diags;
    R_diags << R_leg_diags, R_leg_diags;
    R_rate_diags << R_rate_leg_diags, R_rate_leg_diags;

    // initialize the cost matrices
    this->params.Q_com  = Q_com_diags.asDiagonal();
    this->params.Q_leg  = Q_leg_diags.asDiagonal();
    this->params.Qf_com = Qf_com_diags.asDiagonal();
    this->params.Qf_leg = Qf_leg_diags.asDiagonal();
    this->params.R = R_diags.asDiagonal();
    this->params.R_rate = R_rate_diags.asDiagonal();

    // kinematic limits cost
    this->params.limits_enabled = config_file["COST"]["limits_enabled"].as<bool>();
    this->params.w_limits = config_file["COST"]["w_limits"].as<double>();

    // frcition cone cost
    this->params.friction_enabled = config_file["COST"]["friction_enabled"].as<bool>();
    this->params.w_friction = config_file["COST"]["w_friction"].as<double>();

    // gait cycel costs
    this->params.gait_enabled = config_file["COST"]["gait_enabled"].as<bool>();
    this->params.w_gait = config_file["COST"]["w_gait"].as<double>();
}


// construct the reference trajectory
void Controller::initialize_reference_trajectories(YAML::Node config_file)
{
    // set the desired COM trajectory
    Vector_1d_Traj time = config_file["REFERENCE"]["time"].as<std::vector<double>>();
    Vector_1d_Traj px_des_temp = config_file["REFERENCE"]["px_des"].as<std::vector<double>>();
    Vector_1d_Traj py_des_temp = config_file["REFERENCE"]["py_des"].as<std::vector<double>>();
    Vector_1d_Traj pz_des_temp = config_file["REFERENCE"]["pz_des"].as<std::vector<double>>();
    
    // ensure they are all the same size
    int N_ref = time.size();
    if (N_ref != px_des_temp.size() || N_ref != py_des_temp.size() || N_ref != pz_des_temp.size()) {
        throw std::runtime_error("px_des, py_des, pz_des, and time must be the same size");
    }

    // set the com reference trajectory
    Vector_3d_Traj p_com_ref(N_ref);
    for (int i = 0; i < N_ref; i++) {
        p_com_ref[i] << px_des_temp[i], py_des_temp[i], pz_des_temp[i];
    }

    // build the leg reference vector
    double r_des, theta_x_des, theta_y_des;
    r_des = config_file["REFERENCE"]["r_des"].as<double>();
    theta_x_des = config_file["REFERENCE"]["theta_x_des"].as<double>();
    theta_y_des = config_file["REFERENCE"]["theta_x_des"].as<double>();

    Vector_6d X_leg_ref_L, X_leg_ref_R;
    X_leg_ref_L << r_des,  theta_x_des, theta_y_des, 0.0, 0.0, 0.0;
    X_leg_ref_R << r_des, -theta_x_des, theta_y_des, 0.0, 0.0, 0.0;
    Vector_12d X_leg_ref;

    // set the gait cycle references
    double T_cycle_static, T_cycle_min, T_cycle_max, contact_static, contact_min, contact_max, phase_static, v_max;
    bool dynamic_cycle = config_file["REFERENCE"]["dynamic_cycle"].as<bool>();
    T_cycle_static = config_file["REFERENCE"]["T_cycle_static"].as<double>();
    T_cycle_min = config_file["REFERENCE"]["T_cycle_min"].as<double>();
    T_cycle_max = config_file["REFERENCE"]["T_cycle_max"].as<double>();
    contact_static = config_file["REFERENCE"]["contact_static"].as<double>();
    contact_min = config_file["REFERENCE"]["contact_min"].as<double>();
    contact_max = config_file["REFERENCE"]["contact_max"].as<double>();
    phase_static = config_file["REFERENCE"]["phase_static"].as<double>();
    v_max = config_file["REFERENCE"]["v_max"].as<double>();

    // set the leg reference trajectory
    this->ref_sim.t_ref = time;
    this->ref_sim.p_com_ref = p_com_ref;
    this->ref_sim.X_leg_ref << X_leg_ref_L, X_leg_ref_R;
    this->ref_sim.N_ref = N_ref;
    this->ref_sim.dynamic_cycle = dynamic_cycle;
    this->ref_sim.T_cycle_static = T_cycle_static;
    this->ref_sim.T_cycle_min = T_cycle_min;
    this->ref_sim.T_cycle_max = T_cycle_max;
    this->ref_sim.contact_static = contact_static;
    this->ref_sim.contact_min = contact_min;
    this->ref_sim.contact_max = contact_max;
    this->ref_sim.phase_static = phase_static;
    this->ref_sim.v_max = v_max;

    // set the horizon reference
    this->ref_horizon.X_com_ref.resize(this->params.N_x);
    this->ref_horizon.X_leg_ref.resize(this->params.N_x);
    this->ref_horizon.domain_ref.resize(this->params.N_x);
}


// construct the intial distribution
void Controller::initialize_distribution(YAML::Node config_file)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg; // number of legs
    int Nu = this->params.N_u;        // number of control knots

    // set the epsilon for numerical stability of covariance matrix
    this->dist.epsilon = config_file["DIST_PARAMS"]["epsilon"].as<double>();

    // compute the theoretical lower bound (lower bound of Frobenius norm, assuming epsilon eigs)
    this->min_cov_norm = std::sqrt(3 * n_leg * Nu ) * this->dist.epsilon;

    // initialize the matrices
    this->dist.mean.resize(3 * n_leg * Nu);
    this->dist.cov.resize(3 * n_leg * Nu, 3 * n_leg * Nu);
    this->dist.mean.setZero();
    this->dist.cov.setZero();

    // set the initial mean
    std::vector<double> mean_temp = config_file["DIST_PARAMS"]["mu"].as<std::vector<double>>();
    Vector_6d mean;
    mean << mean_temp[0], mean_temp[1], mean_temp[2], 
            mean_temp[0], mean_temp[1], mean_temp[2];
    for (int i = 0; i < Nu; i++) {
        this->dist.mean.segment<6>(i * 3 * n_leg) = mean;
    }

    // set the initial covariance
    std::vector<double> cov_temp = config_file["DIST_PARAMS"]["sigma"].as<std::vector<double>>();
    Vector_6d cov_diags;
    cov_diags << cov_temp[0], cov_temp[1], cov_temp[2], cov_temp[0], cov_temp[1], cov_temp[2];
    Matrix_6d cov = cov_diags.asDiagonal();
    for (int i = 0; i < Nu; i++) {
        this->dist.cov.block<6, 6>(3 * i * n_leg, 3 * i * n_leg) = cov;
    }

    // set if covariance should be strictly diagonal
    this->dist.diag_cov = config_file["DIST_PARAMS"]["diag_cov"].as<bool>();

    // set the random seed
    this->dist.seed_enabled = config_file["DIST_PARAMS"]["seed_enabled"].as<bool>();
    this->dist.seed = config_file["DIST_PARAMS"]["seed"].as<int>();

    // set if the control inputs should be saturated
    this->dist.saturate = config_file["DIST_PARAMS"]["saturate"].as<bool>();

    // create random device
    std::random_device rand_device;

    // use the random device to seed Mersenne Twister generator
    std::mt19937 rand_generator(rand_device());

    // set the seed if enabled
    if (this->dist.seed_enabled) {
        rand_generator.seed(this->dist.seed);
    }

    // Create a normal distribution
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    
    // set the random number generator and normal distribution
    this->rand_generator = rand_generator;
    this->normal_dist = normal_dist;
}


// sample input trajectories from the distribution
void Controller::sample_input_trajectory()
{
    // perform cholesky decomposition 
    Eigen::LLT<Matrix_d> llt(this->dist.cov);  
    this->L = llt.matrixL();

    // check if the covariance is positive definite
    if (llt.info() == Eigen::NumericalIssue) {
        throw std::runtime_error("Covariance matrix is possibly not positive definite");
    }

    // for saturating hte inputs
    double rdot_lim = this->dynamics.params.rdot_lim;
    double thetadot_lim = this->dynamics.params.thetadot_lim;

    // U ~ N(mu, Sigma) <=> U = L * Z + mu; Z ~ N(0, I)
    // Sample the input trajectories
    Vector_6d u;                                    
    for (int i = 0; i < this->params.K; i++) {
        
        // populate the Z vector; Z ~ N(0, I)
        for (int j = 0; j < this->dist.mean.size(); j++) {
            this->Z_vec_sample(j) = this->normal_dist(this->rand_generator);
        }

        // generate the vectorized input trajectory
        this->U_vec_sample = L * Z_vec_sample + this->dist.mean;

        // unvectorize U_vec into U_traj
        for (int k = 0; k < this->params.N_u; k++) {
            
            // grab the k'th control input
            u = this->U_vec_sample.segment<6>(6 * k);

            // limit the input to the bounds
            if (this->dist.saturate == true) {
                u(0) = std::max(-rdot_lim, std::min(rdot_lim, u(0)));
                u(1) = std::max(-thetadot_lim, std::min(thetadot_lim, u(1)));
                u(2) = std::max(-thetadot_lim, std::min(thetadot_lim, u(2)));
                u(3) = std::max(-rdot_lim, std::min(rdot_lim, u(3)));
                u(4) = std::max(-thetadot_lim, std::min(thetadot_lim, u(4)));
                u(5) = std::max(-thetadot_lim, std::min(thetadot_lim, u(5)));
            }

            // insert into the trajectory
            this->U_traj_sample[k] = u;
        }

        // store the input trajectory
        this->U_mc[i] = U_traj_sample;
    }
}


// compute mean and covariance from a bundle of control inputs
void Controller::update_distribution_params(const Vector_6d_Traj_Bundle& U_elite)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg;
    int Nu = this->params.N_u;

     // initialize the meanto zero
    this->mu.setZero();

    // compute the mean
    Vector_6d u;    
    for (int i = 0; i < this->params.N_elite; i++) {

        // grab the i'th elite trajectory
        this->U_elite_traj = U_elite[i];

        // vectorize the i'th elite trajectory
        for (int j = 0; j < Nu; j++) {
            u = this->U_elite_traj[j];
            this->U_elite_vec.segment<6>(6 * j) = u;
        }   

        // update the mean
        this->mu += this->U_elite_vec;

        // insert the i'th elite trajectory into the data matrix
        this->U_elite_matrix.col(i) = this->U_elite_vec;
    }
    // compute the mean
    this->mu /= this->params.N_elite;

    // compute the sample covariance (K-1 b/c Bessel correction)
    this->Sigma = (1.0 / (this->params.N_elite-1)) * (this->U_elite_matrix.colwise() - this->mu) * (this->U_elite_matrix.colwise() - this->mu).transpose();
    
    // compute the eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<Matrix_d> eig(this->Sigma);
    if (eig.info() == Eigen::NumericalIssue) {
        throw std::runtime_error("Covariance matrix is possibly not positive definite");
    }
    this->eigvec = eig.eigenvectors();
    this->eigval = eig.eigenvalues();
    this->eigvec_inv = eigvec.inverse();

    // modify eigenvalues with epsilon if it gets too low
    for (int i = 0; i < eigval.size(); i++) {
        this->eigval(i) = std::max(this->eigval(i), this->dist.epsilon);
    }

    // rebuild the covariance matrix with the eigenvalue decomposition, add epsilon to eigenvalues
    this->Sigma = this->eigvec * this->eigval.asDiagonal() * this->eigvec_inv;

    //  if stricly diagonal covariance option
    if (this->dist.diag_cov == true) {
        // set the covariance to be diagonal, (Hadamard Product is coeff-wise product, cov * I = diag(cov))
        this->Sigma = this->Sigma.cwiseProduct(this->I);
    }

    // update the distribution
    this->dist.mean = mu;
    this->dist.cov = Sigma;    
}


// generate a reference trajectory for the predictive control to track
void Controller::generate_reference_trajectory(double t_sim, const Vector_6d& x0_com)
{
    // grab the center of mass velocity
    Vector_3d v_com = x0_com.segment<3>(3);

    // build the COM reference trajectory
    Vector_6d xi_com_ref;
    Vector_3d pi_com_ref, p0_com_ref, pf_com_ref;
    Vector_3d vi_com_ref;
    Vector_2i d_ref = Vector_2i::Zero();
    double t, t0, tf;
    for (int i = 0; i < this->params.N_x; i++) {
    
        // compute the time given the global reference
        t = i * this->params.dt_x + t_sim;

        // find where t is in the time reference
        auto it = std::upper_bound(this->ref_sim.t_ref.begin(), this->ref_sim.t_ref.end(), t);
        int idx = std::distance(this->ref_sim.t_ref.begin(), it) - 1; // return -1 if before the first element
                                                                      // return N_ref - 1 if after the last element

        // beyond last element
        if (idx == this->ref_sim.N_ref - 1) {
            // set to last point with zero velocity
            pi_com_ref = this->ref_sim.p_com_ref[this->ref_sim.N_ref - 1];
            vi_com_ref << 0.0, 0.0, 0.0;
        }
        else{
            // linear interpolation
            t0 = this->ref_sim.t_ref[idx];
            tf = this->ref_sim.t_ref[idx + 1];
            p0_com_ref = this->ref_sim.p_com_ref[idx];
            pf_com_ref = this->ref_sim.p_com_ref[idx + 1];
            vi_com_ref = (pf_com_ref - p0_com_ref) / (tf - t0);
            pi_com_ref = p0_com_ref + vi_com_ref * (t - t0);
        }

        // build the COM reference
        xi_com_ref << pi_com_ref, vi_com_ref;

        // TODO: finish this contact reference
        // build the domain reference
        if (this->params.gait_enabled == true) {

            // dynamic gait cycle
            if (this->ref_sim.dynamic_cycle == true) {
                // TODO: implement dynamic gait cycle
                // use v_com
            }

            // static gait cycle
            else {
                // left leg contact reference
                if (std::fmod(t / this->ref_sim.T_cycle_static, 1.0) <= this->ref_sim.contact_static) {
                    d_ref(0) = 1;
                }
                else {
                    d_ref(0) = 0;
                }

                // right leg contact reference (assuming a 50% offset)
                if (std::fmod(t / this->ref_sim.T_cycle_static - this->ref_sim.phase_static, 1.0) <= this->ref_sim.contact_static) {
                    d_ref(1) = 1;
                }
                else {
                    d_ref(1) = 0;
                }
            }
        }

        // insert values into horizon reference trajectory
        this->ref_horizon.X_com_ref[i] = xi_com_ref;               // X_com_ref
        this->ref_horizon.X_leg_ref[i] = this->ref_sim.X_leg_ref;  // X_com_ref
        this->ref_horizon.domain_ref[i] = d_ref;                   // domain_ref
    }
}


// evaluate log barrier function cost on legs
double Controller::cost_limits(const Vector_12d& x_leg) {   
    
    // want to compute the cost of the log barrier function
    double J_limits = 0.0;

    // leg state constraints
    double r_min = this->dynamics.params.r_min;
    double r_max = this->dynamics.params.r_max;
    double theta_x_min = this->dynamics.params.theta_x_min;
    double theta_x_max = this->dynamics.params.theta_x_max;
    double theta_y_min = this->dynamics.params.theta_y_min;
    double theta_y_max = this->dynamics.params.theta_y_max;
    double rdot_lim = this->dynamics.params.rdot_lim;
    double thetadot_lim = this->dynamics.params.thetadot_lim;
    double w_limits = this->params.w_limits;

    // compute the costs
    Vector_6d xi_leg;
    double r, theta_x, theta_y, rdot, thetadot_x, thetadot_y;
    for (int i = 0; i < this->dynamics.n_leg; i++) {
        
        // unpack the leg state
        xi_leg = x_leg.segment<6>(6 * i);
        r = xi_leg(0);
        theta_x = xi_leg(1);
        theta_y = xi_leg(2);
        rdot = xi_leg(3);
        thetadot_x = xi_leg(4);
        thetadot_x = xi_leg(5);

        // compute kinematic limits cost
        J_limits += (r < r_min) ? w_limits * (r - r_min) * (r - r_min) :
                    (r > r_max) ? w_limits * (r - r_max) * (r - r_max) : 0.0;
        if (i == 0) {       // Left leg
            J_limits += (theta_x < theta_x_min) ? w_limits * (theta_x - theta_x_min) * (theta_x - theta_x_min) :
                        (theta_x > theta_x_max) ? w_limits * (theta_x - theta_x_max) * (theta_x - theta_x_max) : 0.0;
        }
        else if (i == 1) {  // Right leg
            J_limits += (theta_x < -theta_x_max) ? w_limits * (theta_x + theta_x_max) * (theta_x + theta_x_max) :
                        (theta_x > -theta_x_min) ? w_limits * (theta_x + theta_x_min) * (theta_x + theta_x_min) : 0.0;
        }
        J_limits += (theta_y < theta_y_min) ? w_limits * (theta_y - theta_y_min) * (theta_y - theta_y_min) :
                    (theta_y > theta_y_max) ? w_limits * (theta_y - theta_y_max) * (theta_y - theta_y_max) : 0.0;

        // compute the velocity limits cost
        J_limits += (rdot < -rdot_lim) ? w_limits * (rdot + rdot_lim) * (rdot + rdot_lim) :
                    (rdot >  rdot_lim) ? w_limits * (rdot - rdot_lim) * (rdot - rdot_lim) : 0.0;
        J_limits += (thetadot_x < -thetadot_lim) ? w_limits * (thetadot_x + thetadot_lim) * (thetadot_x + thetadot_lim) :
                    (thetadot_x >  thetadot_lim) ? w_limits * (thetadot_x - thetadot_lim) * (thetadot_x - thetadot_lim) : 0.0;
        J_limits += (thetadot_y < -thetadot_lim) ? w_limits * (thetadot_y + thetadot_lim) * (thetadot_y + thetadot_lim) :
                    (thetadot_y >  thetadot_lim) ? w_limits * (thetadot_y - thetadot_lim) * (thetadot_y - thetadot_lim) : 0.0;
    }

    return J_limits;
}


// evaluate the cost function given a solution
double Controller::cost_function(const ReferenceLocal& ref, const Solution& Sol, const Vector_6d_Traj& U)
{
    // number of knot points 
    int Nx = this->params.N_x; // number of state knots
    int Nu = this->params.N_u; // number of control knots

    // ************************************ COM COST ************************************

    // compute the COM cost
    double J_com = 0.0;
    Vector_6d ei_com;
    // integrated cost
    for (int i = 0; i < Nx-1; i++) {
        ei_com = Sol.x_sys_t[i].head<6>() - ref.X_com_ref[i];
        J_com += ei_com.transpose() * this->params.Q_com * ei_com;
    }
    //terminal cost
    ei_com = Sol.x_sys_t[Nx-1].head<6>() - ref.X_com_ref[Nx-1];
    J_com += ei_com.transpose() * this->params.Qf_com * ei_com;

    // ************************************ LEG COST ************************************

    // compute the LEG cost
    double J_legs = 0.0;
    Vector_12d ei_leg;
    for (int i = 0; i < Nx-1; i++) {
        ei_leg = Sol.x_leg_t[i] - ref.X_leg_ref[i];
        J_legs += ei_leg.transpose() * this->params.Q_leg * ei_leg;
    }
    // terminal cost
    ei_leg = Sol.x_leg_t[Nx-1] - ref.X_leg_ref[Nx-1];
    J_legs += ei_leg.transpose() * this->params.Qf_leg * ei_leg;

    // ************************************ KINEMATIC LIMITS COST ************************************
    
    // compute the kinematic limit vioalation cost
    double J_limits = 0.0;
    if (this->params.limits_enabled) {
        for (int i = 0; i < Nx; i++) {
            J_limits += this->cost_limits(Sol.x_leg_t[i]);
        }
    }

    // ************************************ INPUT COST ************************************

    // penalize the velocity input (smaller velocities)
    double J_input = 0.0;
    Vector_6d ui;
    for (int i = 0; i < Nu; i++) {        
        // left and right leg input vector
        ui << U[i];

        // compute quadratic cost
        J_input += ui.transpose() * this->params.R * ui;
    }

    // penalize the rate of change of the input (smaller accelerations)
    Vector_6d u_delta;
    for (int i = 0; i < Nu - 1; i++) {
        // left and right leg input vector
        u_delta = (U[i + 1] - U[i]) / this->params.dt_u;

        // compute quadratic cost
        J_input += (u_delta).transpose() * this->params.R_rate * (u_delta);
    }

    // ************************************ FRICTION COST ************************************

    // compute the friction cone cost
    // lambda_t in {lambda in R^3 | lambda_z >= 0, ||lambda_xy|| <= mu * lambda_z}
    double J_friction = 0.0;
    if (this->params.friction_enabled) {

        // useful variables
        Domain d;
        Vector_6d lambdas;
        Vector_3d lam;
        double lam_z, lam_xy_norm;
        double coeff = this->dynamics.params.friction_coeff;
        double w_friction = this->params.w_friction;
        
        // loop through the trajectory
        for (int i = 0; i < Nx; i++) {
            
            // grab the current domain at time t
            d = Sol.domain_t[i];

            // if a leg is in stance, compute the friction cost
            if (d[0] == Contact::STANCE || d[1] == Contact::STANCE) {    
                
                // get the lambda vector at time t
                lambdas = Sol.lambda_t[i];

                // for each leg compute the friction cone violation cost
                for (int j = 0; j < this->dynamics.n_leg; j++) {

                    if (d[j] == Contact::STANCE) {
                        // get the i'th leg lambda vector
                        lam = lambdas.segment<3>(3 * j);
                        lam_z = lam(2);

                        // check that lambda_z >= 0 (negative is possible due to coarse switching surface approximation)
                        if (lam_z > 0.0) {
                            // check that ||lambda_xy|| <= mu * lambda_z
                            lam_xy_norm = lam.head<2>().norm();
                            J_friction += (lam_xy_norm > coeff * lam_z) ? w_friction * (lam_xy_norm - coeff * lam_z) * (lam_xy_norm - coeff * lam_z) : 0.0;                            
                        }
                    }
                }
            }
        }
    }

    // ************************************ GAIT COST ************************************

    // compute the gait cost
    double J_gait = 0.0;
    if (this->params.gait_enabled) {
        
        // useful variables
        Vector_2i d_ref, d;
        for (int i = 0; i < Nx; i++) {
            
            // grab the current domain at time t
            d_ref = ref.domain_ref[i];
            
            // convert to binary
            d(0) = Sol.domain_t[i][0] == Contact::STANCE ? 1 : 0;
            d(1) = Sol.domain_t[i][1] == Contact::STANCE ? 1 : 0;

            if (d != d_ref) {
                J_gait += this->params.w_gait;
            }
        }
    }

    // ************************************ TOTAL COST ************************************

    return J_com + J_legs + J_limits + J_input + J_friction + J_gait;
}


// perform open loop rollouts
void Controller::monte_carlo(double t_sim, Vector_12d x0_sys, Vector_6d p0_feet, Domain d0)
{
    // generate new bundle of input trajectories
    auto t0 = std::chrono::high_resolution_clock::now();
    this-> sample_input_trajectory(); 
    auto tf = std::chrono::high_resolution_clock::now();
    if (this->verbose == true) {
        std::cout << "Time to sample input trajectories: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
    }

    // generate the reference trajectory
    t0 = std::chrono::high_resolution_clock::now();
    this->generate_reference_trajectory(t_sim, x0_sys.head<6>());
    tf = std::chrono::high_resolution_clock::now();
    if (this->verbose == true) {
        std::cout << "Time to generate reference trajectory: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
    }

    // loop over the input trajectories
    Solution sol;
    this->dynamics.resizeSolution(sol, this->params.T_x);
    double cost = 0.0;
    t0 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for private(sol, cost) 
    for (int k = 0; k < this->params.K; k++) {
        
        // initialize the cost
        cost = 0.0;

        // perform the rollout
        auto t1 = std::chrono::high_resolution_clock::now();
        sol = this->dynamics.RK3_rollout(this->params.T_x, this->params.T_u, x0_sys, p0_feet, d0, this->U_mc[k]);
        auto t2 = std::chrono::high_resolution_clock::now();
        if (this->verbose == true) {
            // std::cout << "Time to perform rollout: " << std::chrono::duration<double, std::milli>(t2 - t1).count() << " ms" << std::endl;
        }

        // compute the cost
        t1 = std::chrono::high_resolution_clock::now();
        cost = this->cost_function(this->ref_horizon, sol, this->U_mc[k]);
        t2 = std::chrono::high_resolution_clock::now();
        if (this->verbose == true) {
            // std::cout << "Time to compute cost: " << std::chrono::duration<double, std::milli>(t2 - t1).count() << " ms" << std::endl;
        }
    
        // store the results (use critical sections to avoid race conditions if necessary)
        t1 = std::chrono::high_resolution_clock::now();
        #pragma omp critical
        {
            this->S_mc[k] = sol;
            this->J_mc[k] = cost;
        }
        t2 = std::chrono::high_resolution_clock::now();
        if (this->verbose == true) {
            // std::cout << "Time to store results: " << std::chrono::duration<double, std::milli>(t2 - t1).count() << " ms" << std::endl;
        }
    }
    tf = std::chrono::high_resolution_clock::now();
    if (this->verbose == true) {
        std::cout << "Time to perform rollouts: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
    }

}


// select solutions based on cost
void Controller::sort_trajectories()
{
    // Create an index vector
    Vector_1i_List idx(this->params.K);
    std::iota(idx.begin(), idx.end(), 0);

    // Use nth_element to bring the top N_elite elements to the front in O(n) time
    std::nth_element(idx.begin(), idx.begin() + this->params.N_elite, idx.end(),
                    [this](int i1, int i2) { return this->J_mc[i1] < this->J_mc[i2]; });

    // Populate elite solutions
    for (int i = 0; i < this->params.N_elite; i++) {
        this->S_elite[i] = this->S_mc[idx[i]];
        this->U_elite[i] = this->U_mc[idx[i]];
        this->J_elite[i] = this->J_mc[idx[i]];
    }
}


// perform sampling predictive control of your choice here
RHC_Result Controller::sampling_predictive_control(double t_sim, Vector_12d x0_sys, Vector_6d p0_foot, Domain d0)
{
    // Monte Carlo Result to return
    MC_Result mc;

    // perform the CEM iterations
    auto t0_total = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < this->params.CEM_iters; i++) {

        // perform monte carlo simulation
        auto t0 = std::chrono::high_resolution_clock::now();
        this->monte_carlo(t_sim, x0_sys, p0_foot, d0);
        auto tf = std::chrono::high_resolution_clock::now();
        if (this->verbose == true) {
            std::cout << "Time to perform Monte Carlo: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
        }

        // sort the cost vector in ascending order
        t0 = std::chrono::high_resolution_clock::now();
        this->sort_trajectories();
        tf = std::chrono::high_resolution_clock::now();
        if (this->verbose == true) {
            std::cout << "Time to sort trajectories: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
        }

        // update the distribution parameters
        t0 = std::chrono::high_resolution_clock::now();
        this->update_distribution_params(U_elite);
        tf = std::chrono::high_resolution_clock::now();
        if (this->verbose == true) {
            std::cout << "Time to update distribution parameters: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
        }

        // print some info
        if (this->verbose == true) 
        {
            std::cout << "-----------------------------------" << std::endl;
            std::cout << "CEM Iteration: " << i+1 << std::endl;
            std::cout << "-----------------------------------" << std::endl;
            std::cout << "Time for iteration: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
            std::cout << "Smallest cost: " << J_elite[0] << std::endl;   
            std::cout << "Norm of covariance: " << this->dist.cov.norm() << ", Theoretical min: " << this->min_cov_norm << std::endl;
            std::cout << std::endl;
        }

    }
    auto tf_total = std::chrono::high_resolution_clock::now();
    double T_tot = std::chrono::duration<double>(tf_total - t0_total).count();

    if (this->verbose == true) 
    {
        std::cout << "CEM complete" << std::endl;
        std::cout << "Total time: " << T_tot << " sec" << std::endl;
        std::cout << "Average time per iteration: " << T_tot / this->params.CEM_iters << " [sec], " << this->params.CEM_iters/T_tot << " [Hz]" << std::endl;
        std::cout << "Average Rollout time: " << T_tot / (this->params.CEM_iters * this->params.K) * 1000000.0 << " [us]" << std::endl << std::endl;
    }

    // do a final sort to get the best solution
    Vector_1i_List idx(this->params.N_elite);
    std::iota(idx.begin(), idx.end(), 0);
    
    // get only the best solution
    std::nth_element(idx.begin(), idx.begin() + 1, idx.end(),
                    [this](int i1, int i2) { return this->J_elite[i1] < this->J_elite[i2]; });

    // return the best solution
    RHC_Result rhc;
    rhc.S = this->S_elite[idx[0]];
    rhc.U = this->U_elite[idx[0]];

    return rhc;
}
