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


void Controller::initialize_variables()
{
    int n_leg = this->dynamics.n_leg;
    int N_u = this->params.N_u;

    // vectors and matrices that we use to sample from the distribution
    this->L.resize(3 * N_u * n_leg, 3 * N_u * n_leg);
    this->Z_vec_sample.resize(3 * N_u * n_leg);
    this->U_vec_sample.resize(3 * N_u * n_leg);

    // trajectory bundles and samples
    this->U_traj_samples.resize(this->params.K);
    this->U_traj_sample.resize(N_u);
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

    // log barrier function
    this->params.log_barrier_enabled = config_file["COST"]["log_barrier_enabled"].as<bool>();
    this->params.J_barrier = config_file["COST"]["J_barrier"].as<double>();
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

    // set the leg reference trajectory
    ReferenceGlobal ref;
    ref.t_ref = time;
    ref.p_com_ref = p_com_ref;
    ref.X_leg_ref << X_leg_ref_L, X_leg_ref_R;
    ref.N_ref = N_ref;

    // set the reference
    this->ref = ref;
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

    // set the random 
    this->dist.seed_enabled = config_file["DIST_PARAMS"]["seed_enabled"].as<bool>();
    this->dist.seed = config_file["DIST_PARAMS"]["seed"].as<int>();

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
Vector_6d_Traj_Bundle Controller::sample_input_trajectory()
{
    // perform cholesky decomposition 
    Eigen::LLT<Matrix_d> llt(this->dist.cov);  
    this->L = llt.matrixL();

    // check if the covariance is positive definite
    if (llt.info() == Eigen::NumericalIssue) {
        throw std::runtime_error("Covariance matrix is possibly not positive definite");
    }

    // U ~ N(mu, Sigma) <=> U = L * Z + mu; Z ~ N(0, I)
    // initialize the input trajectory bundle
    Vector_6d u;                                    

    // loop over the number of samples
    for (int i = 0; i < this->params.K; i++) {
        
        // populate the Z vector; Z ~ N(0, I)
        for (int j = 0; j < this->dist.mean.size(); j++) {
            this->Z_vec_sample(j) = this->normal_dist(this->rand_generator);
        }

        // generate the vectorized input trajectory
        this->U_vec_sample = L * Z_vec_sample + this->dist.mean;

        // unvectorize U_vec into U_traj
        for (int k = 0; k < this->params.N_u; k++) {
            u = this->U_vec_sample.segment<6>(6 * k);
            this->U_traj_sample[k] = u;
        }

        // store the input trajectory
        U_traj_samples[i] = U_traj_sample;
    }

    return U_traj_samples; // TODO: still need to return nothing here
}


// TODO: there is code optimziation work to be done here
// compute mean and covariance from a bundle of control inputs
void Controller::update_distribution_params(const Vector_6d_Traj_Bundle& U_eilte)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg;
    int Nu = this->params.N_u;

    // initialize the mean and covariance
    Vector_d mean; // TODO: have this as a class variable
    Matrix_d cov;  // TODO: have this as a class variable
    mean.resize(3 * Nu * n_leg);
    cov.resize(3 * Nu * n_leg, 3 * Nu * n_leg);

    // used for computing the mean
    Matrix_d U_data; // TODO: have this as a class variable
    U_data.resize(3 * Nu * n_leg, this->params.N_elite);

    // compute the mean
    Vector_6d_Traj U_traj(Nu); // TODO: have this as a class variable
    Vector_6d u;
    Vector_d U_vec;
    U_vec.resize(3 * Nu * n_leg);
    for (int i = 0; i < this->params.N_elite; i++) {

        // vectorize the input trajectory
        U_traj = U_eilte[i];
        for (int j = 0; j < Nu; j++) {
            u = U_traj[j];
            U_vec.segment<6>(6 * j) = u;
        }   
        mean += U_vec;

        // insert into data matrix to use later
        U_data.col(i) = U_vec;
    }
    mean /= this->params.N_elite;

    // compute the sample covariance (K-1 b/c Bessel correction)
    cov = (1.0 / (this->params.N_elite-1)) * (U_data.colwise() - mean) * (U_data.colwise() - mean).transpose();
    
    // compute the eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<Matrix_d> eig(cov);
    if (eig.info() == Eigen::NumericalIssue) {
        throw std::runtime_error("Covariance matrix is possibly not positive definite");
    }
    Matrix_d eigvec = eig.eigenvectors();
    Vector_d eigval = eig.eigenvalues();
    Matrix_d eigvec_inv = eigvec.inverse();

    // std::cout << "c" << std::endl;

    // modify eigenvalues with epsilon if it gets too low
    for (int i = 0; i < eigval.size(); i++) {
        eigval(i) = std::max(eigval(i), this->dist.epsilon);
    }

    // rebuild the covariance matrix with the eigenvalue decomposition, add epsilon to eigenvalues
    cov = eigvec * eigval.asDiagonal() * eigvec_inv;

    //  if stricly diagonal covariance option
    if (this->dist.diag_cov == true) {
        // set the covariance to be diagonal, (Hadamard Product, cov * I = diag(cov))
        cov = cov.cwiseProduct(Matrix_d::Identity(3 * Nu * n_leg, 3 * Nu * n_leg));
    }

    // std::cout << "d" << std::endl;

    // update the distribution
    this->dist.mean = mean;
    this->dist.cov = cov;    
}


// TODO: there is code optimziation work to be done here
// generate a reference trajectory for the predictive control to track
ReferenceLocal Controller::generate_reference_trajectory(double t_sim, const Vector_6d& x0_com)
{
    // pass reference to the dynamics
    Vector_6d_Traj X_com_ref(this->params.N_x);  // TODO: have this as a class variable
    Vector_12d_Traj X_leg_ref(this->params.N_x); // TODO: have this as a class variable

    // build the COM reference trajectory
    Vector_6d xi_com_ref;
    Vector_3d pi_com_ref, p0_com_ref, pf_com_ref;
    Vector_3d vi_com_ref;
    double t, t0, tf;
    for (int i = 0; i < this->params.N_x; i++) {
    
        // compute the time given the global reference
        t = i * this->params.dt_x + t_sim;

        // find where t is in the time reference
        auto it = std::upper_bound(this->ref.t_ref.begin(), this->ref.t_ref.end(), t);
        int idx = std::distance(this->ref.t_ref.begin(), it) - 1; // return -1 if before the first element
                                                                  // return N_ref - 1 if after the last element

        // beyond last element
        if (idx == this->ref.N_ref - 1) {
            // set to last point with zero velocity
            pi_com_ref = this->ref.p_com_ref[this->ref.N_ref - 1];
            vi_com_ref << 0.0, 0.0, 0.0;
        }
        else{
            // linear interpolation
            t0 = this->ref.t_ref[idx];
            tf = this->ref.t_ref[idx + 1];
            p0_com_ref = this->ref.p_com_ref[idx];
            pf_com_ref = this->ref.p_com_ref[idx + 1];
            vi_com_ref = (pf_com_ref - p0_com_ref) / (tf - t0);
            pi_com_ref = p0_com_ref + vi_com_ref * (t - t0);
        }

        // build the reference
        xi_com_ref << pi_com_ref, vi_com_ref;

        // insert into trajectory
        X_com_ref[i] = xi_com_ref;
    }

    // build the LEG reference trajectory
    for (int i = 0; i < this->params.N_x; i++) {
        // insert into trajectory
        X_leg_ref[i] = this->ref.X_leg_ref;
    }

    // insert the references
    ReferenceLocal ref_horizon; // TODO: have this as a class variable
    ref_horizon.X_com_ref = X_com_ref;
    ref_horizon.X_leg_ref = X_leg_ref;

    return ref_horizon;
}


// evaluate log barrier function cost on legs
double Controller::cost_log_barrier(const Vector_12d& x_leg) {   
    
    // want to compute the cost of the log barrier function
    double J_log = 0.0;

    // leg state constraints
    double r_min = this->dynamics.params.r_min;
    double r_max = this->dynamics.params.r_max;
    double theta_x_min = this->dynamics.params.theta_x_min;
    double theta_x_max = this->dynamics.params.theta_x_max;
    double theta_y_min = this->dynamics.params.theta_y_min;
    double theta_y_max = this->dynamics.params.theta_y_max;
    double rdot_lim = this->dynamics.params.rdot_lim;
    double thetadot_lim = this->dynamics.params.thetadot_lim;

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

        // compare the state to the limits (this alternative is a simple heuristic for barrier-like behavior)
        J_log += (r < r_min || r > r_max) ? this->params.J_barrier : 0;
        J_log += (theta_x < theta_x_min || theta_x > theta_x_max) ? this->params.J_barrier : 0;
        if (i == 0) {
            J_log += (theta_y < theta_y_min || theta_y > theta_y_max) ? this->params.J_barrier : 0;
        }
        else if (i == 1) {
            J_log += (theta_y > -theta_y_min || theta_y < -theta_y_max) ? this->params.J_barrier : 0;
        }
        J_log += (std::abs(rdot) > rdot_lim) ? this->params.J_barrier : 0;
        J_log += (std::abs(thetadot_x) > thetadot_lim) ? this->params.J_barrier : 0;
        J_log += (std::abs(thetadot_y) > thetadot_lim) ? this->params.J_barrier : 0;
    }

    return J_log;
}


// TODO: there is code optimziation work to be done here. MAJOR source of inefficiency
// evaluate the cost function given a solution
double Controller::cost_function(const ReferenceLocal& ref, const Solution& Sol, const Vector_6d_Traj& U)
{
    // number of knot points 
    int Nx = this->params.N_x; // number of state knots
    int Nu = this->params.N_u; // number of control knots

    // ************************************ COM COST ************************************

    // build the COM trajectory
    Vector_6d_Traj X_com(Nx);
    for (int i = 0; i < Nx; i++) {
        X_com[i] = Sol.x_sys_t[i].head<6>();
    }

    // compute the COM cost
    Vector_6d ei_com;
    double J_com = 0.0;
    // integrated cost
    for (int i = 0; i < Nx-1; i++) {
        ei_com = X_com[i] - ref.X_com_ref[i];
        J_com += ei_com.transpose() * this->params.Q_com * ei_com;
    }
    //terminal cost
    ei_com = X_com[Nx-1] - ref.X_com_ref[Nx-1];
    J_com += ei_com.transpose() * this->params.Qf_com * ei_com;

    // ************************************ LEG COST ************************************

    // compute the LEG cost
    Vector_12d ei_leg;
    double J_legs = 0.0;
    // integrated cost
    for (int i = 0; i < Nx-1; i++) {
        ei_leg = Sol.x_leg_t[i] - ref.X_leg_ref[i];
        J_legs += ei_leg.transpose() * this->params.Q_leg * ei_leg;
    }
    // terminal cost
    ei_leg = Sol.x_leg_t[Nx-1] - ref.X_leg_ref[Nx-1];
    J_legs += ei_leg.transpose() * this->params.Qf_leg * ei_leg;

    // log barrier function cost
    double J_legs_barrier = 0.0;
    
    // compute the log barrier function cost
    if (this->params.log_barrier_enabled) {
        for (int i = 0; i < Nx; i++) {
            J_legs_barrier += this->cost_log_barrier(Sol.x_leg_t[i]);
        }
        J_legs += J_legs_barrier;
    }
    // if not enabled, set to zero
    else {
        J_legs_barrier = 0.0;
        J_legs += J_legs_barrier;
    }

    // ************************************ INPUT COST ************************************

    // penalize the velocity input (smaller velocities)
    Vector_6d ui;
    double J_input = 0.0;
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

    // ************************************ TOTAL COST ************************************

    // total cost
    double J_total; 
    J_total = J_com + J_legs + J_input;

    return J_total;
}


// perform open loop rollouts
MC_Result Controller::monte_carlo(double t_sim, Vector_12d x0_sys, Vector_6d p0_feet, Domain d0)
{
    // generate bundle of input trajectories
    Vector_6d_Traj_Bundle U_bundle(this->params.K);            // TODO: have this as a class variable
    U_bundle = this->sample_input_trajectory();  

    // initialize the containers for the solutions
    Solution_Bundle Sol_bundle(this->params.K); // TODO: have this as a class variable
    Vector_1d_List J(this->params.K);

    // generate the reference trajectory
    ReferenceLocal ref;
    ref = this->generate_reference_trajectory(t_sim, x0_sys.head<6>());

    // loop over the input trajectories
    Solution sol;
    this->dynamics.resizeSolution(sol, this->params.T_x);
    double cost = 0.0;
    #pragma omp parallel for private(sol, cost) // TODO: ok, here, how do I handle parallized regions if I want to be computationally efficeint?
    for (int k = 0; k < this->params.K; k++) {
        
        // initialize the cost
        cost = 0.0;

        // perform the rollout
        sol = this->dynamics.RK3_rollout(this->params.T_x, this->params.T_u, x0_sys, p0_feet, d0, U_bundle[k]);

        // compute the cost
        cost = this->cost_function(ref, sol, U_bundle[k]);
    
        // store the results (use critical sections to avoid race conditions if necessary)
        #pragma omp critical
        {
            Sol_bundle[k] = sol;
            J[k] = cost;
        }
    }

    // pack solutions into a tuple
    MC_Result mc;
    mc.S = Sol_bundle;
    mc.U = U_bundle;
    mc.J = J;

    // return the solutions
    return mc;
}


// select solutions based on cost
void Controller::sort_trajectories(const Solution_Bundle&  S,      const Vector_6d_Traj_Bundle& U,      const Vector_1d_List& J,
                                         Solution_Bundle&  S_elite,      Vector_6d_Traj_Bundle& U_elite,      Vector_1d_List& J_elite)
{
    // Create an index vector
    Vector_1i_List idx(this->params.K);
    std::iota(idx.begin(), idx.end(), 0);

    // Use nth_element to bring the top N_elite elements to the front in O(n) time
    std::nth_element(idx.begin(), idx.begin() + this->params.N_elite, idx.end(),
                     [&J](int i1, int i2) { return J[i1] < J[i2]; });

    // Populate elite solutions
    for (int i = 0; i < this->params.N_elite; i++) {
        S_elite[i] = S[idx[i]];
        U_elite[i] = U[idx[i]];
        J_elite[i] = J[idx[i]];
    }
}


// TODO: there is code optimziation work to be done here
// perform sampling predictive control of your choice here
RHC_Result Controller::sampling_predictive_control(double t_sim, Vector_12d x0_sys, Vector_6d p0_foot, Domain d0)
{
    // Monte Carlo Result to return
    MC_Result mc;

    // monte carlo variables
    Solution_Bundle S(this->params.N_elite), S_elite(this->params.N_elite);
    Vector_6d_Traj_Bundle U(this->params.N_elite), U_elite(this->params.N_elite);
    Vector_1d_List J(this->params.N_elite), J_elite(this->params.N_elite);

    // perform the CEM iterations
    auto t0_total = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < this->params.CEM_iters; i++) {

        // perform monte carlo simulation
        auto t0 = std::chrono::high_resolution_clock::now();
        mc = this->monte_carlo(t_sim, x0_sys, p0_foot, d0);

        // store monte carlo results
        S = mc.S;
        U = mc.U;
        J = mc.J;

        // sort the cost vector in ascending order
        this->sort_trajectories(S, U, J, S_elite, U_elite, J_elite);

        // update the distribution parameters
        this->update_distribution_params(U_elite);
        auto tf = std::chrono::high_resolution_clock::now();

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
                     [&J_elite](int i1, int i2) { return J_elite[i1] < J_elite[i2]; });

    // get the best solution
    Solution S_opt;
    Vector_6d_Traj U_opt;
    S_opt = S_elite[idx[0]];
    U_opt = U_elite[idx[0]];
    
    // return the best solution
    RHC_Result rhc;
    rhc.S = S_opt;
    rhc.U = U_opt;

    return rhc;
}
