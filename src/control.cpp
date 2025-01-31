#include "../inc/control.h"


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
        throw std::runtime_error("N_x must be greater than N_u");
    }
    if (this->params.N_elite > this->params.K) {
        throw std::runtime_error("N_elite must be less than or equal to K");
    }

    // build the cost matrices from the diagonal elements
    Vector_1d_List Q_com_diags_temp = config_file["COST"]["Q_com_diags"].as<std::vector<double>>();
    Vector_1d_List Q_leg_diags_temp = config_file["COST"]["Q_leg_diags"].as<std::vector<double>>();
    Vector_1d_List Qf_com_diags_temp = config_file["COST"]["Qf_com_diags"].as<std::vector<double>>();
    Vector_1d_List Qf_leg_diags_temp = config_file["COST"]["Qf_leg_diags"].as<std::vector<double>>();
    Vector_1d_List R_diags_temp = config_file["COST"]["R_diags"].as<std::vector<double>>();
    Vector_1d_List R_rate_diags_temp = config_file["COST"]["R_rate_diags"].as<std::vector<double>>();

    Vector_4d Q_com_diags, Qf_com_diags;
    Vector_4d Q_leg_diag, Qf_leg_diag;
    for (int i = 0; i < Q_com_diags_temp.size(); i++) {
        Q_com_diags[i] = Q_com_diags_temp[i];
        Qf_com_diags[i] = Qf_com_diags_temp[i];
        Q_leg_diag[i] = Q_leg_diags_temp[i];
        Qf_leg_diag[i] = Qf_leg_diags_temp[i];
    }
    Vector_8d Q_leg_diags, Qf_leg_diags;
    Q_leg_diags << Q_leg_diag, Q_leg_diag;
    Qf_leg_diags << Qf_leg_diag, Qf_leg_diag;
    
    Vector_2d R_leg_diags, R_rate_leg_diags;
    for (int i = 0; i < R_diags_temp.size(); i++) {
        R_leg_diags[i] = R_diags_temp[i];
        R_rate_leg_diags[i] = R_rate_diags_temp[i];
    }
    Vector_4d R_diags, R_rate_diags;
    R_diags << R_leg_diags, R_leg_diags;
    R_rate_diags << R_rate_leg_diags, R_rate_leg_diags;

    // initialize the cost matrices
    this->params.Q_com  = Q_com_diags.asDiagonal();
    this->params.Q_leg  = Q_leg_diags.asDiagonal();
    this->params.Qf_com = Qf_com_diags.asDiagonal();
    this->params.Qf_leg = Qf_leg_diags.asDiagonal();
    this->params.R = R_diags.asDiagonal();
    this->params.R_rate = R_rate_diags.asDiagonal();

    // gait cycle weight
    this->params.gait_cycle_weight = config_file["COST"]["gait_cycle_weight"].as<double>();

    // log barrier function
    this->params.log_barrier_enabled = config_file["COST"]["log_barrier_enabled"].as<bool>();
    this->params.r_weight = config_file["COST"]["r_weight"].as<double>();
    this->params.theta_weight = config_file["COST"]["theta_weight"].as<double>();
    this->params.rdot_weight = config_file["COST"]["rdot_weight"].as<double>();
    this->params.thetadot_weight = config_file["COST"]["thetadot_weight"].as<double>();

    // compute u(t) dt (N of the integration is not necessarily equal to the number of control points)
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


// construct the reference trajecotry
void Controller::initialize_reference_trajectories(YAML::Node config_file)
{
    // gait cycle parameters
    this->T_cycle = config_file["REFERENCE"]["T_cycle"].as<double>();
    this->T_SSP = config_file["REFERENCE"]["T_SSP"].as<double>();

    // reference parameters
    this->r_des = config_file["REFERENCE"]["r_des"].as<double>();
    this->theta_des = config_file["REFERENCE"]["theta_des"].as<double>();
    this->vx_des = config_file["REFERENCE"]["vx_des_"].as<double>();
    this->pz_des = config_file["REFERENCE"]["pz_des_"].as<double>();

    // set the desired COM trajectory
    Vector_1d_Traj time = config_file["REFERENCE"]["time"].as<std::vector<double>>();
    Vector_1d_Traj px_des_temp = config_file["REFERENCE"]["px_des"].as<std::vector<double>>();
    Vector_1d_Traj pz_des_temp = config_file["REFERENCE"]["pz_des"].as<std::vector<double>>();
    
    // ensure they are all the same size
    int N_ref = time.size();
    if (N_ref != px_des_temp.size() || N_ref != pz_des_temp.size()) {
        throw std::runtime_error("px_des, pz_des, and time must be the same size");
    }

    // set the reference trajectory for use inside the class
    Vector_2d_Traj p_com_ref(N_ref);
    for (int i = 0; i < N_ref; i++) {
        p_com_ref[i] << px_des_temp[i], pz_des_temp[i];
    }

    // set the reference trajectory for use inside the class
    this->t_ref = time;
    this->p_com_ref = p_com_ref;
}


// construct the intial distribution
void Controller::initialize_distribution(YAML::Node config_file)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg;
    int Nu = this->params.N_u;

    // initialize the matrices
    this->dist.mean.resize(2 * Nu * n_leg);
    this->dist.cov.resize(2 * Nu * n_leg, 2 * Nu * n_leg);
    this->dist.mean.setZero();
    this->dist.cov.setZero();

    // set the epsilon for numerical stability of covariance matrix
    this->dist.epsilon = config_file["DIST_PARAMS"]["epsilon"].as<double>();

    // compute the theoretical lower bound (lower bound of Frobenius norm, assume epsilon eigs)
    this->min_cov_norm = std::sqrt(2 * Nu * n_leg) * this->dist.epsilon;

    // set the initial mean
    std::vector<double> mean_temp = config_file["DIST_PARAMS"]["mu"].as<std::vector<double>>();
    Vector_4d mean;
    mean << mean_temp[0], mean_temp[1], mean_temp[0], mean_temp[1];
    for (int i = 0; i < Nu; i++) {
        this->dist.mean.segment<4>(2 * i * n_leg) = mean;
    }

    // set the initial covariance
    std::vector<double> cov_temp = config_file["DIST_PARAMS"]["sigma"].as<std::vector<double>>();
    Vector_4d cov_diags;
    cov_diags << cov_temp[0], cov_temp[1], cov_temp[0], cov_temp[1];
    Matrix_4d cov = cov_diags.asDiagonal();
    for (int i = 0; i < Nu; i++) {
        this->dist.cov.block<4, 4>(2 * i * n_leg, 2 * i * n_leg) = cov;
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
Vector_4d_Traj_Bundle Controller::sample_input_trajectory(int K)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg;
    int Nu = this->params.N_u;

    // sample the input trajectories
    Vector_d mu = this->dist.mean;
    Matrix_d Sigma = this->dist.cov;

    // perform cholesky decomposition 
    Eigen::LLT<Matrix_d> llt(Sigma);  
    Matrix_d L = llt.matrixL();  

    // check if the covariance is positive definite
    if (llt.info() == Eigen::NumericalIssue) {
        throw std::runtime_error("Covariance matrix is possibly not positive definite");
    }

    // Generate random input trajectories and store them in the input bundle
    Vector_d Z_vec, U_vec;
    Z_vec.resize(mu.size());
    U_vec.resize(mu.size());

    // U ~ N(mu, Sigma) <=> U = L * Z + mu; Z ~ N(0, I)
    // initialize the input trajectory bundle
    Vector_4d_Traj_Bundle U_bundle(K);
    Vector_4d_Traj U_traj(Nu);
    Vector_4d u;

    // loop over the number of samples
    for (int i = 0; i < K; i++) {
        
        // populate the Z vector; Z ~ N(0, I)
        for (int j = 0; j < mu.size(); j++) {
            Z_vec(j) = this->normal_dist(this->rand_generator);
        }

        // generate the vectorized input trajectory
        U_vec = L * Z_vec + mu;

        // unvectorize U_vec into U_traj
        for (int k = 0; k < this->params.N_u; k++) {
            u = U_vec.segment<4>(4 * k);
            U_traj[k] = u;
        }

        // store the input trajectory
        U_bundle[i] = U_traj;
    }

    return U_bundle;
}


// compute mean and covariance from a bundle of control inputs
void Controller::update_distribution_params(Vector_4d_Traj_Bundle U_bundle)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg;
    int Nu = this->params.N_u;

    // initialize the mean and covariance
    Vector_d mean;
    Matrix_d cov;
    mean.resize(2 * Nu * n_leg);
    cov.resize(2 * Nu * n_leg, 2 * Nu * n_leg);

    // size of the bundle
    int K = U_bundle.size(); // not necceseraly equal to this->params.K

    // used for computing the mean
    Matrix_d U_data;
    U_data.resize(2 * Nu * n_leg, K);

    // compute the mean
    Vector_4d_Traj U_traj(Nu);
    Vector_4d u;
    Vector_d U_vec;
    U_vec.resize(2 * Nu * n_leg);
    for (int i = 0; i < K; i++) {

        // vectorize the input trajectory
        U_traj = U_bundle[i];
        for (int j = 0; j < Nu; j++) {
            u = U_traj[j];
            U_vec.segment<4>(4 * j) = u;
        }   
        mean += U_vec;

        // insert into data matrix to use later
        U_data.col(i) = U_vec;
    }
    mean /= K;

    // compute the sample covariance (K-1 b/c Bessel correction)
    cov = (1.0 / (K-1)) * (U_data.colwise() - mean) * (U_data.colwise() - mean).transpose();
    
    // compute the eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<Matrix_d> eig(cov);
    if (eig.info() == Eigen::NumericalIssue) {
        throw std::runtime_error("Covariance matrix is possibly not positive definite");
    }
    Matrix_d eigvec = eig.eigenvectors();
    Vector_d eigval = eig.eigenvalues();
    Matrix_d eigvec_inv = eigvec.inverse();

    // modify eigenvalues with epsilon if it gets too low
    for (int i = 0; i < eigval.size(); i++) {
        eigval(i) = std::max(eigval(i), this->dist.epsilon);
    }

    // rebuild the covariance matrix with the eigenvalue decomposition, add epsilon to eigenvalues
    cov = eigvec * eigval.asDiagonal() * eigvec_inv;

    //  if stricly diagonal covariance option
    if (this->dist.diag_cov == true) {
        // set the covariance to be diagonal, (Hadamard Product, cov * I = diag(cov))
        cov = cov.cwiseProduct(Matrix_d::Identity(2 * Nu * n_leg, 2 * Nu * n_leg));
    }

    // update the distribution
    this->dist.mean = mean;
    this->dist.cov = cov;    
}


// generate a reference trajectory for the predictive control to track
Reference Controller::generate_reference_trajectory(double t_sim, Vector_4d x0_com)
{
    // pass reference to the dynamics
    Vector_4d_Traj X_com_ref(this->params.N_x);
    Vector_8d_Traj X_leg_ref(this->params.N_x);
    Vector_2i_Traj D_ref(this->params.N_x);

    // FIXME: when on single thread this segfaults. WHY?
    // get the "global time" reference trajectory
    Vector_1d_Traj t_ref_global = this->t_ref;
    Vector_2d_Traj p_com_ref_global = this->p_com_ref;
    int N_ref = t_ref_global.size();

    // build the COM reference trajectory
    Vector_4d xi_com_ref;
    Vector_2d pi_com_ref, p0_com_ref, pf_com_ref;
    Vector_2d vi_com_ref;
    double t, t0, tf;
    for (int i = 0; i < this->params.N_x; i++) {
    
        // compute the time given the global reference
        t = i * this->params.dt_x + t_sim;

        // find where t is in the time reference
        auto it = std::upper_bound(t_ref.begin(), t_ref.end(), t);
        int idx = std::distance(t_ref.begin(), it) - 1; // return -1 if before the first element
                                                        // return N_ref - 1 if after the last element

        // beyond last element
        if (idx == N_ref - 1) {
            // set to last point with zero velocity
            pi_com_ref = p_com_ref_global[N_ref - 1];
            vi_com_ref << 0.0, 0.0;
        }
        else{
            // linear interpolation
            t0 = t_ref_global[idx];
            tf = t_ref_global[idx + 1];
            p0_com_ref = p_com_ref_global[idx];
            pf_com_ref = p_com_ref_global[idx + 1];
            vi_com_ref = (pf_com_ref - p0_com_ref) / (tf - t0);
            pi_com_ref = p0_com_ref + vi_com_ref * (t - t0);
        }

        // build the reference
        xi_com_ref << pi_com_ref, vi_com_ref;

        // build the COM reference (only for hardcoded references)
        // t = i * this->params.dt_x;
        // xi_com_ref << x0_com(0) + this->vx_des * t, 
        //               this->pz_des, 
        //               this->vx_des, 
        //               0.0;

        // insert into trajectory
        X_com_ref[i] = xi_com_ref;
    }

    // build the LEG reference trajectory
    Vector_4d xi_leg_left_ref, xi_leg_right_ref;
    Vector_8d x_leg_ref;
    xi_leg_left_ref <<  this->r_des,  this->theta_des, 0.0, 0.0;
    xi_leg_right_ref << this->r_des, -this->theta_des, 0.0, 0.0;
    x_leg_ref << xi_leg_left_ref, xi_leg_right_ref;
    for (int i = 0; i < this->params.N_x; i++) {
        // insert into trajectory
        X_leg_ref[i] = x_leg_ref;
    }

    // build the domain reference trajectory
    Vector_2i d_ref;
    Vector_2i d;
    d_ref << 1, 1;                  // both feet in stance
    d << 1, 0;                      // both feet in stance
    double T_cycle = this->T_cycle; // time for SSP
    double T_SSP = this->T_SSP;     // time for SSP
    double T_reset = 0.0;           // time for reset
    for (int i = 0; i < this->params.N_x; i++) {
        
        // global time
        t = i * this->params.dt_x;

        // in one domain
        if ((t - T_reset) < T_SSP) {
            D_ref[i] = d;
        }
        // reset the domain
        else{
            d = d_ref - d;
            D_ref[i] = d;
            T_reset = t;
        }
    }

    // insert the references
    Reference ref;
    ref.X_com_ref = X_com_ref;
    ref.X_leg_ref = X_leg_ref;
    ref.D_ref = D_ref;

    return ref;
}


// evaluate log barrier function cost on legs
double Controller::cost_log_barrier(Vector_8d x_leg) {   
    
    // want to compute the cost of the log barrier function
    double J_log = 0.0;

    // log barrier function parameters
    double r_weight = this->params.r_weight;
    double rdot_weight = this->params.rdot_weight;
    double theta_weight = this->params.theta_weight;
    double thetadot_weight = this->params.thetadot_weight;

    // leg state constraints
    double r_min = this->dynamics.params.r_min;
    double r_max = this->dynamics.params.r_max;
    double theta_min = this->dynamics.params.theta_min;
    double theta_max = this->dynamics.params.theta_max;
    double rdot_lim = this->dynamics.params.rdot_lim;
    double thetadot_lim = this->dynamics.params.thetadot_lim;

    // compute the costs
    Vector_4d xi_leg;
    double r, theta, rdot, thetadot;
    double r_cost, rdot_cost, theta_cost, thetadot_cost;
    for (int i = 0; i < this->dynamics.n_leg; i++) {
        // unpack the leg state
        xi_leg = x_leg.segment<4>(4 * i);
        r = xi_leg(0);
        theta = xi_leg(1);
        rdot = xi_leg(2);
        thetadot = xi_leg(3);
        
        // compare the state to the limits

        // // compute the costs
        // r_cost = -r_weight * std::log(r - r_min) - r_weight * std::log(r_max - r);
        // rdot_cost = -rdot_weight * std::log(rdot + rdot_lim) - rdot_weight * std::log(rdot_lim - rdot);
        // theta_cost = -theta_weight * std::log(theta - theta_min) - theta_weight * std::log(theta_max - theta);
        // thetadot_cost = -thetadot_weight * std::log(thetadot + thetadot_lim) - thetadot_weight * std::log(thetadot_lim - thetadot);

        // add to the total cost
        // J_log += r_cost + rdot_cost + theta_cost + thetadot_cost;
    }

    return J_log =0.0;
}


// evaluate the cost function given a solution
double Controller::cost_function(const Reference& ref, const Solution& Sol, const Vector_4d_Traj& U)
{
    // trajectory length 
    int N = this->params.N_x;
    int Nu = this->params.N_u;

    // ************************************ COM COST ************************************

    // build the COM trajectory
    Vector_4d_Traj X_com(N);
    for (int i = 0; i < N; i++) {
        X_com[i] = Sol.x_sys_t[i].head<4>();
    }

    // compute the COM cost
    Vector_4d ei_com;
    double J_com = 0.0;
    // integrated cost
    for (int i = 0; i < N-1; i++) {
        ei_com = X_com[i] - ref.X_com_ref[i];
        J_com += ei_com.transpose() * this->params.Q_com * ei_com;
    }
    //terminal cost
    ei_com = X_com[N-1] - ref.X_com_ref[N-1];
    J_com += ei_com.transpose() * this->params.Qf_com * ei_com;

    // ************************************ LEG COST ************************************

    // compute the LEG cost
    Vector_8d ei_leg;
    double J_legs = 0.0;
    // integrated cost
    for (int i = 0; i < N-1; i++) {
        ei_leg = Sol.x_leg_t[i] - ref.X_leg_ref[i];
        J_legs += ei_leg.transpose() * this->params.Q_leg * ei_leg;
    }
    // terminal cost
    ei_leg = Sol.x_leg_t[N-1] - ref.X_leg_ref[N-1];
    J_legs += ei_leg.transpose() * this->params.Qf_leg * ei_leg;

    // log barrier function cost
    double J_legs_barrier = 0.0;
    if (this->params.log_barrier_enabled) {
        for (int i = 0; i < N; i++) {
            J_legs_barrier += this->cost_log_barrier(Sol.x_leg_t[i]);
        }
    }

    // ************************************ CONTACT CYCLE COST ************************************

    // TODO: Consider using a cost based on foot velocity and leg force
    // convert domain to a binary contact
    Vector_2i di;
    Domain d;
    double J_cycle = 0.0;
    double c;
    for (int i = 0; i < N; i++) {
        // convert contact type to binary contact
        d = Sol.domain_t[i];
        d[0] == Contact::STANCE ? di(0) = 1 : di(0) = 0;
        d[1] == Contact::STANCE ? di(1) = 1 : di(0) = 0;
        
        // if same as reference, no cost
        di == ref.D_ref[i] ? c = 0.0 : c = 1.0 * this->params.gait_cycle_weight;
        J_cycle += c;
    }

    // ************************************ INPUT COST ************************************

    // penalize the velocity input (want it to be closer to zero)
    Vector_4d ui;
    double J_input = 0.0;
    for (int i = 0; i < Nu; i++) {
        // left and right leg input vector
        ui << U[i];

        // compute quadratic cost
        J_input += ui.transpose() * this->params.R * ui;
    }

    // penalize the rate of change of the input (smaller rates of change)
    Vector_4d u_delta;
    for (int i = 0; i < Nu - 1; i++) {
        // left and right leg input vector
        u_delta = (U[i + 1] - U[i]) / this->params.dt_u;

        // compute quadratic cost
        J_input += (u_delta).transpose() * this->params.R_rate * (u_delta);
    }

    // ************************************ TOTAL COST ************************************

    // total cost
    double J_total; 
    J_total = J_com + J_legs + J_legs_barrier + J_input + J_cycle;

    return J_total;
}


// perform open loop rollouts
MC_Result Controller::monte_carlo(double t_sim, Vector_8d x0_sys, Vector_4d p0_feet, Domain d0, int K)
{
    // generate bundle of input trajectories
    Vector_4d_Traj_Bundle U_bundle(K);
    U_bundle = this->sample_input_trajectory(K);

    // initialize the containers for the solutions
    Solution_Bundle Sol_bundle(K);
    Vector_1d_List J(K);

    // generate the reference trajectory
    Reference ref;
    ref = this->generate_reference_trajectory(t_sim, x0_sys.head<4>());

    // loop over the input trajectories
    Solution sol;
    this->dynamics.resizeSolution(sol, this->params.T_x);
    double cost = 0.0;
    // #pragma omp parallel for private(sol, cost)
    for (int k = 0; k < K; k++) {
        
        // initialize the cost
        cost = 0.0;

        // perform the rollout
	
        // auto t0 = std::chrono::high_resolution_clock::now();
        this->dynamics.RK3_rollout(this->params.T_x, this->params.T_u, x0_sys, p0_feet, d0, U_bundle[k], sol);
        // auto tf = std::chrono::high_resolution_clock::now();
            // std::cout << "Time for rollout: " << std::chrono::duration<double, std::micro>(tf - t0).count() << " micros" << std::endl;

        // compute the cost
        // t0 = std::chrono::high_resolution_clock::now();
        cost = this->cost_function(ref, sol, U_bundle[k]);
        // tf = std::chrono::high_resolution_clock::now();
            // std::cout << "Time for cost: " << std::chrono::duration<double, std::micro>(tf - t0).count() << " micros" << std::endl;

        // store the results (use critical sections to avoid race conditions if necessary)
        // #pragma omp critical
        // {
            Sol_bundle[k] = sol;
            J[k] = cost;
        // }
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
void Controller::sort_trajectories(const Solution_Bundle&  S,       const Vector_4d_Traj_Bundle& U,        const Vector_1d_Traj& J,
                                   Solution_Bundle& S_elite, Vector_4d_Traj_Bundle& U_elite, Vector_1d_Traj& J_elite)
{
    int K = this->params.K;
    int N_elite = this->params.N_elite;

    // Create an index vector
    Vector_1i_List idx(K);
    std::iota(idx.begin(), idx.end(), 0);

    // Use nth_element to bring the top N_elite elements to the front in O(n) time
    std::nth_element(idx.begin(), idx.begin() + N_elite, idx.end(),
                     [&J](int i1, int i2) { return J[i1] < J[i2]; });

    // Populate elite solutions
    for (int i = 0; i < N_elite; i++) {
        S_elite[i] = S[idx[i]];
        U_elite[i] = U[idx[i]];
        J_elite[i] = J[idx[i]];
    }
}


// perform sampling predictive control of your choice here
RHC_Result Controller::sampling_predictive_control(double t_sim, Vector_8d x0_sys, Vector_4d p0_foot, Domain d0)
{
    // Monte Carlo Result to return
    MC_Result mc;

    // monte carlo variables
    Solution_Bundle S, S_elite(this->params.N_elite);
    Vector_4d_Traj_Bundle U, U_elite(this->params.N_elite);
    Vector_1d_List J, J_elite(this->params.N_elite);

    // perform the CEM iterations
    auto t0_total = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < this->params.CEM_iters; i++) {

        // perform monte carlo simulation
        auto t0 = std::chrono::high_resolution_clock::now();
        auto ts = std::chrono::high_resolution_clock::now();
        mc = this->monte_carlo(t_sim, x0_sys, p0_foot, d0, this->params.K);
        auto te = std::chrono::high_resolution_clock::now();
        // std::cout << "Time for monte carlo: " << std::chrono::duration<double, std::micro>(te - ts).count() << " micros" << std::endl;

        // store monte carlo results
        S = mc.S;
        U = mc.U;
        J = mc.J;

        // sort the cost vector in ascending order
        ts = std::chrono::high_resolution_clock::now();
        this->sort_trajectories(S, U, J, S_elite, U_elite, J_elite);
        te = std::chrono::high_resolution_clock::now();
        // std::cout << "Time for sort: " << std::chrono::duration<double, std::micro>(te - ts).count() << " micros" << std::endl;

        // update the distribution parameters
        ts = std::chrono::high_resolution_clock::now();
        this->update_distribution_params(U_elite);
        te = std::chrono::high_resolution_clock::now();
        // std::cout << "Time for update: " << std::chrono::duration<double, std::micro>(te - ts).count() << " micros" << std::endl;
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

    // return the best solution
    RHC_Result rhc;
    rhc.S = S_elite[0];
    rhc.U = U_elite[0];

    return rhc;
}

