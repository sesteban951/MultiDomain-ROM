#include "../inc/control.h"


// Constructor
Controller::Controller(YAML::Node config_file) : dynamics(config_file)
{
    // set the control parameters
    this->params.N = config_file["CTRL_PARAMS"]["N"].as<int>();
    this->params.dt = config_file["CTRL_PARAMS"]["dt"].as<double>();
    this->params.K = config_file["CTRL_PARAMS"]["K"].as<int>();
    this->params.Nu = config_file["CTRL_PARAMS"]["Nu"].as<int>();
    this->params.N_elite = config_file["CTRL_PARAMS"]["N_elite"].as<int>();
    this->params.CEM_iters = config_file["CTRL_PARAMS"]["CEM_iters"].as<int>();
    
    // build the cost matrices from the diagonal elements
    Vector_1d_List Q_com_diags_temp = config_file["COST"]["Q_com_diags"].as<std::vector<double>>();
    Vector_1d_List Q_leg_diags_temp = config_file["COST"]["Q_leg_diags"].as<std::vector<double>>();
    Vector_1d_List Qf_com_diags_temp = config_file["COST"]["Qf_com_diags"].as<std::vector<double>>();
    Vector_1d_List Qf_leg_diags_temp = config_file["COST"]["Qf_leg_diags"].as<std::vector<double>>();
    Vector_1d_List R_diags_temp = config_file["COST"]["R_diags"].as<std::vector<double>>();

    Vector_4d Q_com_diags, Q_leg_diags, Qf_com_diags, Qf_leg_diags;
    for (int i = 0; i < Q_com_diags_temp.size(); i++) {
        Q_com_diags[i] = Q_com_diags_temp[i];
        Q_leg_diags[i] = Q_leg_diags_temp[i];
        Qf_com_diags[i] = Qf_com_diags_temp[i];
        Qf_leg_diags[i] = Qf_leg_diags_temp[i];
    }
    Vector_12d Q_diags, Qf_diags;
    Q_diags << Q_com_diags, Q_leg_diags, Q_leg_diags;
    Qf_diags << Qf_com_diags, Qf_leg_diags, Qf_leg_diags;
    
    Vector_2d R_leg_diags;
    for (int i = 0; i < R_diags_temp.size(); i++) {
        R_leg_diags[i] = R_diags_temp[i];
    }
    Vector_4d R_diags;
    R_diags << R_leg_diags, R_leg_diags;
    
    // initialize the cost matrices
    this->params.Q  = Q_diags.asDiagonal();
    this->params.Qf = Qf_diags.asDiagonal();
    this->params.R  = R_diags.asDiagonal();

    // construct the reference trajectory
    this->initialize_reference_trajectories(config_file);

    // construct the initial distribution
    this->initialize_distribution(config_file);

    // set the number of parallel threads to use
    bool threading_enabled = config_file["THREADING"]["enabled"].as<bool>();
    if (threading_enabled == true) {

        // set the number of threads
        int num_threads = config_file["THREADING"]["num_threads"].as<int>();
        omp_set_num_threads(num_threads); // Set the number of threads for parallel regions
        
        // enable nested parallelism
        bool nested = config_file["THREADING"]["nested"].as<bool>();
        if (nested == true) { omp_set_nested(1); }
    }
}


// construct the reference trajecotry
void Controller::initialize_reference_trajectories(YAML::Node config_file)
{
    // reference parameters
    this->r_des = config_file["REFERENCE"]["r_des"].as<double>();
    this->theta_des = config_file["REFERENCE"]["theta_des"].as<double>();
    this->vx_des = config_file["REFERENCE"]["vx_des_"].as<double>();
    this->pz_des = config_file["REFERENCE"]["pz_des_"].as<double>();

    // set the desired COM trajectory
    Vector_1d_List px_des_temp = config_file["REFERENCE"]["px_des"].as<std::vector<double>>();
    Vector_1d_List pz_des_temp = config_file["REFERENCE"]["pz_des"].as<std::vector<double>>();
    Vector_1d_Traj time = config_file["REFERENCE"]["time"].as<std::vector<double>>();

    Vector_2d_Traj p_com_traj(time.size());
    Vector_2d p_com_des_;
    int N_ref = time.size();
    for (int i = 0; i < N_ref; i++) {
        p_com_des_ << px_des_temp[i], pz_des_temp[i];
        p_com_traj[i] = p_com_des_;
    }

    // some time stuff for the reference trajectory
    int N = this->params.N;
    double dt = this->params.dt;
    double rhc_time = (N-1) * dt;         // receding horizon total time
    double t_end = time[N_ref - 1];       // reference traj total time
    int N_ref_dt = std::ceil(rhc_time / dt); // number of knots in the reference trajectory

    // build the reference trajectory
    Vector_1d_Traj time_ref(N_ref_dt);
    Vector_4d_Traj X_com_ref(N_ref_dt);
    Vector_2d p_com_des, v_com_des;
    Vector_4d x_com_des;
    double t, t0, tf;
    for (int i = 0; i < N_ref_dt; i++) {
        // build the time array
        t = i * dt;
        time_ref[i] = t;

        // find where the time is in the reference trajectory
        auto it = std::upper_bound(time.begin(), time.end(), t);
        int idx = std::distance(time.begin(), it) - 1; // return -1 if before the first element
                                                       // return N_ref - 1 if after the last element

        // beyond last element
        if (idx == N_ref - 1) {
            p_com_des = p_com_traj[N_ref - 1]; // position is the last element
            v_com_des << 0.0, 0.0;            // velocity is zero
        }
        else{
            t0 = time[idx];
            tf = time[idx + 1];
            Vector_2d p0 = p_com_traj[idx];
            Vector_2d pf = p_com_traj[idx + 1];
            v_com_des = (pf - p0) / (tf - t0);
            p_com_des = p0 + v_com_des * (t - t0);
        }

        // populate the desired COM state
        x_com_des << p_com_des, v_com_des;
        X_com_ref[i] = x_com_des;
    }
}


// construct the intial distribution
void Controller::initialize_distribution(YAML::Node config_file)
{
    // some useful ints to use
    int n_leg = this->dynamics.n_leg;
    int Nu = this->params.Nu;

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
    int Nu = this->params.Nu;

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
        for (int k = 0; k < this->params.Nu; k++) {
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
    int Nu = this->params.Nu;

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
Vector_12d_Traj Controller::generate_reference_trajectory(Vector_4d x0_com)
{
    // pass reference to the dynamics
    Vector_12d_Traj X_ref;
    X_ref.resize(this->params.N);

    // populate the reference trajectory
    Vector_12d xi_ref;
    Vector_4d xi_com_ref, xi_leg_left_ref, xi_leg_right_ref;
    xi_leg_left_ref <<  this->r_des,  this->theta_des, 0.0, 0.0;
    xi_leg_right_ref << this->r_des, -this->theta_des, 0.0, 0.0;
    
    for (int i = 0; i < this->params.N; i++) {
        
        // build the COM reference    
        xi_com_ref << x0_com(0) + this->vx_des * i * this->params.dt, 
                      this->pz_des, 
                      this->vx_des, 
                      0.0;

        // build full state reference
        xi_ref << xi_com_ref, xi_leg_left_ref, xi_leg_right_ref;

        // insert into trajectory
        X_ref[i] = xi_ref;
    }

    return X_ref;
}


// evaluate the cost function given a solution
double Controller::cost_function(Vector_12d_Traj X_ref, Solution Sol, Vector_4d_Traj U)
{
    // trajectory length 
    int N = this->params.N;
    int Nu = this->params.Nu;

    // upack the relevant variables
    Vector_8d_Traj X_sys = Sol.x_sys_t;
    Vector_8d_Traj X_leg = Sol.x_leg_t;
    // Vector_8d_Traj X_foot = Sol.x_foot_t; // TODO: add some kind of cost for the foot state

    // convert to the trajectories into big matrices
    Matrix_d X_sys_mat, X_leg_mat;
    X_sys_mat.resize(8, N);
    X_leg_mat.resize(8, N);
    for (int i = 0; i < N; i++) {
        // populate the system matrix
        X_sys_mat.col(i) = X_sys[i];

        // populate the leg matrix
        X_leg_mat.col(i) = X_leg[i];
    }

    // combine the COM state with the LEG state
    Matrix_d X_com;
    X_com.resize(4, N);
    X_com = X_sys_mat.block(0, 0, 4, N);
    
    Matrix_d X;
    X.resize(12, N);
    X << X_com, X_leg_mat;

    // state cost
    Vector_12d xi, xi_des, ei;
    double J_state, cost;
    J_state = 0.0;
    for (int i = 0; i < N; i++) {
        // compute error state
        xi = X.col(i);
        xi_des = X_ref[i];
        ei = (xi - xi_des);

        // compute the stage cost
        cost = ei.transpose() * this->params.Q * ei;
        J_state += cost;
    }

    // input cost
    Vector_4d ui;
    double J_input = 0.0;
    for (int i = 0; i < Nu; i++) {
        // left and right leg input vector
        ui << U[i];

        // compute quadratic cost
        J_input += ui.transpose() * this->params.R * ui;
    }

    // terminal cost
    xi = X.col(N - 1);
    xi_des = X_ref[N - 1];
    ei = (xi - xi_des);
    cost = ei.transpose() * this->params.Qf * ei;
    J_state += cost;

    // total cost
    double J_total; 
    J_total = J_state + J_input;

    return J_total;
}


// perform open loop rollouts
MC_Result Controller::monte_carlo(Vector_8d x0_sys, Vector_4d p0_feet, Domain d0, int K)
{
    // compute u(t) dt (N of the integration is not necessarily equal to the number of control points)
    double T = (this->params.N-1) * this->params.dt;
    double dt_u = T / (this->params.Nu-1);

    // generate the time arrays
    Vector_1d_Traj T_x;
    Vector_1d_Traj T_u;
    T_x.resize(this->params.N);
    T_u.resize(this->params.Nu);
    for (int i = 0; i < this->params.N; i++) {
        T_x[i] = i * this->params.dt;
    }
    for (int i = 0; i < this->params.Nu; i++) {
        T_u[i] = i * dt_u;
    }

    // generate bundle of input trajectories
    Vector_4d_Traj_Bundle U_bundle(K);
    U_bundle = this->sample_input_trajectory(K);

    // initialize the containers for the solutions
    Solution_Bundle Sol_bundle(K);
    Vector_1d_List J(K);

    // generate the reference trajectory
    Vector_12d_Traj X_ref;
    X_ref = this->generate_reference_trajectory(x0_sys.head<4>());

    // loop over the input trajectories
    Solution sol;
    double cost = 0.0;
    #pragma omp parallel for private(sol, cost)
    for (int k = 0; k < K; k++) {
        
        // initialize the solution and cost
        sol = Solution();
        cost = 0.0;

        // perform the rollout
        sol = this->dynamics.RK3_rollout(T_x, T_u, x0_sys, p0_feet, d0, U_bundle[k]);

        // compute the cost
        cost = this->cost_function(X_ref, sol, U_bundle[k]);

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
void Controller::sort_trajectories(Solution_Bundle  S,       Vector_4d_Traj_Bundle U,        Vector_1d_Traj J,
                                   Solution_Bundle& S_elite, Vector_4d_Traj_Bundle& U_elite, Vector_1d_Traj& J_elite)
{
    // sort the cost vector in ascending order
    Vector_1i_List idx(this->params.K);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&J](int i1, int i2) {return J[i1] < J[i2];});

    // elite solution containers
    Solution_Bundle S_elite_(this->params.N_elite);
    Vector_4d_Traj_Bundle U_elite_(this->params.N_elite);
    Vector_1d_Traj J_elite_(this->params.N_elite);

    // populate the elite solutions
    for (int i = 0; i < this->params.N_elite; i++) {
        S_elite_[i] = S[idx[i]];  // top solutions
        U_elite_[i] = U[idx[i]];  // top control inputs
        J_elite_[i] = J[idx[i]];  // top costs
    }

    // update the elite solutions
    S_elite = S_elite_;
    U_elite = U_elite_;
    J_elite = J_elite_;
}


// perform sampling predictive control of your choice here
Solution Controller::sampling_predictive_control(Vector_8d x0_sys, Vector_4d p0_foot, Domain d0)
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
        mc = this->monte_carlo(x0_sys, p0_foot, d0, this->params.K);

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
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "CEM Iteration: " << i+1 << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "Time for iteration: " << std::chrono::duration<double, std::milli>(tf - t0).count() << " ms" << std::endl;
        std::cout << "Smallest cost: " << J_elite[0] << std::endl;   
        std::cout << "Norm of covariance: " << this->dist.cov.norm() << ", Theoretical min: " << this->min_cov_norm << std::endl;
        std::cout << std::endl;

    }
    auto tf_total = std::chrono::high_resolution_clock::now();
    
    double T_tot = std::chrono::duration<double>(tf_total - t0_total).count();

    std::cout << "CEM complete" << std::endl;
    std::cout << "Total time: " << T_tot << " sec" << std::endl;
    std::cout << "Average time per iteration: " << T_tot / this->params.CEM_iters << " [sec], " << this->params.CEM_iters/T_tot << " [Hz]" << std::endl;
    std::cout << "Average Rollout time: " << T_tot / (this->params.CEM_iters * this->params.K) * 1000000.0 << " [us]" << std::endl << std::endl;

    // return the best solution
    return S_elite[0];
}

