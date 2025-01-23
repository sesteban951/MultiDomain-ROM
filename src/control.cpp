#include "../inc/control.h"


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
    std::vector<double> Q_com_diags_temp = config_file["CTRL_PARAMS"]["Q_com_diags"].as<std::vector<double>>();
    std::vector<double> Q_leg_diags_temp = config_file["CTRL_PARAMS"]["Q_leg_diags"].as<std::vector<double>>();
    std::vector<double> Qf_com_diags_temp = config_file["CTRL_PARAMS"]["Qf_com_diags"].as<std::vector<double>>();
    std::vector<double> Qf_leg_diags_temp = config_file["CTRL_PARAMS"]["Qf_leg_diags"].as<std::vector<double>>();
    std::vector<double> R_diags_temp = config_file["CTRL_PARAMS"]["R_diags"].as<std::vector<double>>();

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
    
    // desired state variables
    this->pz_des = config_file["STATE"]["pz_des"].as<double>();
    this->vx_des = config_file["STATE"]["vx_des"].as<double>();
    this->r_des = config_file["STATE"]["r_des"].as<double>();

    // construct the initial distribution
    this->initialize_distribution(config_file);
}


// construct the intial distribution
void Controller::initialize_distribution(YAML::Node config_file)
{
    // some useful ints to use
    const int n_leg = this->dynamics.n_leg;
    const int Nu = this->params.Nu;

    // initialize the matrices
    this->dist.mean.resize(2 * Nu * n_leg);
    this->dist.cov.resize(2 * Nu * n_leg, 2 * Nu * n_leg);
    this->dist.mean.setZero();
    this->dist.cov.setZero();

    // set the epsilon for numerical stability of covariance matrix
    this->dist.epsilon = config_file["DIST_PARAMS"]["epsilon"].as<double>();

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
Vector_2d_Traj_Bundle Controller::sample_input_trajectory(int K)
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
    Vector_4d U_t;
    Vector_2d_List U_t_list(n_leg);
    Vector_2d_Traj U_traj(Nu);
    Z_vec.resize(mu.size());
    U_vec.resize(mu.size());

    // U ~ N(mu, Sigma) <=> U = L * Z + mu; Z ~ N(0, I)
    // initialize the input trajectory bundle
    Vector_2d_Traj_Bundle U_bundle(K);
    for (int i = 0; i < K; i++) {
        
        // populate the Z vector; Z ~ N(0, I)
        for (int i = 0; i < mu.size(); i++) {
            Z_vec(i) = this->normal_dist(this->rand_generator);
        }

        // generate the vectorized input trajectory
        U_vec = L * Z_vec + mu;

        // unvectorize U_vec into U_traj
        for (int j = 0; j < this->params.Nu; j++) {
            U_t = U_vec.segment<4>(4 * j);
            U_t_list[0] = U_t.segment<2>(0);
            U_t_list[1] = U_t.segment<2>(2);
            U_traj[j] = U_t_list;
        }

        // store the input trajectory
        U_bundle[i] = U_traj;
    }

    return U_bundle;
}


// compute mean and covariance from a bundle of control inputs
void Controller::update_distribution_params(Vector_2d_Traj_Bundle U_bundle)
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
    Vector_2d_Traj U_traj;
    Vector_2d_List U_t;
    Vector_4d U_t_vec;
    Vector_d U_traj_vec;
    U_traj_vec.resize(2 * Nu * n_leg);
    for (int i = 0; i < K; i++) {

        // vectorize the input trajectory
        U_traj = U_bundle[i];
        for (int j = 0; j < Nu; j++) {
            U_t = U_traj[j];
            U_t_vec << U_t[0], U_t[1];
            U_traj_vec.segment<4>(4 * j) = U_t_vec;
        }   
        mean += U_traj_vec;

        // insert into data matrix to use later
        U_data.col(i) = U_traj_vec;
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
        cov = cov.cwiseProduct(Matrix_d::Identity(this->params.Nu * 2, this->params.Nu * 2));
    }

    // update the distribution
    this->dist.mean = mean;
    this->dist.cov = cov;    
}

