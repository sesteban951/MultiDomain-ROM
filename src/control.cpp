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

    std::cout << "Q:\n" << this->params.Q << std::endl;
    std::cout << "Qf:\n" << this->params.Qf << std::endl;
    std::cout << "R:\n" << this->params.R << std::endl;

    // construct the initial distribution
    this->initialize_distribution(config_file);
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

    // set the initial mean
    std::vector<double> mean_temp = config_file["DIST_PARAMS"]["mu"].as<std::vector<double>>();
    Vector_2d mean_;
    Vector_4d mean;
    mean_ << mean_temp[0], mean_temp[1];
    mean << mean_, mean_;
    for (int i = 0; i < Nu; i++) {
        this->dist.mean.segment<n_leg * Nu>(2 * i * n_leg) = mean;
    }

    // set the initial covariance
    std::vector<double> cov_temp = config_file["DIST_PARAMS"]["sigma"].as<std::vector<double>>();
    Matrix_4d cov;
    cov << cov_temp[0], 0.0, 0.0, 0.0,
           0.0, cov_temp[1], ;
    for (int i = 0; i < this->params.Nu * this->dynamics.n_leg; i++) {
        this->dist.cov.block<2, 2 * this>(2 * i, 2 * i) = cov;
    }

    // set if covariance should be strictly diagonal
    this->dist.diag_cov = config_file["DIST_PARAMS"]["diag_cov"].as<bool>();

    // set the random 
    this->dist.seed = config_file["DIST_PARAMS"]["seed"].as<int>();
    this->dist.seed_enabled = config_file["DIST_PARAMS"]["seed_enabled"].as<bool>();

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
