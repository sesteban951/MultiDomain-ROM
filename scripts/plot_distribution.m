%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SOME SIM DISTRIBUTION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Load data
data_folder = "../data/3D/";
t = load(data_folder + 'time.csv');
mean_t = load(data_folder + 'mean.csv');
cov_t = load(data_folder + 'covariance.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove every even index
t = t(1:2:end);
mean_t = mean_t(1:2:end,:);
cov_t = cov_t(1:2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unapack varaibles from the YAML config
config_file = '../config/config_3D.yaml';
config = yaml.loadFile(config_file);

% get horizon and number of samples
N_sim = length(t);
N_u = config.CTRL_PARAMS.N_u;
n_legs = 2;
n_inputs = 3;

% shape of covariance matrix
cov_rows = n_legs * n_inputs * N_u;
cov_cols = cov_rows;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create tensor of covariances
cov = zeros(cov_rows, cov_cols, N_sim);
Sigma = zeros(cov_rows, cov_cols);
for i = 1:N_sim
    
    % build the covariance matrix
    cov_vec = cov_t(i,:);
    for j = 1:cov_rows
        for k = 1:cov_cols
            Sigma(j,k) = cov_vec((j-1)*cov_cols + k);
        end
    end

    % store the covariance matrix
    cov(:,:,i) = Sigma;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the first covariance matrix with heat map
figure;
colorbar;
c_min = min(cov(:));
c_max = max(cov(:));
caxis([c_min, c_max]); % Set constant color limits

idx = 1;
tic;
while idx < length(t)
    
    % Plot the covariance via heat map
    sigma_plot = imagesc(cov(:,:,idx));
    colorbar;
    caxis([c_min, c_max]); % Ensure the color limits stay the same

    msg = sprintf('Time = %.3f [sec]', t(idx));
    title(msg);
    drawnow;
    
    idx = idx + 1;
    while toc < t(idx)
        % wait
    end

    if idx < length(t)
        delete(sigma_plot);
    end
end