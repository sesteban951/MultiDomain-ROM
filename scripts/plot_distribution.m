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

% remove every other data point
% t = t(1:2:end);
% mean_t = mean_t(1:2:end,:);
% cov_t = cov_t(1:2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unapack varaibles from the YAML config
config_file = '../config/config_3D.yaml';
config = yaml.loadFile(config_file);

% get horizon and number of samples
N_sim = length(t);
N_x = config.CTRL_PARAMS.N_x;
N_u = config.CTRL_PARAMS.N_u;
dt_x = config.CTRL_PARAMS.dt_x;
dt_u = (N_x - 1) / (N_u - 1) * dt_x;
n_legs = 2;
n_inputs = 3;

% shape of covariance matrix
cov_rows = n_legs * n_inputs * N_u;
cov_cols = cov_rows;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animate_mean = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build the input time trajectory
t_input = 0 : dt_u : (N_u-1)*dt_u;

% build the mean trajectories
rdot_L = zeros(1, N_u, N_sim);
thetadot_x_L = zeros(1, N_u, N_sim);
thetadot_y_L = zeros(1, N_u, N_sim);
rdot_R = zeros(1, N_u, N_sim);
thetadot_x_R = zeros(1, N_u, N_sim);
thetadot_y_R = zeros(1, N_u, N_sim);
for i = 1:N_sim
    u_traj = mean_t(i,:);
    for j = 1:N_u
        u_t = u_traj(6*(j-1)+1 : 6*j);
        rdot_L(1,j,i) = u_t(1);
        thetadot_x_L(1,j,i) = u_t(2);
        thetadot_y_L(1,j,i) = u_t(3);
        rdot_R(1,j,i) = u_t(4);
        thetadot_x_R(1,j,i) = u_t(5);
        thetadot_y_R(1,j,i) = u_t(6);
    end    
end

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

if animate_mean == 1

    figure("Name", "Mean Trajectories");

    % create the subplots
    rdot_L_plot = subplot(2,3,1);
    thetadot_x_L_plot = subplot(2,3,2);
    thetadot_y_L_plot = subplot(2,3,3);
    rdot_R_plot = subplot(2,3,4);
    thetadot_x_R_plot = subplot(2,3,5);
    thetadot_y_R_plot = subplot(2,3,6);

    % initialize the plots
    rdot_L_lin = plot(rdot_L_plot, t_input, rdot_L(:,:,1), 'LineWidth', 1.5);
    thetadot_x_L_lin = plot(thetadot_x_L_plot, t_input, thetadot_x_L(:,:,1), 'LineWidth', 1.5);
    thetadot_y_L_lin = plot(thetadot_y_L_plot, t_input, thetadot_y_L(:,:,1), 'LineWidth', 1.5);
    rdot_R_lin = plot(rdot_R_plot, t_input, rdot_R(:,:,1), 'LineWidth', 1.5);
    thetadot_x_R_lin = plot(thetadot_x_R_plot, t_input, thetadot_x_R(:,:,1), 'LineWidth', 1.5);
    thetadot_y_R_lin = plot(thetadot_y_R_plot, t_input, thetadot_y_R(:,:,1), 'LineWidth', 1.5);

    % set the axis limits once
    rdot_L_plot.YLim = [min(rdot_L(:)), max(rdot_L(:))];
    thetadot_x_L_plot.YLim = [min(thetadot_x_L(:)), max(thetadot_x_L(:))];
    thetadot_y_L_plot.YLim = [min(thetadot_y_L(:)), max(thetadot_y_L(:))];
    rdot_R_plot.YLim = [min(rdot_R(:)), max(rdot_R(:))];
    thetadot_x_R_plot.YLim = [min(thetadot_x_R(:)), max(thetadot_x_R(:))];
    thetadot_y_R_plot.YLim = [min(thetadot_y_R(:)), max(thetadot_y_R(:))];

    % turn the grid on
    grid(rdot_L_plot, 'on');
    grid(thetadot_x_L_plot, 'on');
    grid(thetadot_y_L_plot, 'on');
    grid(rdot_R_plot, 'on');
    grid(thetadot_x_R_plot, 'on');
    grid(thetadot_y_R_plot, 'on');

    % set the labels with LaTeX interpreter
    xlabel(rdot_L_plot, 'Time [sec]', 'Interpreter', 'latex');
    ylabel(rdot_L_plot, '$\dot{r}_L$ [m/s]', 'Interpreter', 'latex');
    xlabel(thetadot_x_L_plot, 'Time [sec]', 'Interpreter', 'latex');
    ylabel(thetadot_x_L_plot, '$\dot{\theta}_{x,L}$ [rad/s]', 'Interpreter', 'latex');
    xlabel(thetadot_y_L_plot, 'Time [sec]', 'Interpreter', 'latex');
    ylabel(thetadot_y_L_plot, '$\dot{\theta}_{y,L}$ [rad/s]', 'Interpreter', 'latex');
    xlabel(rdot_R_plot, 'Time [sec]', 'Interpreter', 'latex');
    ylabel(rdot_R_plot, '$\dot{r}_L$ [m/s]', 'Interpreter', 'latex');
    xlabel(thetadot_x_R_plot, 'Time [sec]', 'Interpreter', 'latex');
    ylabel(thetadot_x_R_plot, '$\dot{\theta}_{x,R}$ [rad/s]', 'Interpreter', 'latex');
    xlabel(thetadot_y_R_plot, 'Time [sec]', 'Interpreter', 'latex');
    ylabel(thetadot_y_R_plot, '$\dot{\theta}_{y,R}$ [rad/s]', 'Interpreter', 'latex');

    % put an yline = 0
    yline(rdot_L_plot, 0, 'k');
    yline(thetadot_x_L_plot, 0, 'k');
    yline(thetadot_y_L_plot, 0, 'k');
    yline(rdot_R_plot, 0, 'k');
    yline(thetadot_x_R_plot, 0, 'k');
    yline(thetadot_y_R_plot, 0, 'k');

    % set to painters
    set(gcf, 'Renderer', 'painters');

    % Loop over time and update plots
    idx = 1;
    tic;
    while idx < length(t)
        
        % Update the plot data
        set(rdot_L_lin, 'YData', rdot_L(:,:,idx));
        set(thetadot_x_L_lin, 'YData', thetadot_x_L(:,:,idx));
        set(thetadot_y_L_lin, 'YData', thetadot_y_L(:,:,idx));
        set(rdot_R_lin, 'YData', rdot_R(:,:,idx));
        set(thetadot_x_R_lin, 'YData', thetadot_x_R(:,:,idx));
        set(thetadot_y_R_lin, 'YData', thetadot_y_R(:,:,idx));

        % Update the title with the current time
        msg = sprintf('Time = %.3f [sec]', t(idx));
        sgtitle(msg);

        drawnow;

        % Increment index and wait if necessary
        idx = idx + 1;
        if idx < length(t)
            while toc < t(idx)
                % Wait until the next time step
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if animate_mean == 0
    % Plot the first covariance matrix with heat map
    figure;
    colorbar;
    c_min = min(cov(:,:,2:end), [], 'all');
    c_max = max(cov(:,:,2:end), [], 'all');
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
end