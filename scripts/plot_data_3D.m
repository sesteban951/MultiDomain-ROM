%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SOME SIM RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Load data
data_folder = "../data/3D/";

t = load(data_folder + 'time.csv');
x_sys = load(data_folder + 'state_sys.csv');
x_leg = load(data_folder + 'state_leg.csv');
x_foot = load(data_folder + 'state_foot.csv');
u = load(data_folder + 'input.csv');
lambd = load(data_folder + 'lambda.csv');
tau = load(data_folder + 'tau.csv');
d = load(data_folder + 'domain.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% downsample
downsample = 0;
hz = 40;

% plotting / animation
animate = 1;   % animation = 1; plot states = 0
rt = 1.0;      % realtime rate
perspective = '';    % 'T'op, 'F'ront, 'S'ide
replays = 2;   % how many times to replay the animation
plot_com = 1;  % plot the foot trajectory
plot_foot = 1; % plot the foot trajectory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% segment the time
t_interval = [t(1) t(end)];
% t_interval = [0 0.25];

% apply time window
idx = find(t >= t_interval(1) & t <= t_interval(2));
t = t(idx);
x_sys = x_sys(idx,:);
x_leg = x_leg(idx,:);
x_foot = x_foot(idx,:);
u = u(idx,:);
lambd = lambd(idx,:);
tau = tau(idx,:);
d = d(idx,:);

% frequency of the data
if downsample == 1
    dt_data = t(2) - t(1);
    nth_sample = round(1/(hz * dt_data));

    % downsample the data
    t = downsample(t, nth_sample);
    x_sys = downsample(x_sys, nth_sample);
    x_leg = downsample(x_leg, nth_sample);
    x_foot = downsample(x_foot, nth_sample);
    u = downsample(u, nth_sample);
    lambd = downsample(lambd, nth_sample);
    tau = downsample(tau, nth_sample);
    d = downsample(d, nth_sample);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unapack varaibles from the YAML config
config_file = '../config/config_3D.yaml';
config = yaml.loadFile(config_file);

% some time parameters
dt = config.CTRL_PARAMS.dt_x;

% some kineamtic limits and velocities
r_min = config.SYS_PARAMS.r_min;
r_max = config.SYS_PARAMS.r_max;
theta_x_min = config.SYS_PARAMS.theta_x_min;
theta_x_max = config.SYS_PARAMS.theta_x_max;
theta_y_min = config.SYS_PARAMS.theta_y_min;
theta_y_max = config.SYS_PARAMS.theta_y_max;
rdot_lim = config.SYS_PARAMS.rdot_lim;
thetadot_lim = config.SYS_PARAMS.thetadot_lim;

% extract some variables
interp = config.SYS_PARAMS.interp;
friction_coeff = config.SYS_PARAMS.friction_coeff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% system state
p_com = x_sys(:,1:3);
v_com = x_sys(:,4:6);
x_leg_commands_L = x_sys(:,7:9);
x_leg_commands_R = x_sys(:,10:12);

% leg states
x_leg_L = x_leg(:,1:6);
x_leg_R = x_leg(:,7:12);

% foot states
x_foot_L = x_foot(:,1:6);
x_foot_R = x_foot(:,7:12);

% inputs
u_L = u(:,1:3);
u_R = u(:,4:6);
udot_L = diff(u_L) / dt;
udot_R = diff(u_R) / dt;

% lambda leg forces
lambd_L = lambd(:,1:3);
lambd_R = lambd(:,4:6);
lambd_L_norm = zeros(length(t), 1);
lambd_R_norm = zeros(length(t), 1);
for i = 1:length(t)
    lambd_L_norm(i) = norm(lambd_L(i,:));
    lambd_R_norm(i) = norm(lambd_R(i,:));
end

% ankle torques
tau_L = tau(:,1:3);
tau_R = tau(:,4:6);

% domain
d_L = d(:,1);
d_R = d(:,2);

% friction cones
lambda_x_L = lambd_L(:,1);
lambda_y_L = lambd_L(:,2);
lambda_z_L = lambd_L(:,3);
lambda_x_R = lambd_R(:,1);
lambda_y_R = lambd_R(:,2);
lambda_z_R = lambd_R(:,3);

z_max = max([lambda_z_L; lambda_z_R]);
height = z_max;
r = friction_coeff * z_max;
num_pts = 100;
theta = linspace(0, 2*pi, num_pts);
x_cone = r * cos(theta);
y_cone = r * sin(theta);
z_cone = height * ones(1, num_pts);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if animate == 0

    fig = figure('Name', 'States', 'WindowState', 'maximized');
    set(gcf,'renderer','painters');
    tabgroup = uitabgroup(fig);

    tab1 = uitab(tabgroup, 'Title', 'States');
    axes('Parent', tab1);

    % COM STATES
    subplot(4,6,1);
    hold on; grid on;
    plot(t, p_com(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_x$ [m]', 'Interpreter', 'latex');
    title('COM x-pos');

    subplot(4,6,2);
    hold on; grid on;
    plot(t, p_com(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_y$ [m]', 'Interpreter', 'latex');
    title('COM y-pos');

    subplot(4,6,3); 
    hold on; grid on;
    plot(t, p_com(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_z$ [m/s]', 'Interpreter', 'latex');
    title('COM z-pos');

    subplot(4,6,7);
    hold on; grid on;
    plot(t, v_com(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_x$ [m/s]', 'Interpreter', 'latex');
    title('COM x-vel');

    subplot(4,6,8);
    hold on; grid on;
    plot(t, v_com(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_y$ [m/s]', 'Interpreter', 'latex');
    title('COM y-vel');

    subplot(4,6,9);
    hold on; grid on;
    plot(t, v_com(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_z$ [m/s]', 'Interpreter', 'latex');
    title('COM z-vel');

    % FOOT STATES
    subplot(4,6,4);
    hold on; grid on;
    plot(t, x_foot_L(:,1), 'LineWidth', 2);
    plot(t, x_foot_R(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_{x}$ [m]', 'Interpreter', 'latex');
    title('FOOT x-pos');
    % legend('L', 'R');

    subplot(4,6,5);
    hold on; grid on;
    plot(t, x_foot_L(:,2), 'LineWidth', 2);
    plot(t, x_foot_R(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_{y}$ [m]', 'Interpreter', 'latex');
    title('FOOT y-pos');
    % legend('L', 'R');

    subplot(4,6,6);
    hold on; grid on;
    plot(t, x_foot_L(:,3), 'LineWidth', 2);
    plot(t, x_foot_R(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_{z}$ [m]', 'Interpreter', 'latex');
    title('FOOT z-pos');
    % legend('L', 'R');

    subplot(4,6,10);
    hold on; grid on;
    plot(t, x_foot_L(:,4), 'LineWidth', 2);
    plot(t, x_foot_R(:,4), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_{x}$ [m/s]', 'Interpreter', 'latex');
    title('FOOT x-vel');
    % legend('L', 'R');

    subplot(4,6,11);
    hold on; grid on;
    plot(t, x_foot_L(:,5), 'LineWidth', 2);
    plot(t, x_foot_R(:,5), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_{y}$ [m/s]', 'Interpreter', 'latex');
    title('FOOT y-vel');
    % legend('L', 'R');

    subplot(4,6,12);
    hold on; grid on;
    plot(t, x_foot_L(:,6), 'LineWidth', 2);
    plot(t, x_foot_R(:,6), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_{z}$ [m/s]', 'Interpreter', 'latex');
    title('FOOT z-vel');
    % legend('L', 'R');

    % LEG STATES
    subplot(4,6,13);
    hold on; grid on;
    plot(t, x_leg_L(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    % plot(t, x_leg_commands_L(:,1), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,1), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    % plot(t, x_leg_commands_R(:,1), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    yline(r_min, '--', 'Min');
    yline(r_max, '--', 'Max');
    xlabel('Time [sec]');
    ylabel('$r$ [m]', 'Interpreter', 'latex');
    title('LEG Length');
    % legend('$r_L$', '$\hat{r}_L$', '$r_R$', '$\hat{r}_R$', 'Interpreter', 'latex');

    subplot(4,6,14);
    hold on; grid on;
    plot(t, x_leg_L(:,2), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    % plot(t, x_leg_commands_L(:,2), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    % plot(t, x_leg_commands_R(:,2), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    yline(theta_x_min, '--', 'Min');
    yline(theta_x_max, '--', 'Max');
    yline(-theta_x_min, '--', 'Min');
    yline(-theta_x_max, '--', 'Max');
    xlabel('Time [sec]');
    ylabel('$\theta_x$ [rad]', 'Interpreter', 'latex');
    title('LEG Angle x-axis');
    % legend('$\theta_{x,L}$', '$\hat{\theta}_{x,L}$', '$\theta_{x,R}$', '$\hat{\theta}_{x,R}$', 'Interpreter', 'latex');

    subplot(4,6,15);
    hold on; grid on;
    plot(t, x_leg_L(:,3), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    % plot(t, x_leg_commands_L(:,3), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,3), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    yline(theta_y_min, '--', 'Min');
    yline(theta_y_max, '--', 'Max');
    % plot(t, x_leg_commands_R(:,3), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    xlabel('Time [sec]');
    ylabel('$\theta_y$ [rad]', 'Interpreter', 'latex');
    title('LEG Angle y-axis');
    % legend('$\theta_{y,L}$', '$\hat{\theta}_{y,L}$', '$\theta_{y,R}$', '$\hat{\theta}_{y,R}$', 'Interpreter', 'latex');

    subplot(4,6,[19:20]);
    hold on; grid on;
    plot(t, x_leg_L(:,4), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    % plot(t, u_L(:,1), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,4), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    % plot(t, u_R(:,1), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    yline(-rdot_lim, '--', 'Min');
    yline(rdot_lim, '--', 'Max');
    xlabel('Time [sec]');
    ylabel('$\dot{r}$ [m/s]', 'Interpreter', 'latex');
    title('LEG Length vel');
    % legend('$\dot{r}_L$', '$\hat{\dot{r}}_L$', '$\dot{r}_R$', '$\hat{\dot{r}}_R$', 'Interpreter', 'latex');

    subplot(4,6,[21:22]);
    hold on; grid on;
    plot(t, x_leg_L(:,5), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    % plot(t, u_L(:,2), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,5), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    % plot(t, u_R(:,2), 'LineWideth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    yline(-thetadot_lim, '--', 'Min');
    yline(thetadot_lim, '--', 'Max');
    xlabel('Time [sec]');
    ylabel('$\dot{\theta}_x$ [rad/s]', 'Interpreter', 'latex');
    title('LEG Vel x-axis');
    % legend('$\dot{\theta}_{x,L}$', '$\hat{\dot{\theta}}_{x,L}$', '$\dot{\theta}_{x,R}$', '$\hat{\dot{\theta}}_{x,R}$', 'Interpreter', 'latex');

    subplot(4,6,[23:24]);
    hold on; grid on;
    plot(t, x_leg_L(:,6), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    % plot(t, u_L(:,3), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,6), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    % plot(t, u_R(:,3), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    yline(-thetadot_lim, '--', 'Min');
    yline(thetadot_lim, '--', 'Max');
    xlabel('Time [sec]');
    ylabel('$\dot{\theta}_y$ [rad/s]', 'Interpreter', 'latex');
    title('LEG Vel y-axis');
    % legend('$\dot{\theta}_{y,L}$', '$\hat{\dot{\theta}}_{y,L}$', '$\dot{\theta}_{y,R}$', '$\hat{\dot{\theta}}_{y,R}$', 'Interpreter', 'latex');

    % DOMAIN
    subplot(4,6,[16:18]);
    hold on; grid on;
    stairs(t, d_L, 'LineWidth', 2);
    stairs(t, d_R, 'LineWidth', 1.5);
    xlabel('Time [sec]');
    title('Domain');
    ylim([-0.5, 1.5]);
    yticks([0, 1]);
    yticklabels({'F', 'G'});
    % legend('L', 'R');

    tab2 = uitab(tabgroup, 'Title', 'Inputs');
    axes('Parent', tab2);

    % INPUT
    subplot(3,2,1); 
    hold on; grid on;
    if interp == 'Z'
        stairs(t, u_L(:,1), 'LineWidth', 2);
        stairs(t, u_R(:,1), 'LineWidth', 2);
    elseif interp == 'L'
        plot(t, u_L(:,1), 'LineWidth', 2);
        plot(t, u_R(:,1), 'LineWidth', 2);
    end
    xlabel('Time [sec]');
    ylabel('$\dot{l_0}$', 'interpreter', 'latex');
    title('LEG l_0 ');
    legend('L', 'R');

    subplot(3,2,3); 
    hold on; grid on;
    if interp == 'Z'
        stairs(t, u_L(:,2), 'LineWidth', 2);
        stairs(t, u_R(:,2), 'LineWidth', 2);
    elseif interp == 'L'
        plot(t, u_L(:,2), 'LineWidth', 2);
        plot(t, u_R(:,2), 'LineWidth', 2);
    end
    xlabel('Time [sec]');
    ylabel('$\dot{\theta}_x$', 'interpreter', 'latex');
    title('LEG angle x');
    legend('L', 'R');

    subplot(3,2,5); 
    hold on; grid on;
    if interp == 'Z'
        stairs(t, u_L(:,3), 'LineWidth', 2);
        stairs(t, u_R(:,3), 'LineWidth', 2);
    elseif interp == 'L'
        plot(t, u_L(:,3), 'LineWidth', 2);
        plot(t, u_R(:,3), 'LineWidth', 2);
    end
    xlabel('Time [sec]');
    ylabel('$\dot{\theta}_y$', 'interpreter', 'latex');
    title('LEG angle y');
    legend('L', 'R');

    subplot(3,2,2);
    hold on; grid on;
    plot(udot_L(:,1), 'LineWidth', 2);
    plot(udot_R(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\ddot{l_0}$', 'interpreter', 'latex');
    title('LEG l_0 accel');
    legend('L', 'R');

    subplot(3,2,4);
    hold on; grid on;
    plot(udot_L(:,2), 'LineWidth', 2);
    plot(udot_R(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\ddot{\theta}_x$', 'interpreter', 'latex');
    title('LEG angle x accel');
    legend('L', 'R');

    subplot(3,2,6);
    hold on; grid on;
    plot(udot_L(:,3), 'LineWidth', 2);
    plot(udot_R(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\ddot{\theta}_y$', 'interpreter', 'latex');
    title('LEG angle y accel');
    legend('L', 'R');

    % FORCES AND TORQUES
    tab3 = uitab(tabgroup, 'Title', 'Forces and Torques');
    axes('Parent', tab3);

    subplot(2, 3, 1)
    hold on; grid on;
    plot(t, lambd_L(:,1), 'LineWidth', 2);
    plot(t, lambd_R(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\lambda_x$', 'Interpreter', 'latex');
    title('Lambda x-axis');

    subplot(2, 3, 2)
    hold on; grid on;
    plot(t, lambd_L(:,2), 'LineWidth', 2);
    plot(t, lambd_R(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\lambda_y$', 'Interpreter', 'latex');
    title('Lambda y-axis');

    subplot(2, 3, 3)
    hold on; grid on;
    plot(t, lambd_L(:,3), 'LineWidth', 2);
    plot(t, lambd_R(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\lambda_z$', 'Interpreter', 'latex');
    title('Lambda z-axis');

    subplot(2, 3, 4)
    hold on; grid on;
    plot(t, tau_L(:,1), 'LineWidth', 2);
    plot(t, tau_R(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\tau_x$', 'Interpreter', 'latex');
    title('Torque x-axis');

    subplot(2, 3, 5)
    hold on; grid on;
    plot(t, tau_L(:,2), 'LineWidth', 2);
    plot(t, tau_R(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\tau_y$', 'Interpreter', 'latex');
    title('Torque y-axis');

    subplot(2, 3, 6)
    hold on; grid on;
    plot(t, tau_L(:,3), 'LineWidth', 2);
    plot(t, tau_R(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\tau_z$', 'Interpreter', 'latex');
    title('Torque z-axis');

    % FRICITON CONES
    tab4 = uitab(tabgroup, 'Title', 'Friction Cones');
    axes('Parent', tab4);

    subplot(1,2,1);
    hold on; grid on; axis equal;
    quiver3(0, 0, 0, 0.4*r, 0, 0, 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0.4*r, 0, 'g', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0, 0.4*r, 'b', 'LineWidth', 2);
    fill3(x_cone, y_cone, z_cone, 'b', 'FaceAlpha', 0.3);
    fill3([0, x_cone], [0, y_cone], [0, z_cone], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot3(lambda_x_L, lambda_y_L, lambda_z_L, 'LineWidth', 2, 'Color', 'b');
    xlabel('$\lambda_x$', 'Interpreter', 'latex');
    ylabel('$\lambda_y$', 'Interpreter', 'latex');
    zlabel('$\lambda_z$', 'Interpreter', 'latex');
    title('Left Leg Friction Cone');
    view(45, 30);

    subplot(1,2,2);
    hold on; grid on; axis equal;
    quiver3(0, 0, 0, 0.4*r, 0, 0, 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0.4*r, 0, 'g', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0, 0.4*r, 'b', 'LineWidth', 2);
    fill3(x_cone, y_cone, z_cone, 'r', 'FaceAlpha', 0.3);
    fill3([0, x_cone], [0, y_cone], [0, z_cone], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot3(lambda_x_R, lambda_y_R, lambda_z_R, 'LineWidth', 2, 'Color', 'r');
    xlabel('$\lambda_x$', 'Interpreter', 'latex');
    ylabel('$\lambda_y$', 'Interpreter', 'latex');
    zlabel('$\lambda_z$', 'Interpreter', 'latex');
    title('Right Leg Friction Cone');
    view(45, 30);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% animate the com trajectory
if animate == 1

    % Create the figure and set its position to the second monitor
    fig = figure('Name', 'Animation', 'WindowState', 'maximized');
    set(gcf,'renderer','painters');

    % plot the z,y,z axis
    hold on; grid on; axis equal;
    quiver3(0, 0, 0, 1, 0, 0, 'r', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 1, 0, 'g', 'LineWidth', 2);
    quiver3(0, 0, 0, 0, 0, 1, 'b', 'LineWidth', 2);
    xlabel('$p_x$ [m]', 'Interpreter', 'latex');
    ylabel('$p_y$ [m]', 'Interpreter', 'latex');
    zlabel('$p_z$ [m]', 'Interpreter', 'latex');
    
    % set the axis limits
    px_min = min([p_com(:,1); x_foot_L(:,1); x_foot_R(:,1)]);
    px_max = max([p_com(:,1); x_foot_L(:,1); x_foot_R(:,1)]);
    py_min = min([p_com(:,2); x_foot_L(:,2); x_foot_R(:,2)]);
    py_max = max([p_com(:,2); x_foot_L(:,2); x_foot_R(:,2)]);
    pz_min = min([p_com(:,3); x_foot_L(:,3); x_foot_R(:,3)]);
    pz_max = max([p_com(:,3); x_foot_L(:,3); x_foot_R(:,3)]);
    xlim([px_min-0.5, px_max+0.5]);
    ylim([py_min-0.5, py_max+0.5]);
    zlim([min(0, pz_min)-0.25, pz_max+0.25]);

    % put a ground plane grid
    x_grid = linspace(px_min-0.5, px_max+0.5, 10);
    y_grid = linspace(py_min-0.5, py_max+0.5, 10);
    [X, Y] = meshgrid(x_grid, y_grid);
    Z = zeros(size(X));
    surf(X, Y, Z, 'FaceAlpha', 0.5, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', [0.7, 0.7, 0.7]);

    % set the view
    if perspective == 'T'      % (T) top view
        azimuth = 0;
        elevation = 90;
    elseif perspective == 'S'  % (S) side view
        azimuth = 0;
        elevation = 0;
    elseif perspective == 'F'  % (F) front view
        azimuth = 90;
        elevation = 0;
    else
        azimuth = 45;
        elevation = 25;
    end
    view(azimuth, elevation);

    % play the animation 
    t  = t * (1/rt);
    for i = 1:replays
        pause(0.5);
        tic;
        ind = 1;
        com_pts = [];
        foot_pts_L = [];
        foot_pts_R = [];
        while true

            % get COM position 
            px_com = p_com(ind,1);
            py_com = p_com(ind,2);
            pz_com = p_com(ind,3);

            % draw the legs
            px_foot_L = x_foot_L(ind,1);
            py_foot_L = x_foot_L(ind,2);
            pz_foot_L = x_foot_L(ind,3);
            px_foot_R = x_foot_R(ind,1);
            py_foot_R = x_foot_R(ind,2);
            pz_foot_R = x_foot_R(ind,3);
            if d_L(ind) == 1
                leg_L = plot3([px_com, px_foot_L], [py_com, py_foot_L], [pz_com, pz_foot_L], 'LineWidth', 3, 'Color', [0 0.4470 0.7410]);
            else
                leg_L = plot3([px_com, px_foot_L], [py_com, py_foot_L], [pz_com, pz_foot_L], 'LineWidth', 3, 'Color', [0 0.4470 0.7410], 'LineStyle', '-.');
            end
            if d_R(ind) == 1
                leg_R = plot3([px_com, px_foot_R], [py_com, py_foot_R], [pz_com, pz_foot_R], 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840]);
            else
                leg_R = plot3([px_com, px_foot_R], [py_com, py_foot_R], [pz_com, pz_foot_R], 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-.');
            end
            ball_foot_L = plot3(px_foot_L, py_foot_L, pz_foot_L, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0 0.4470 0.7410]);
            ball_foot_R = plot3(px_foot_R, py_foot_R, pz_foot_R, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
            
            % draw the mass
            mass = plot3(px_com, py_com, pz_com, 'ko', 'MarkerSize', 35, 'MarkerFaceColor', [0.3, 0.3, 0.3], 'LineWidth', 1.5, 'MarkerEdgeColor', 'k');
            
            %  draw trajectory trail
            if plot_foot == 1
                foot_L = plot3(px_foot_L, py_foot_L, pz_foot_L, 'bo', 'MarkerSize', 1, 'MarkerFaceColor', 'b');
                foot_R = plot3(px_foot_R, py_foot_R, pz_foot_R, 'ro', 'MarkerSize', 1, 'MarkerFaceColor', 'r');
                foot_pts_L = [foot_pts_L; foot_L];
                foot_pts_R = [foot_pts_R; foot_R];
            end
            if plot_com == 1
                pt_pos = plot3(px_com, py_com, pz_com, 'k.', 'MarkerSize', 5);
                com_pts = [com_pts; pt_pos];
            end

            drawnow;
            
            % title
            msg = sprintf('Time: %0.3f [sec]\n px = %0.3f, py = %0.3f, pz = %0.3f\n vx = %0.3f, vy = %0.3f, vz = %0.3f',...
                           t(ind) * rt, p_com(ind,1), p_com(ind,2), p_com(ind,3), v_com(ind,1), v_com(ind,2), v_com(ind,3));
            title(msg);
            
            % wait until the next time step
            while toc< t(ind+1)
                % wait
            end
            
            % increment the index
            if ind+1 >= length(t)
                break;
            else
                ind = ind + 1;
                delete(mass);
                delete(leg_L);
                delete(leg_R);
                delete(ball_foot_L);
                delete(ball_foot_R);
            end
        end

        pause(0.1);

        % clean the plot if still replaying
        if i < replays
            delete(mass);
            delete(leg_L);
            delete(leg_R);
            delete(ball_foot_L);
            delete(ball_foot_R);
            for j = 1:length(com_pts)
                if plot_com == 1
                    delete(com_pts(j));
                end
                if plot_foot == 1
                    delete(foot_pts_L(j));
                    delete(foot_pts_R(j));
                end
            end
        end
    end
end
