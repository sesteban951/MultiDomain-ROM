%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SOME SIM RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Load data
data_folder = "../data/";

t = load(data_folder + 'time.csv');
x_sys = load(data_folder + 'state_sys.csv');
x_leg = load(data_folder + 'state_leg.csv');
x_foot = load(data_folder + 'state_foot.csv');
u = load(data_folder + 'input.csv');
lambd = load(data_folder + 'lambda.csv');
tau = load(data_folder + 'tau.csv');
d = load(data_folder + 'domain.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unapack varaibles from the YAML config
config_file = '../config/config.yaml';
config = yaml.loadFile(config_file);

% extract some variables
r_des = config.REFERENCE.r_des;
theta_des = config.REFERENCE.theta_des;
interp = config.SYS_PARAMS.interp;

% trajecotry variables
px_traj = config.REFERENCE.pz_des;
pz_traj = config.REFERENCE.pz_des;
t_traj = config.REFERENCE.time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% segment the time
t_interval = [t(1) t(end)];
% t_interval = [0 2.8];

% plotting / animation
animate = 1;   % animation = 1; plot states = 0
rt = 1.0;      % realtime rate
replays = 1;   % how many times to replay the animation
plot_com = 1;  % plot the foot trajectory
plot_foot = 1; % plot the foot trajectory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% system state
p_com = x_sys(:,1:2);
v_com = x_sys(:,3:4);
x_leg_commands_L = x_sys(:,5:6);
x_leg_commands_R = x_sys(:,7:8);

% leg states
x_leg_L = x_leg(:,1:4);
x_leg_R = x_leg(:,5:8);

% foot states
x_foot_L = x_foot(:,1:4);
x_foot_R = x_foot(:,5:8);

% inputs
u_L = u(:,1:2);
u_R = u(:,3:4);

% lambda leg forces
lambd_L = lambd(:,1:2);
lambd_R = lambd(:,3:4);
lambd_L_norm = zeros(length(t), 1);
lambd_R_norm = zeros(length(t), 1);
for i = 1:length(t)
    lambd_L_norm(i) = norm(lambd_L(i,:));
    lambd_R_norm(i) = norm(lambd_R(i,:));
end

% ankle torques
tau_L = tau(:,1);
tau_R = tau(:,2);

% domain
d_L = d(:,1);
d_R = d(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if animate == 0
    % plot all states
    figure('Name', 'COM States', 'WindowState', 'maximized');

    % COM STATES
    subplot(3,6,1);
    hold on; grid on;
    plot(t, p_com(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_x$ [m]', 'Interpreter', 'latex');
    title('COM x-pos');

    subplot(3,6,2);
    hold on; grid on;
    plot(t, p_com(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_z$ [m]', 'Interpreter', 'latex');
    title('COM z-pos');

    subplot(3,6,7); 
    hold on; grid on;
    plot(t, v_com(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_x$ [m/s]', 'Interpreter', 'latex');
    title('COM x-vel');

    subplot(3,6,8);
    hold on; grid on;
    plot(t, v_com(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_z$ [m/s]', 'Interpreter', 'latex');
    title('COM z-vel');

    % LEG STATES
    subplot(3,6,3);
    hold on; grid on;
    plot(t, x_leg_L(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    plot(t, x_leg_commands_L(:,1), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,1), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    plot(t, x_leg_commands_R(:,1), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    xlabel('Time [sec]');
    ylabel('$r$ [m]', 'Interpreter', 'latex');
    title('LEG Length');
    % legend('$r_L$', '$\hat{r}_L$', '$r_R$', '$\hat{r}_R$', 'Interpreter', 'latex');

    subplot(3,6,4);
    hold on; grid on;
    plot(t, x_leg_L(:,2), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    plot(t, x_leg_commands_L(:,2), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    plot(t, x_leg_commands_R(:,2), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    xlabel('Time [sec]');
    ylabel('$\theta$ [rad]', 'Interpreter', 'latex');
    title('LEG Angle');
    % legend('$\theta_L$', '$\hat{\theta}_L$', '$\theta_R$', '$\hat{\theta}_R$', 'Interpreter', 'latex');

    subplot(3,6,9);
    hold on; grid on;
    plot(t, x_leg_L(:,3), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    plot(t, u_L(:,1), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,3), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    plot(t, u_R(:,1), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    xlabel('Time [sec]');
    ylabel('$\dot{r}$ [m/s]', 'Interpreter', 'latex');
    title('LEG Length Vel');
    % legend('$\dot{r}_L$', '$\dot{\hat{r}}_L$', '$\dot{r}_R$', '$\dot{\hat{r}}_R$', 'Interpreter', 'latex');

    subplot(3,6,10);
    hold on; grid on;
    plot(t, x_leg_L(:,4), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    plot(t, u_L(:,2), 'LineWidth', 1.0, 'Color', [0.8500 0.3250 0.0980]);
    plot(t, x_leg_R(:,4), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    plot(t, u_R(:,2), 'LineWidth', 1.0, 'Color', [0.4660 0.6740 0.1880]);
    xlabel('Time [sec]');
    ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');
    title('LEG Angle Vel');
    % legend('$\dot{\theta_L}$', '$\dot{\hat{\theta}}_L$', '$\dot{\theta}_R$', '$\dot{\hat{\theta}}_R$', 'Interpreter', 'latex');

    % FOOT STATES
    subplot(3,6,5);
    hold on; grid on;
    plot(t, x_foot_L(:,1), 'LineWidth', 2);
    plot(t, x_foot_R(:,1), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_{foot,x}$ [m]', 'Interpreter', 'latex');
    title('FOOT x-pos');
    legend('L', 'R');

    subplot(3,6,6);
    hold on; grid on;
    plot(t, x_foot_L(:,2), 'LineWidth', 2);
    plot(t, x_foot_R(:,2), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$p_{foot,z}$ [m]', 'Interpreter', 'latex');
    title('FOOT z-pos');
    legend('L', 'R');

    subplot(3,6,11);
    hold on; grid on;
    plot(t, x_foot_L(:,3), 'LineWidth', 2);
    plot(t, x_foot_R(:,3), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_{foot,x}$ [m/s]', 'Interpreter', 'latex');
    title('FOOT x-vel');
    legend('L', 'R');

    subplot(3,6,12);
    hold on; grid on;
    plot(t, x_foot_L(:,4), 'LineWidth', 2);
    plot(t, x_foot_R(:,4), 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$v_{foot,z}$ [m/s]', 'Interpreter', 'latex');
    title('FOOT z-vel');
    legend('L', 'R');

    % INPUT
    subplot(3,6,13); 
    hold on; grid on;
    if interp == 'Z'
        stairs(t, u_L(:,1), 'LineWidth', 2);
        stairs(t, u_R(:,1), 'LineWidth', 2);
    elseif interp == 'L'
        plot(t, u_L(:,1), 'LineWidth', 2);
        plot(t, u_R(:,1), 'LineWidth', 2);
    end
    xlabel('Time [sec]');
    ylabel('$\hat{\dot{l_0}}$', 'interpreter', 'latex');
    title('leg vel input');
    legend('L', 'R');

    subplot(3,6,14); 
    hold on; grid on;
    if interp == 'Z'
        stairs(t, u_L(:,2), 'LineWidth', 2);
        stairs(t, u_R(:,2), 'LineWidth', 2);
    elseif interp == 'L'
        plot(t, u_L(:,2), 'LineWidth', 2);
        plot(t, u_R(:,2), 'LineWidth', 2);
    end
    xlabel('Time [sec]');
    ylabel('$\hat{\dot{\theta}}$', 'interpreter', 'latex');
    title('angle vel input');
    legend('L', 'R');

    % LAMBDA LEG FORCES
    subplot(3,6,[15:16]);
    hold on; grid on;
    plot(t, lambd_L_norm, 'LineWidth', 2);
    plot(t, lambd_R_norm, 'LineWidth', 2);
    xlabel('Time [sec]');
    ylabel('$\|\lambda\|$', 'interpreter', 'latex');
    title('Leg Force');
    legend('L', 'R');

    % DOMAIN
    subplot(3,6,[17:18]);
    hold on; grid on;
    stairs(t, d_L, 'LineWidth', 2);
    stairs(t, d_R, 'LineWidth', 1.5);
    xlabel('Time [sec]');
    title('Domain');
    ylim([-0.5, 1.5]);
    yticks([0, 1]);
    yticklabels({'F', 'G'});
    legend('L', 'R');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% animate the com trajectory
if animate == 1

    figure('Name', 'Animation');
    hold on;

    xline(0);
    yline(0);
    xlabel('$p_x$ [m]', 'Interpreter', 'latex');
    ylabel('$p_z$ [m]', 'Interpreter', 'latex');
    grid on; axis equal;
    px_min = min([p_com(:,1); x_foot_L(:,1); x_foot_R(:,1)]);
    px_max = max([p_com(:,1); x_foot_L(:,1); x_foot_R(:,1)]);
    pz_min = min([p_com(:,2); x_foot_L(:,2); x_foot_R(:,2)]);
    pz_max = max([p_com(:,2); x_foot_L(:,2); x_foot_R(:,2)]);
    xlim([px_min-0.5, px_max+0.5]);
    ylim([min(0, pz_min)-0.25, pz_max+0.25]);

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
            px = p_com(ind,1);
            pz = p_com(ind,2);

            % draw the legs
            px_foot_L = x_foot_L(ind,1);
            pz_foot_L = x_foot_L(ind,2);
            px_foot_R = x_foot_R(ind,1);
            pz_foot_R = x_foot_R(ind,2);
            if d_L(ind) == 1
                leg_L = plot([px, px_foot_L], [pz, pz_foot_L], 'LineWidth', 3, 'Color', [0 0.4470 0.7410]);
            else
                leg_L = plot([px, px_foot_L], [pz, pz_foot_L], 'LineWidth', 3, 'Color', [0 0.4470 0.7410], 'LineStyle', '-.');
            end
            if d_R(ind) == 1
                leg_R = plot([px, px_foot_R], [pz, pz_foot_R], 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840]);
            else
                leg_R = plot([px, px_foot_R], [pz, pz_foot_R], 'LineWidth', 3, 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-.');
            end
            ball_foot_L = plot(px_foot_L, pz_foot_L, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0 0.4470 0.7410]);
            ball_foot_R = plot(px_foot_R, pz_foot_R, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
            
            % draw the mass
            mass = plot(px, pz, 'ko', 'MarkerSize', 35, 'MarkerFaceColor', [0.3, 0.3, 0.3], 'LineWidth', 1.5, 'MarkerEdgeColor', 'k');
            
            %  draw trajectory trail
            if plot_foot == 1
                foot_L = plot(px_foot_L, pz_foot_L, 'bo', 'MarkerSize', 1, 'MarkerFaceColor', 'b');
                foot_R = plot(px_foot_R, pz_foot_R, 'ro', 'MarkerSize', 1, 'MarkerFaceColor', 'r');
                foot_pts_L = [foot_pts_L; foot_L];
                foot_pts_R = [foot_pts_R; foot_R];
            end
            if plot_com == 1
                pt_pos = plot(px, pz, 'k.', 'MarkerSize', 5);
                com_pts = [com_pts; pt_pos];
            end

            drawnow;
            
            % title
            msg = sprintf('Time: %0.3f [sec]\n vx = %0.3f, px = %0.3f\n vz = %0.3f, pz = %0.3f',...
                         t(ind) * rt, v_com(ind,1), p_com(ind,1), v_com(ind,2), p_com(ind,2));
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
