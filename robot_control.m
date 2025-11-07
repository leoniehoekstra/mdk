function robot_control()

%% Robot parameters
a = 100; b = 800; c = 800; d = 200; e = 100; f = d/2;
q = zeros(6, 1);  % Initialize joint variables

% >>>> TO MODIFY >>>>>>>>
% Unit twists of the joints T_i^{(i-1), (i-1)}
T1 = [0; 0; 1; 0; 0; 0];
T2 = [0; -1; 0; 0; 0; 0];
T3 = [0; 1; 0; 0; 0; 0];
T4 = [0; -1; 0; 0; 0; 0];
T5 = [0; 0; -1; 0; 0; 0];
T6 = [0; -1; 0; 0; 0; 0];

unit_twists = {T1, T2, T3, T4, T5, T6};

%% Test Jacobian and Tbb_ee
format shortG
T = 0.2;

% All joints have velocity 0
dq = zeros(6, 1);

% Joint 1 has velocity 1
dq(1) = pi/2;



% Simulate 1 second with timestep T
for t = 0:T:1

    % Move joints by dq*t
    q = q + T*dq;

    % Configurations of the joints H_i^{i-1}(0)
    H0_1 = [cos(q(1)), -sin(q(1)), 0, (d/2)*sin(q(1)); sin(q(1)), cos(q(1)), 0, (-d/2)*cos(q(1)); 0, 0, 1, a; 0, 0, 0, 1];
    H1_2 = [cos(q(2)), 0, -sin(q(2)), -b*sin(q(2)); 0, 1, 0, 0; sin(q(2)), 0, cos(q(2)), b*cos(q(2)); 0, 0, 0, 1];
    H2_3 = [cos(q(3)), 0, sin(q(3)), c*sin(q(3)); 0, 1, 0, 0; -sin(q(3)), 0, cos(q(3)), c*cos(q(3)); 0, 0, 0, 1];
    H3_4 = [cos(q(4)), 0, -sin(q(4)), (-e/2)*sin(q(4)); 0, 1, 0, (-d/2); sin(q(4)), 0, cos(q(4)), (e/2)*cos(q(4)); 0, 0, 0, 1];
    H4_5 = [cos(q(5)), -sin(q(5)), 0, (d/2)*sin(q(5)); sin(q(5)), cos(q(5)), 0, (d/2)*cos(q(5)); 0, 0, 1, (e/2); 0, 0, 0, 1];
    H5_6 = [cos(q(6)), 0, sin(q(6)), 0; 0, 1, 0, f; -sin(q(6)), 0, cos(q(6)), 0; 0, 0, 0, 1];
    H0s = {H0_1, H1_2, H2_3, H3_4, H4_5, H5_6};

    fprintf('Joint 1 = %.2f pi\n', q(1)/pi);

    disp("Jacobian:");
    disp(round(get_jacobian(unit_twists, H0s, q), 2));

    disp("Twist:");
    Tbb_ee = get_jacobian(unit_twists, H0s, q)*dq;
    disp(Tbb_ee);
end

% <<<<<<<<<<<<<<<<<<<<

%% Control Part (not needed for us)

%%
%{
% Control Loop Parameters
T = 0.02;       % Time step
t_end = 20;     % Total simulation time
N = round(t_end / T);  % Number of iterations
t = 0;          % Initialize time vector

% Create figure and axes for the control loop
fig2 = figure('Name', 'Control 6-DOF Robot Manipulator', 'NumberTitle', 'off');
ax2 = axes('Parent', fig2, 'Projection', 'perspective'); % Set projection to perspective
hold(ax2, 'on');
grid(ax2, 'on');
axis(ax2, 'equal');
xlabel(ax2, 'X');
ylabel(ax2, 'Y');
zlabel(ax2, 'Z');
title(ax2, 'Robot Trajectory Tracking');
view(ax2, [1, 1, 1]);
xlim(ax2, [-500, 500]);
ylim(ax2, [-500, 500]);
zlim(ax2, [-100, 1000]);
axis(ax2, 'manual');  % Fix axis limits

% Before the control loop
trajectory = [];  % Initialize trajectory array

% Main control loop

for i = 1:N
    % Current time
    t = (i - 1) * T;

    % Generate setpoint trajectory
    x = 600;
    y = 200 * cos(t);
    z = 400 + 200 * sin(2*t);
    setpoint = [x; y; z];

    % Compute desired velocities (if applicable)
    dx = 0;
    dy = -200 * sin(t);
    dz = 400 * cos(2*t);
    dsetpoint = [dx; dy; dz];

    % Store the trajectory
    trajectory = [trajectory, setpoint];

    % Calculate desired joint velocities
    %dq = calculate_dq(q, setpoint, dsetpoint, H0s, unit_twists);

    % Update joint variables
    q = q + T * dq;

    % Compute the transformations for the current q
    Hs = direct_kinematics(unit_twists, H0s, q);

    % Plot the robot
    plot_robot(ax2, Hs);

    % Plot the desired end-effector position
    %plot3(ax2, setpoint(1), setpoint(2), setpoint(3), 'kx', 'MarkerSize', 10, 'LineWidth', 2);

    % Plot the trajectory
    %plot3(ax2, trajectory(1, :), trajectory(2, :), trajectory(3, :), 'b--');

    % Update the title with the current time
    title(ax2, sprintf('Robot Trajectory Tracking (t = %.2f s)', t));

    drawnow;
    pause(T);  % Control the loop timing


end

hold(ax2, 'off');
%}
end