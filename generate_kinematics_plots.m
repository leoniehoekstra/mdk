function generate_kinematics_plots()
%GENERATE_KINEMATICS_PLOTS Produce validation figures for direct and
% differential kinematics using preset joint configurations.
%
% This helper runs two automated experiments:
%   1) Joint sweeps to visualise the end-effector pose returned by
%      DIRECT_KINEMATICS for each single-joint motion.
%   2) Jacobian-based differential checks that compare GET_JACOBIAN
%      against finite-difference twists for representative configurations.
%
% The resulting figures are saved under ./figures for convenient reuse in
% reports.

    out_dir = fullfile(pwd, 'figures');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    [unit_twists, H0s] = load_robot_definition();
    run_direct_kinematics_sweeps(unit_twists, H0s, out_dir);
    run_differential_checks(unit_twists, H0s, out_dir);

    fprintf('Validation figures saved to %s\n', out_dir);
end

function run_direct_kinematics_sweeps(unit_twists, H0s, out_dir)
% Sweep each joint independently and record the resulting end-effector
% pose. Generates position and orientation response plots.

    n_joints = numel(unit_twists);
    sweep_range = linspace(-pi/2, pi/2, 75);  % 75 evenly spaced samples

    position_data = zeros(numel(sweep_range), 3, n_joints);
    euler_data = zeros(numel(sweep_range), 3, n_joints);

    base_config = zeros(n_joints, 1);

    for joint_idx = 1:n_joints
        for sample_idx = 1:numel(sweep_range)
            q = base_config;
            q(joint_idx) = sweep_range(sample_idx);

            Hs = direct_kinematics(unit_twists, H0s, q);
            H_ee = Hs{end};

            position_data(sample_idx, :, joint_idx) = H_ee(1:3, 4).';
            euler_data(sample_idx, :, joint_idx) = rotm_to_eul_zyx(H_ee(1:3, 1:3));
        end
    end

    % Position plots
    figure('Name', 'Direct Kinematics - Position Responses', 'NumberTitle', 'off');
    for joint_idx = 1:n_joints
        subplot(2, 3, joint_idx);
        plot(rad2deg(sweep_range), squeeze(position_data(:, :, joint_idx)), 'LineWidth', 1.2);
        xlabel(sprintf('q_%d (deg)', joint_idx));
        ylabel('Position (mm)');
        legend({'x', 'y', 'z'}, 'Location', 'best');
        grid on;
        title(sprintf('Joint %d sweep', joint_idx));
    end
    position_path = fullfile(out_dir, 'direct_kinematics_position.png');
    saveas(gcf, position_path);

    % Orientation plots
    figure('Name', 'Direct Kinematics - Orientation Responses', 'NumberTitle', 'off');
    for joint_idx = 1:n_joints
        subplot(2, 3, joint_idx);
        plot(rad2deg(sweep_range), rad2deg(squeeze(euler_data(:, :, joint_idx))), 'LineWidth', 1.2);
        xlabel(sprintf('q_%d (deg)', joint_idx));
        ylabel('ZYX Euler (deg)');
        legend({'\psi (Z)', '\theta (Y)', '\phi (X)'}, 'Location', 'best');
        grid on;
        title(sprintf('Joint %d sweep', joint_idx));
    end
    orientation_path = fullfile(out_dir, 'direct_kinematics_orientation.png');
    saveas(gcf, orientation_path);
end

function run_differential_checks(unit_twists, H0s, out_dir)
% Evaluate the geometric Jacobian over a grid of joint 1/2 values and
% compare analytic twists with finite-difference estimates for random
% configurations.

    n_joints = numel(unit_twists);

    % Manipulability heatmap using q1/q2 sweeps
    q1 = linspace(-pi/3, pi/3, 40);
    q2 = linspace(-pi/3, pi/3, 40);
    manipulability = zeros(numel(q1), numel(q2));

    for idx1 = 1:numel(q1)
        for idx2 = 1:numel(q2)
            q = zeros(n_joints, 1);
            q(1) = q1(idx1);
            q(2) = q2(idx2);

            J = get_jacobian(unit_twists, H0s, q);
            manipulability(idx1, idx2) = det(J(1:3, 1:3));
        end
    end

    figure('Name', 'Differential Kinematics - Manipulability', 'NumberTitle', 'off');
    surf(rad2deg(q2), rad2deg(q1), manipulability, 'EdgeColor', 'none');
    xlabel('q_2 (deg)'); ylabel('q_1 (deg)'); zlabel('det(J_{pos})');
    title('Position Jacobian determinant');
    colorbar; grid on; view(45, 35);
    manipulability_path = fullfile(out_dir, 'jacobian_manipulability.png');
    saveas(gcf, manipulability_path);

    % Random configuration checks against finite differences
    num_trials = 150;
    step = 1e-4;
    errors = zeros(num_trials, 1);

    for trial = 1:num_trials
        q = (rand(n_joints, 1) - 0.5) * (pi / 2);
        dq = (rand(n_joints, 1) - 0.5);

        J = get_jacobian(unit_twists, H0s, q);
        twist_analytic = J * dq;

        Hs = direct_kinematics(unit_twists, H0s, q);
        H_plus = direct_kinematics(unit_twists, H0s, q + step * dq);

        T0 = Hs{end};
        T1 = H_plus{end};

        R0 = T0(1:3, 1:3); p0 = T0(1:3, 4);
        R1 = T1(1:3, 1:3); p1 = T1(1:3, 4);

        Rdot = (R1 - R0) / step;
        omega_fd = vee(Rdot * R0');
        v_fd = (p1 - p0) / step;

        twist_fd = [omega_fd; v_fd];
        errors(trial) = norm(twist_analytic - twist_fd);
    end

    figure('Name', 'Differential Kinematics - Velocity Agreement', 'NumberTitle', 'off');
    histogram(errors, 24);
    xlabel('||J dq - twist_{FD}||'); ylabel('Frequency');
    title('Jacobian vs. finite-difference twist');
    grid on;
    velocity_path = fullfile(out_dir, 'jacobian_velocity_error.png');
    saveas(gcf, velocity_path);
end

function [unit_twists, H0s] = load_robot_definition()
% Replicates the parameter set used by robot.m so the validation runs with
% the same serial-chain description.

    a = 100; b = 800; c = 800; d = 200; e = 100; f = d / 2;

    unit_twists = {
        [0; 0; 1; 0; 0; 0];
        [0; -1; 0; 0; 0; 0];
        [0; 1; 0; 0; 0; 0];
        [0; -1; 0; 0; 0; 0];
        [0; 0; -1; 0; 0; 0];
        [0; -1; 0; 0; 0; 0] };

    H0s = {
        [1 0 0 0; 0 1 0 -d/2; 0 0 1 a;   0 0 0 1];
        [1 0 0 0; 0 1 0 0;    0 0 1 b;   0 0 0 1];
        [1 0 0 0; 0 1 0 0;    0 0 1 c;   0 0 0 1];
        [1 0 0 0; 0 1 0 -d/2; 0 0 1 e/2; 0 0 0 1];
        [1 0 0 0; 0 1 0 d/2;  0 0 1 e/2; 0 0 0 1];
        [1 0 0 0; 0 1 0 f;    0 0 1 0;   0 0 0 1] };
end

function eul = rotm_to_eul_zyx(R)
% Convert rotation matrix to ZYX Euler angles (psi, theta, phi).
    if abs(R(3,1)) < 1
        theta = atan2(-R(3,1), hypot(R(3,2), R(3,3)));
        psi = atan2(R(2,1), R(1,1));
        phi = atan2(R(3,2), R(3,3));
    else
        % Gimbal lock fallback
        theta = atan2(-R(3,1), 0);
        psi = atan2(-R(1,2), R(2,2));
        phi = 0;
    end
    eul = [psi, theta, phi];
end

function v = vee(S)
% Extract the vector corresponding to a 3x3 skew-symmetric matrix.
    v = [S(3,2); S(1,3); S(2,1)];
end

function deg = rad2deg(rad)
% Local rad2deg to avoid toolbox dependencies.
    deg = (180/pi) * rad;
end

function rad = deg2rad(deg)
% Local deg2rad for completeness (unused but kept for symmetry).
    rad = (pi/180) * deg;
end
