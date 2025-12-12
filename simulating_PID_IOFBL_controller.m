%% ========================================================================
% ROTARY INVERTED PENDULUM (RIP)
% Input–Output Feedback Linearisation (IOFBL) with PID Controller
% -------------------------------------------------------------------------
% This script:
%   1. Derives the nonlinear equations of motion symbolically.
%   2. Computes the IOFBL terms a(x) and b(x).
%   3. Builds a PID-based control law.
%   4. Simulates the nonlinear closed-loop dynamics using ode45.
%   5. Plots pendulum and arm position/velocity responses.
% ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Define system parameters and symbolic variables
% ------------------------------------------------------------------------
syms J1 m2 L1 l2 g C1 J2 real          % physical constants
syms x1 x2 x3 x4 real                  % states: θ1, θ̇1, θ2, θ̇2
syms x2_dot x4_dot tau real            % accelerations and control torque

%% ------------------------------------------------------------------------
% 2. Equations of motion from Lagrange formulation
% ------------------------------------------------------------------------
eq1 = (J1 + m2*L1^2)*x2_dot + ...
       m2*L1*l2*sin(x3)*x4^2 - ...
       m2*L1*l2*cos(x3)*x4_dot == tau - C1*x2;

eq2 = (J2 + m2*l2^2)*x4_dot + ...
       m2*L1*l2*cos(x3)*x2_dot - ...
       m2*g*l2*sin(x3) == 0;

% Solve for accelerations θ̈1 and θ̈2
solution   = solve([eq1, eq2], [x2_dot, x4_dot]);
x2_dot_expr = solution.x2_dot;
x4_dot_expr = solution.x4_dot;

%% ------------------------------------------------------------------------
% 3. Derive IOFBL components a(x) and b(x)
% ------------------------------------------------------------------------
y = x3;                                      % output: pendulum angle
y_ddot_expr = x4_dot_expr;                   % second derivative of output
b_x = jacobian(y_ddot_expr, tau);            % input gain term
a_x = subs(y_ddot_expr, tau, 0);             % nonlinear dynamics without input

% IOFBL control law: tau = (v - a(x)) / b(x)
syms v real
control_law = (v - a_x) / b_x;

%% ------------------------------------------------------------------------
% 4. Zero-dynamics analysis (internal hidden dynamics)
% ------------------------------------------------------------------------
constraint_vars = {x3, x4};
constraint_vals = {0, 0};
a_x_zero = subs(a_x, constraint_vars, constraint_vals);
b_x_zero = subs(b_x, constraint_vars, constraint_vals);
tau_star = simplify(-a_x_zero / b_x_zero);

syms omega_1 omega_2 real
omega_1_dot = omega_2;
x2_dot_constrained = subs(x2_dot_expr, constraint_vars, constraint_vals);
omega_2_dot_expr = subs(x2_dot_constrained, tau, tau_star);
omega_2_dot = simplify(subs(omega_2_dot_expr, x2, omega_2));
omega_2_dot
%% ------------------------------------------------------------------------
% 5. PID control law design (virtual input v)
% ------------------------------------------------------------------------
syms Kp Kd Ki y_d e_int real
e     = y_d - y;
e_dot = -x4;
v_pid_law = Kp*e + Ki*e_int + Kd*e_dot;
final_control_law = subs(control_law, v, v_pid_law);

%% ------------------------------------------------------------------------
% 6. Numerical simulation setup
% ------------------------------------------------------------------------
disp('--- Running Numerical Simulation ---');

% Physical parameters (numeric values)
param_values.J1 = 0.01;
param_values.m2 = 0.032;
param_values.L1 = 0.19;
param_values.l2 = 0.11/2;
param_values.g  = 9.81;
param_values.C1 = 0.02;
param_values.J2 = 0.00014;

% Convert symbolic expressions to numeric functions for speed
vars_plant   = [x1; x2; x3; x4; tau; J1; m2; L1; l2; g; C1; J2];
vars_control = [x1; x2; x3; x4; J1; m2; L1; l2; g; C1; J2];
state_model_func = matlabFunction([x2_dot_expr; x4_dot_expr], 'Vars', {vars_plant});
a_func = matlabFunction(a_x, 'Vars', {vars_control});
b_func = matlabFunction(b_x, 'Vars', {vars_control});

% Simulation parameters
tspan   = [0 10];
x0      = [0; 0; 0.03; 0];          % initial conditions (0.03-rad disturbance)
y_d_val = 0;                        % desired pendulum angle
Kp_val  = 120;  Ki_val = 0.02;  Kd_val = 25;  % PID gains

% Integrate system dynamics
[t, state_history] = ode45(@(t, s) pendulum_ode(t, s, param_values, ...
                            state_model_func, a_func, b_func, ...
                            y_d_val, Kp_val, Kd_val, Ki_val), ...
                            tspan, [x0; 0]);  % includes integral error state

%% ------------------------------------------------------------------------
% 7. Plot results: pendulum and arm positions & velocities
% ------------------------------------------------------------------------
figure('Name','IOFBL with PID Controller Performance');

subplot(2,2,1);
plot(t, state_history(:,3), 'b', 'LineWidth', 2); hold on;
plot(t, y_d_val*ones(size(t)), 'r--', 'LineWidth', 1.5);
grid on;
title('Pendulum Angle (\theta_2)','FontSize',15);
xlabel('Time (s)','FontSize',12); ylabel('Angle (rad)','FontSize',12);
legend('\theta_2 Actual','\theta_2 Desired','Location','SouthEast','FontSize',12);
set(gca,'FontSize', 12);

subplot(2,2,3);
plot(t, state_history(:,4), 'g', 'LineWidth', 2);
grid on;
title('Pendulum Angular Velocity','FontSize',15);
xlabel('Time (s)','FontSize',12); ylabel('Angular Velocity (rad/s)','FontSize',12);
set(gca,'FontSize',12);

subplot(2,2,2);
plot(t, state_history(:,1), 'm', 'LineWidth', 2);
grid on;
title('Arm Angle (\theta_1)','FontSize',15);
xlabel('Time (s)','FontSize',12); ylabel('Angle (rad)','FontSize',12);
set(gca,'FontSize',12);

subplot(2,2,4);
plot(t, state_history(:,2), 'c', 'LineWidth', 2);
grid on;
title('Arm Angular Velocity','FontSize',15);
xlabel('Time (s)','FontSize',12); ylabel('Angular Velocity (rad/s)','FontSize',12);
set(gca,'FontSize',12);

%% ------------------------------------------------------------------------
% 8. Dynamics function for ode45
% ------------------------------------------------------------------------
function d_state = pendulum_ode(~, state, params, state_model_f, a_f, b_f, y_d, Kp, Kd, Ki)
    % Extract states
    x_vec = state(1:4);    % [theta1, theta1_dot, theta2, theta2_dot]
    e_int = state(5);      % integral of error
    p     = struct2array(params)';

    % Compute feedback linearising control
    e = y_d - x_vec(3);
    e_dot = -x_vec(4);
    v = Kp*e + Ki*e_int + Kd*e_dot;
    a_val = a_f([x_vec; p]);
    b_val = b_f([x_vec; p]);
    tau = (v - a_val) / b_val;

    % System dynamics
    x_dot_vals = state_model_f([x_vec; tau; p]);
    d_e_int = e;   % integral error derivative

    % Return state derivatives
    d_state = [x_vec(2);
               x_dot_vals(1);
               x_vec(4);
               x_dot_vals(2);
               d_e_int];
end
