%% Rotary Inverted Pendulum: IOFBL + State-Feedback Control
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. SYMBOLIC MODEL DERIVATION
% -------------------------------------------------------------------------
syms J1 m2 L1 l2 g C1 J2 real
syms x1 x2 x3 x4 real   % x1=theta1, x2=theta1_dot, x3=theta2, x4=theta2_dot
syms x2_dot x4_dot tau real

% Nonlinear equations of motion
eq1 = (J1 + m2*L1^2)*x2_dot + ...
      m2*L1*l2*sin(x3)*x4^2 - ...
      m2*L1*l2*cos(x3)*x4_dot == tau - C1*x2;

eq2 = (J2 + m2*l2^2)*x4_dot + ...
      m2*L1*l2*cos(x3)*x2_dot - ...
      m2*g*l2*sin(x3) == 0;

% Solve for accelerations
solution = solve([eq1, eq2], [x2_dot, x4_dot]);
x2_dot_expr = simplify(solution.x2_dot);
x4_dot_expr = simplify(solution.x4_dot);

%% ------------------------------------------------------------------------
% 2. INPUT-OUTPUT FEEDBACK LINEARIZATION TERMS
% -------------------------------------------------------------------------
y = x3;
y_ddot_expr = x4_dot_expr;
b_x = jacobian(y_ddot_expr, tau);
a_x = simplify(subs(y_ddot_expr, tau, 0));

disp('ÿ = a(x) + b(x)*τ');
disp('a(x) = '); pretty(a_x);
disp('b(x) = '); pretty(b_x);

control_law = simplify((sym('v') - a_x)/b_x);
disp('Control Law: τ = (v - a(x)) / b(x)');

%% ------------------------------------------------------------------------
% 3. ZERO-DYNAMICS ANALYSIS
% -------------------------------------------------------------------------
constraint_vars = {x3, x4};
constraint_vals = {0, 0};
a_x_zero = subs(a_x, constraint_vars, constraint_vals);
b_x_zero = subs(b_x, constraint_vars, constraint_vals);
tau_star = simplify(-a_x_zero / b_x_zero);

x2_dot_constrained = subs(x2_dot_expr, constraint_vars, constraint_vals);
omega_2_dot = simplify(subs(x2_dot_constrained, tau, tau_star));
disp('Zero Dynamics:');
pretty(sym('omega_1_dot') == sym('omega_2'));
pretty(sym('omega_2_dot') == omega_2_dot);

%% ------------------------------------------------------------------------
% 4. LINEAR STATE-FEEDBACK DESIGN (POLE PLACEMENT)
% -------------------------------------------------------------------------
% The linearised system is: ÿ = v = -K1*y - K2*ẏ
% Choose K1, K2 via desired poles (ζ, ω_n)
zeta = 0.8;     % damping ratio
ts = 0.2;       % desired settling time [s]

% --- Determine Natural frequency ---
omega_n = ( 1/(ts*sqrt(1-zeta^2)) )*(pi - atan((sqrt(1-zeta^2))/zeta));

% --- Calculate Gains ---
K1 = omega_n^2;             % k1 = ω_n^2
K2 = 2*zeta*omega_n;        % k2 = 2ζω_n


% --- Calculate Percentage Overshoot ---
PO = exp(-zeta*pi / sqrt(1 - zeta^2)) * 100;

fprintf('\nChosen damping ζ = %.2f, Settling time t_s = %.2f s\n', zeta, ts);
fprintf('Natural Frequency (exact) : ω_n = %.3f rad/s\n', omega_n);
fprintf('Percentage Overshoot      : %.2f%%\n', PO);
fprintf('Controller Gains: K1 = %.2f, K2 = %.2f\n', K1, K2);

%% ------------------------------------------------------------------------
% 5. NUMERICAL SIMULATION
% -------------------------------------------------------------------------
% Physical parameters
param.J1 = 0.01;  param.m2 = 0.032;  param.L1 = 0.19;
param.l2 = 0.11/2; param.g = 9.81;  param.C1 = 0.02;  param.J2 = 0.00014;

% Convert symbolic to numeric functions
vars_plant = [x1; x2; x3; x4; tau; J1; m2; L1; l2; g; C1; J2];
state_model_f = matlabFunction([x2_dot_expr; x4_dot_expr], 'Vars', {vars_plant});
a_f = matlabFunction(a_x, 'Vars', {[x1; x2; x3; x4; J1; m2; L1; l2; g; C1; J2]});
b_f = matlabFunction(b_x, 'Vars', {[x1; x2; x3; x4; J1; m2; L1; l2; g; C1; J2]});

% Initial conditions
x0 = [0; 0; 0.03; 0]; 
tspan = [0 10];
y_d = 0;

% ODE solver
[t, x_hist] = ode45(@(t,x) rip_ode(t,x,param,state_model_f,a_f,b_f,K1,K2,y_d), tspan, x0);

%% ------------------------------------------------------------------------
% 6. PLOTTING
% -------------------------------------------------------------------------
figure('Name','State-Feedback Controller Performance');

subplot(2,2,1);
plot(t,x_hist(:,3),'b','LineWidth',2); hold on;
plot(t,y_d*ones(size(t)),'r--','LineWidth',1.5);
xlabel('Time (s)','FontSize',12);
ylabel('\theta_2 (rad)','FontSize',12);
title('Pendulum Angle Response','FontSize',15);
legend('\theta_2','Desired','Location','Best');
grid on; set(gca,'FontSize',12);

subplot(2,2,2);
plot(t,x_hist(:,1),'m','LineWidth',2);
xlabel('Time (s)','FontSize',12);
ylabel('\theta_1 (rad)','FontSize',12);
title('Arm Angle','FontSize',15);
grid on; set(gca,'FontSize',12);

subplot(2,2,3);
plot(t,x_hist(:,4),'g','LineWidth',2);
xlabel('Time (s)','FontSize',12);
ylabel('\thetȧ_2 (rad/s)','FontSize',12);
title('Pendulum Angular Velocity','FontSize',15);
grid on; set(gca,'FontSize',12);

subplot(2,2,4);
plot(t,x_hist(:,2),'c','LineWidth',2);
xlabel('Time (s)','FontSize',12);
ylabel('\thetȧ_1 (rad/s)','FontSize',12);
title('Arm Angular Velocity','FontSize',15);
grid on; set(gca,'FontSize',12);




%% ------------------------------------------------------------------------
% 7. ODE FUNCTION (STATE-FEEDBACK)
% -------------------------------------------------------------------------
function dx = rip_ode(~,x,param,state_model_f,a_f,b_f,K1,K2,y_d)
    p = struct2array(param)';
    a_val = a_f([x; p]);
    b_val = b_f([x; p]);

    % State-feedback on (y, ydot) = (x3, x4)
    y = x(3); ydot = x(4);
    v = -K1*(y - y_d) - K2*ydot;

    tau = (v - a_val)/b_val;
    accel = state_model_f([x; tau; p]);
    dx = [x(2); accel(1); x(4); accel(2)];
end
