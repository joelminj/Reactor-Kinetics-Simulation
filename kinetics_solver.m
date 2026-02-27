% RP Activity 5
% [Author]: [Joel Minj]

clear; clc; close all;

% NEUTRONIC PARAMETERS
td = 1e-4; % Prompt neutron lifetime (s)
beta_i = [0.00025, 0.00127, 0.00125, 0.00272, 0.00074, 0.00027]; 
lambda_i = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]; 
beta = sum(beta_i);

% DECAY HEAT PARAMETERS (3-Group Model)
% gamma_j: Fraction of total power, lambda_hj: Decay constants (s^-1)
gamma_j = [0.035, 0.020, 0.015]; % Total decay heat = 7% of nominal power
lambda_hj = [0.1, 0.01, 0.001];  % Short, medium, and long-lived groups
prompt_fraction = 1 - sum(gamma_j);

%THERMAL HYDRAULIC PARAMETERS
P0 = 3000e6;           % Nominal Thermal Power (W)
T_in = 290;            % Coolant inlet temperature (C)
M_F_Cp_F = 3e7;        % Fuel heat capacity (J/K)
M_C_Cp_C = 1e8;        % Coolant heat capacity (J/K)
UA = 1e7;              % Heat transfer coeff Fuel->Coolant (W/K)
m_dot_Cp_C = 1e8;      % Coolant mass flow heat capacity rate (W/K)

alpha_F = -3e-5;       % Fuel Doppler coefficient (-3 pcm/C)
alpha_C = -20e-5;      % Coolant moderator coefficient (-20 pcm/C)

% STEADY STATE INITIALIZATION
f_c = 0.025;
T_C0 = T_in + P0 / (2 * m_dot_Cp_C);  
T_F0 = T_C0 + ((1 - f_c) * P0) / UA;                

% Time parameters
t_span = [0 200]; % Simulate for 200 seconds to watch the decay heat

% INITIAL CONDITIONS
phi_0 = 1;
C_0 = (beta_i' ./ (td * lambda_i')) * phi_0; 
P_d0 = (gamma_j .* P0 * phi_0)'; % Initial decay heat for the 3 groups

% State Vector: [Flux(1); Precursors(6); T_F(1); T_C(1); DecayHeat(3)] -> 12x1
Y_0 = [phi_0; C_0; T_F0; T_C0; P_d0]; 

% ODE SYSTEM DEFINITION
kinetics_ode = @(t, Y) reactor_dynamics(t, Y, beta, beta_i, lambda_i, td, ...
                                        P0, T_in, M_F_Cp_F, M_C_Cp_C, UA, m_dot_Cp_C, ...
                                        alpha_F, alpha_C, T_F0, T_C0, ...
                                        gamma_j, lambda_hj, prompt_fraction);

% Solve using ode15s
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_out, Y_out] = ode15s(kinetics_ode, t_span, Y_0, options);

% EXTRACT RESULTS
phi_res = Y_out(:, 1);
T_F_res = Y_out(:, 8);
T_C_res = Y_out(:, 9);
P_decay_total = sum(Y_out(:, 10:12), 2); % Sum of the 3 decay heat groups

% Calculate Total Thermal Power for plotting
P_total_res = P0 * phi_res * prompt_fraction + P_decay_total;

% Reconstruct rho_ext and rho_net for plotting
rho_ext_res = zeros(size(t_out));
rho_net_res = zeros(size(t_out));
for k = 1:length(t_out)
    rho_ext_res(k) = get_rho_ext(t_out(k));
    rho_net_res(k) = rho_ext_res(k) + alpha_F*(T_F_res(k)-T_F0) + alpha_C*(T_C_res(k)-T_C0);
end

% PLOTS
fig = figure('Position', [100, 50, 1200, 1000], 'Color', 'white');
set(fig, 'DefaultLineLineWidth', 2);

% Plot 1: Reactivity Profile
ax1 = subplot(4, 1, 1);
plot(t_out, rho_ext_res*1e5, 'r--', 'LineWidth', 2.5, 'DisplayName', 'External Insertion', 'Marker', 'none'); hold on;
plot(t_out, rho_net_res*1e5, 'k-', 'LineWidth', 3, 'DisplayName', 'Net Reactivity');
yline(0, 'k:', 'LineWidth', 1.5, 'Alpha', 0.5, 'DisplayName', 'Zero Reactivity');
title('Reactivity Profile', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$\rho$ (pcm)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 11, 'Box', 'on', 'Color', 'white');
grid on; ax1.GridAlpha = 0.3; ax1.MinorGridAlpha = 0.15;
set(ax1, 'FontSize', 11);

% Plot 2: Neutron Flux
ax2 = subplot(4, 1, 2);
plot(t_out, phi_res, 'b-', 'LineWidth', 2.5);
title('Neutron Flux (Fission Chain Reaction)', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$\phi(t)$ / $\phi_0$', 'Interpreter', 'latex', 'FontSize', 12);
grid on; ax2.GridAlpha = 0.3; ax2.MinorGridAlpha = 0.15;
set(ax2, 'FontSize', 11);

% Plot 3: Thermal Power
ax3 = subplot(4, 1, 3);
plot(t_out, P_total_res / P0, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Total Thermal Power'); hold on;
plot(t_out, P_decay_total / P0, 'm--', 'LineWidth', 2.5, 'DisplayName', 'Decay Heat Contribution');
yline(1, 'k:', 'LineWidth', 1.5, 'Alpha', 0.5, 'DisplayName', 'Nominal Power');
title('Reactor Thermal Power', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$P(t)$ / $P_0$', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 11, 'Box', 'on', 'Color', 'white');
grid on; ax3.GridAlpha = 0.3; ax3.MinorGridAlpha = 0.15;
set(ax3, 'FontSize', 11);

% Plot 4: Core Temperatures
ax4 = subplot(4, 1, 4);
plot(t_out, T_F_res, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Fuel ($T_F$)'); hold on;
plot(t_out, T_C_res, 'c-', 'LineWidth', 2.5, 'DisplayName', 'Coolant ($T_C$)');
title('Core Temperatures', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Temperature ($^{\circ}$C)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 11, 'Box', 'on', 'Color', 'white');
grid on; ax4.GridAlpha = 0.3; ax4.MinorGridAlpha = 0.15;
set(ax4, 'FontSize', 11);


% Time-Dependent Reactivity
function rho = get_rho_ext(t)
    if t < 20
        % Ramp insertion: Pulling rods at 5 pcm per second
        rho = 5e-5 * t; 
    else
        % SCRAM: Rods drop at t=20s, inserting -5000 pcm
        rho = -5000e-5; 
    end
end

% ODE FUNCTION DEFINITION
function dYdt = reactor_dynamics(t, Y, beta, beta_i, lambda_i, td, ...
                                 P0, T_in, M_F_Cp_F, M_C_Cp_C, UA, m_dot_Cp_C, ...
                                 alpha_F, alpha_C, T_F0, T_C0, ...
                                 gamma_j, lambda_hj, prompt_fraction)
    % Unpack
    phi = Y(1);  C = Y(2:7);  T_F = Y(8);  T_C = Y(9);  P_d = Y(10:12);
    
    % 1. Reactivity
    rho_ext = get_rho_ext(t);
    rho_net = rho_ext + alpha_F * (T_F - T_F0) + alpha_C * (T_C - T_C0);
    
    % 2. Neutronics (6-group)
    dphi_dt = ((rho_net - beta) / td) * phi + sum(lambda_i' .* C);
    dC_dt   = (beta_i' / td) * phi - lambda_i' .* C;
    
    % 3. Decay Heat (3-group)
    dP_d_dt = lambda_hj' .* (gamma_j' * P0 * phi) - lambda_hj' .* P_d;
    
    % 4. Thermal-Hydraulics
    P_total = P0 * phi * prompt_fraction + sum(P_d);
    f_c = 0.025; % Fraction of heat deposited directly in coolant
    dT_F_dt = ((1 - f_c) * P_total - UA * (T_F - T_C)) / M_F_Cp_F;
    dT_C_dt = (f_c * P_total + UA * (T_F - T_C) - 2 * m_dot_Cp_C * (T_C - T_in)) / M_C_Cp_C;
    
    % Pack derivatives
    dYdt = [dphi_dt; dC_dt; dT_F_dt; dT_C_dt; dP_d_dt];
end
