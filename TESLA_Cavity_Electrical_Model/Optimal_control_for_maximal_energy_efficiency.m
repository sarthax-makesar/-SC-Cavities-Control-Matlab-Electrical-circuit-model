% MATLAB Script for Optimal Cavity Filling Characteristics 

clearvars;
close all;
clc;

%  Parameters 
f_half_bandwidth = 217; % Cavity half-bandwidth (Hz)
omega_half_bandwidth = 2 * pi * f_half_bandwidth; % (rad/s)

Vc_target = 25e6; % Target V_cav amplitude (V), e.g., 25 MV
tf_specific = 500e-6; % Specific filling time (s), e.g., 500 us

% omega0_rho_val: Term from cavity state eq. (Ohm/s or V/(A*s))
omega0_rho_val = 4.247e12;
% R_L: Loaded shunt impedance (Ohms)
R_L = omega0_rho_val / (2 * omega_half_bandwidth);

%  1. Optimal Energy Efficiency vs. Filling Time 
tf_range = linspace(100e-6, 2000e-6, 500);
eta_fo = 1 - exp(-2 * omega_half_bandwidth * tf_range); 

figure('Name', 'Optimal Energy Efficiency');
plot(tf_range * 1e6, eta_fo, 'LineWidth', 1.5);
title('Optimal Energy Efficiency $\eta_{fo}$ vs. Filling Time $t_f$','Interpreter','latex');
xlabel('Filling Time $t_f$ ($\mu$s)','Interpreter','latex');
ylabel('Energy Efficiency $\eta_{fo}$','Interpreter','latex');
grid on;
ylim([0 1]);
text(tf_range(100)*1e6, eta_fo(100)+0.05, sprintf('f_{1/2} = %d Hz', f_half_bandwidth), 'FontSize', 9);

%  2. Specific Optimal Filling Case 
t = linspace(0, tf_specific, 1000); % Time vector

% Optimal cavity voltage magnitude |v(t)| 
v_t_abs = Vc_target * (sinh(omega_half_bandwidth * t) / sinh(omega_half_bandwidth * tf_specific));
v_t_abs(1) = 0;

% Optimal unified input magnitude |u(t)| [V/s]
u_t_abs = omega_half_bandwidth * Vc_target * (exp(omega_half_bandwidth * t) / sinh(omega_half_bandwidth * tf_specific));

% Generator current i_g(t) [A]
i_g_t = u_t_abs / omega0_rho_val;

% Forward Power P_f(t) [W]
P_f_t = 0.5 * R_L * (i_g_t.^2);

% Reflected Power P_r(t) [W]
P_r_t = (abs(v_t_abs - R_L * i_g_t).^2) / (2 * R_L);

% Stored Energy W_tf at end of fill [J]
W_tf = (Vc_target^2) / (2 * omega0_rho_val);
eta_calculated = W_tf / trapz(t, P_f_t); %  efficiency Calculation

% Plotting for specific case
figure('Name', 'Optimal Cavity Filling Details');

% Subplot 1: Cavity Voltage and Current
subplot(2,1,1);
yyaxis left;
plot(t * 1e6, v_t_abs / 1e6, 'b-', 'LineWidth', 1.5);
ylabel('$|v(t)|$ (MV)','Interpreter','latex');
ylim([0, Vc_target/1e6 * 1.05]);

yyaxis right;
plot(t * 1e6, i_g_t * 1e3, 'r-', 'LineWidth', 1.5);
ylabel('$i_g(t)$ (mA)','Interpreter','latex');
min_ig_plot = min(0, floor(min(i_g_t*1e3)*0.9));
ylim([min_ig_plot, ceil(max(i_g_t*1e3)*1.05)]);

title(sprintf('Optimal Filling: $V_{cav}$ & $I_{gen}$ ($t_f = %d \\mu s$)', round(tf_specific*1e6)), 'Interpreter','latex');
xlabel('Time ($\mu$s)','Interpreter','latex');
legend('$|v(t)|$', '$i_g(t)$', 'Location', 'East', 'Interpreter','latex');
grid on;
ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'r';

% Subplot 2: Forward and Reflected Power
subplot(2,1,2);
plot(t * 1e6, P_f_t / 1e3, 'b-', 'LineWidth', 1.5);
hold on;
plot(t * 1e6, P_r_t / 1e3, 'r-', 'LineWidth', 1.5);
hold off;
title(sprintf('Optimal Filling: Power ($t_f = %d \\mu s$)', round(tf_specific*1e6)),'Interpreter','latex');
xlabel('Time ($\mu$s)','Interpreter','latex');
ylabel('Power (kW)','Interpreter','latex');
legend('$P_f(t)$', '$P_r(t)$', 'Location', 'West','Interpreter','latex');
grid on;
max_Pf_plot = max(P_f_t/1e3); if isempty(max_Pf_plot) || max_Pf_plot==0, max_Pf_plot=1; end
text(max(t)*1e6*0.05, max_Pf_plot*0.8, sprintf('$W(t_f) = %.1f$ J', W_tf),'Interpreter','latex', 'FontSize',9);
text(max(t)*1e6*0.05, max_Pf_plot*0.65, sprintf('$\\eta_{fo} = %.3f$', eta_calculated),'Interpreter','latex', 'FontSize',9);

sgtitle('Optimal Cavity Filling Characteristics', 'Interpreter', 'latex', 'FontSize', 14);

disp('Script finished. Plots generated.');
disp(['Parameters used: Vc_target = ' num2str(Vc_target/1e6) ' MV, tf_specific = ' num2str(tf_specific*1e6) ' us.']);
disp(['Calculated efficiency for specific case: ' num2str(eta_calculated, '%.4f')]);
