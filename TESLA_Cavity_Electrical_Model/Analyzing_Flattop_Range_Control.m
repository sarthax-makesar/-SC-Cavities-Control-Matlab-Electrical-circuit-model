% MATLAB Script for Analyzing Flattop Range Control

clearvars;
close all;
clc;

%  Parameters 
f_hb = 217; % Cavity half-bandwidth (Hz)
w_hb = 2 * pi * f_hb; % omega_1/2 (rad/s)

Vc_mag = 25e6; % Target cavity voltage magnitude |Vc| (V)
phi_b = 0; % Relative beam phase (radians, 0 for on-crest)

w0_rho = 4.247e12; % Cavity parameter omega0*rho (Ohm/s or V/(A*s))

disp(' Flattop Range Control Analysis ');
disp(['Target |Vc| = ' num2str(Vc_mag/1e6) ' MV, Beam phase phi_b = ' num2str(phi_b) ' rad.']);
disp(' ');

%  Helper function for Flattop Pf and eta_b 
% Calculation of Pf and eta_b for given |Vc|, |Ib|, r_detuning, phi_b_in
function [Pf_calc, eta_b_calc, Pb_calc, g_factor] = calculate_flattop_metrics(Vc_abs, Ib_abs, r_detuning, phi_b_in, w_1_2_in, w0_rho_in)
    if Vc_abs == 0, Pf_calc=inf; eta_b_calc=0; Pb_calc=0; g_factor=inf; return; end
    
    Ub_mag = w0_rho_in * Ib_abs; % Magnitude of unified beam loading
    
    % Beam power
    Pb_calc = Vc_abs * Ib_abs * cos(phi_b_in); % P_b = |Vc| |Ib| cos(phi_b)
    
    % Factor g (Eq. 4.1-20 / 2.4-4 context)
    g_factor = (Ub_mag * cos(phi_b_in)) / (w_1_2_in * Vc_abs);
    if g_factor == 0 && Pb_calc ~=0 
        Pf_calc = inf; % If g is zero (no beam coupling effectively contributing to g), Pf would be large or undefined by this formula for finite Pb
        eta_b_calc = 0;
        if Pb_calc == 0, Pf_calc = 0; end % No beam, no beam power, Pf for sustaining field only (not covered by this formula)
        return;
    elseif Pb_calc == 0 % No beam power
        Ac_val = -w_1_2_in + 1j * r_detuning * w_1_2_in;
        u_g_complex = Ub_mag * exp(1j*phi_b_in) - Ac_val * Vc_abs; % Assuming Vc_abs is on real axis
        R_L_eff = w0_rho_in / (2*w_1_2_in); % Effective R_L from w0_rho = 2*w_1/2*R_L
        Pf_calc = 0.5 * R_L_eff * (abs(u_g_complex)/w0_rho_in)^2;
        eta_b_calc = 0;
        return;
    end

    % Here r is tan(psi) from Eq. 2.4-4, so r = DeltaOmega/omega_1/2
    Pf_calc = ( (1+g_factor)^2 + (r_detuning - g_factor*tan(phi_b_in))^2 ) / (4*g_factor) * Pb_calc;
    
    if Pf_calc == 0 || isnan(Pf_calc) || isinf(Pf_calc)
        eta_b_calc = 0; % Avoid division by zero or NaN
    else
        eta_b_calc = Pb_calc / Pf_calc;
    end
    if eta_b_calc < 0, eta_b_calc = 0; end % Efficiency cannot be negative
    if eta_b_calc > 1, eta_b_calc = 1; end % Max efficiency is 1
end

%   Pf and eta_b vs. Relative Detuning 'r' 
r_vals_plot1 = linspace(-1.5, 1.5, 201); % Relative detuning r = DeltaOmega / omega_1/2
Ib_plot1 = 8e-3; % Fixed beam current (A) for this plot, e.g., 8mA
Pf_vs_r = zeros(size(r_vals_plot1));
eta_b_vs_r = zeros(size(r_vals_plot1));

for i = 1:length(r_vals_plot1)
    [Pf_vs_r(i), eta_b_vs_r(i), ~, ~] = calculate_flattop_metrics(Vc_mag, Ib_plot1, r_vals_plot1(i), phi_b, w_hb, w0_rho);
end

figure('Name', 'Flattop: Pf and Efficiency vs. Detuning r');
% Subplot A: Forward Power vs. r
subplot(2,1,1);
plot(r_vals_plot1, Pf_vs_r / 1e3, 'b-', 'LineWidth', 1.5);
title({sprintf('Forward Power $P_f$ vs. Relative Detuning $r=\\Delta\\omega/\\omega_{1/2}$'); ...
       sprintf('$|V_c|=%.0f$MV, $I_b=%.0fmA, \\phi_b=%.1f$rad', Vc_mag/1e6, Ib_plot1*1e3, phi_b)}, 'Interpreter','latex');
xlabel('Relative Detuning $r$','Interpreter','latex');
ylabel('$P_f$ (kW)','Interpreter','latex');
grid on;
[min_Pf_r, idx_min_Pf_r] = min(Pf_vs_r);
hold on; plot(r_vals_plot1(idx_min_Pf_r), min_Pf_r/1e3, 'ro', 'MarkerFaceColor','r');
text(r_vals_plot1(idx_min_Pf_r), min_Pf_r/1e3 + 0.05*min_Pf_r/1e3, sprintf('Min $P_f \\approx %.1f$kW @ $r \\approx %.2f$', min_Pf_r/1e3, r_vals_plot1(idx_min_Pf_r)),'Interpreter','latex','FontSize',8); hold off;

% Subplot B: Efficiency vs. r
subplot(2,1,2);
plot(r_vals_plot1, eta_b_vs_r, 'r-', 'LineWidth', 1.5);
title({'Flattop Efficiency $\eta_b$ vs. Relative Detuning $r$'},'Interpreter','latex');
xlabel('Relative Detuning $r$','Interpreter','latex');
ylabel('Efficiency $\eta_b$','Interpreter','latex');
grid on; ylim([0 1.05]);
[max_eta_r, idx_max_eta_r] = max(eta_b_vs_r);
hold on; plot(r_vals_plot1(idx_max_eta_r), max_eta_r, 'bo', 'MarkerFaceColor','b');
text(r_vals_plot1(idx_max_eta_r), max_eta_r - 0.05, sprintf('Max $\\eta_b \\approx %.3f$ @ $r \\approx %.2f$', max_eta_r, r_vals_plot1(idx_max_eta_r)),'Interpreter','latex','FontSize',8); hold off;
disp('PLOT 1: Forward Power and Flattop Efficiency vs. Relative Detuning r.');

%  Pf and eta_b vs. Beam Current Ib 
Ib_vals_plot2 = linspace(0.1e-3, 20e-3, 200); % Beam current from 0.1mA to 20mA
r_plot2 = 0; % Fixed relative detuning for this plot (e.g., on resonance r=0)
% Can add another r value for comparison if needed
% r_plot2_off_res = 0.5; 

Pf_vs_Ib = zeros(size(Ib_vals_plot2));
eta_b_vs_Ib = zeros(size(Ib_vals_plot2));
% Pf_vs_Ib_off_res = zeros(size(Ib_vals_plot2));
% eta_b_vs_Ib_off_res = zeros(size(Ib_vals_plot2));


for i = 1:length(Ib_vals_plot2)
    [Pf_vs_Ib(i), eta_b_vs_Ib(i), ~, ~] = calculate_flattop_metrics(Vc_mag, Ib_vals_plot2(i), r_plot2, phi_b, w_hb, w0_rho);
%   [Pf_vs_Ib_off_res(i), eta_b_vs_Ib_off_res(i), ~, ~] = calculate_flattop_metrics(Vc_mag, Ib_vals_plot2(i), r_plot2_off_res, phi_b, w_hb, w0_rho);
end

figure('Name', 'Flattop: Pf and Efficiency vs. Beam Current');
% Subplot A: Forward Power vs. Ib
subplot(2,1,1);
plot(Ib_vals_plot2 * 1e3, Pf_vs_Ib / 1e3, 'b-', 'LineWidth', 1.5); hold on;
% plot(Ib_vals_plot2 * 1e3, Pf_vs_Ib_off_res / 1e3, 'g--', 'LineWidth', 1.5);
title({sprintf('Forward Power $P_f$ vs. Beam Current $I_b$'); 
       sprintf('$|V_c|=%.0f$MV, $r=%.1f, \\phi_b=%.1f$rad', Vc_mag/1e6, r_plot2, phi_b)},'Interpreter','latex');
xlabel('Beam Current $I_b$ (mA)','Interpreter','latex');
ylabel('$P_f$ (kW)','Interpreter','latex');
grid on; % legend({sprintf('r = %.1f',r_plot2), sprintf('r = %.1f',r_plot2_off_res)},'Location','NorthWest','Interpreter','latex'); 
hold off;

% Subplot B: Efficiency vs. Ib
subplot(2,1,2);
plot(Ib_vals_plot2 * 1e3, eta_b_vs_Ib, 'r-', 'LineWidth', 1.5); hold on;
% plot(Ib_vals_plot2 * 1e3, eta_b_vs_Ib_off_res, 'm--', 'LineWidth', 1.5);
title({'Flattop Efficiency $\eta_b$ vs. Beam Current $I_b$'},'Interpreter','latex');
xlabel('Beam Current $I_b$ (mA)','Interpreter','latex');
ylabel('Efficiency $\eta_b$','Interpreter','latex');
grid on; ylim([0 1.05]);
% legend({sprintf('r = %.1f',r_plot2), sprintf('r = %.1f',r_plot2_off_res)},'Location','NorthEast','Interpreter','latex');
hold off;
disp('PLOT 2: Forward Power and Flattop Efficiency vs. Beam Current.');

%  Note on Time-Varying Detuning 
disp(' ');
disp('NOTE: This script assumes constant detuning during flattop for the parametric plots.');
disp('For time-varying detuning DeltaOmega(t) during flattop, one would calculate u_g(t) ');
disp('over time, then Pf(t), and integrate Pf(t) to find average power for efficiency.');
