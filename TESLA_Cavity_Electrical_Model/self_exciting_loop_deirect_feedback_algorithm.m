% MATLAB Script for Analyzing Self-Exciting Loop (SEL) Control (Cleaner Figure 3 with Energy Display)
% Based on "Complex Envelope Control..." document, Chapter 4.1

clearvars;
close all;
clc;

% --- Parameters ---
f_half_bandwidth = 217; % Cavity half-bandwidth (Hz)
omega_half_bandwidth = 2 * pi * f_half_bandwidth; % (rad/s)

Vc_target = 25e6; % Target V_cav amplitude (V).
tf_specific = 500e-6; % Specific filling time (s)

omega0_rho_val = 4.247e12; % Term from cavity state eq. (Ohm/s or V/(A*s))
R_L = omega0_rho_val / (2 * omega_half_bandwidth); % Loaded shunt impedance (Ohms)

% Time vector
t = linspace(0, tf_specific, 1000);
small_epsilon = 1e-12; % For numerical stability near t=0

disp('--- Analyzing Self-Exciting Loop (SEL) Control ---');
disp(['Target Vc = ' num2str(Vc_target/1e6) ' MV, Filling time tf = ' num2str(tf_specific*1e6) ' us.']);
disp('This script visualizes concepts related to SEL as described in the document.');
disp(' ');

% --- 1. Optimal Energy Efficiency vs. Filling Time (Target for SEL) ---
tf_range = linspace(100e-6, 2000e-6, 500);
eta_fo = 1 - exp(-2 * omega_half_bandwidth * tf_range); % Eq. 4.1-6

figure('Name', 'Optimal Energy Efficiency (Target for SEL)');
plot(tf_range * 1e6, eta_fo, 'LineWidth', 1.5);
title({'Optimal Energy Efficiency $\eta_{fo}$ vs. Filling Time $t_f$';'(This is the target efficiency SEL aims for)'},'Interpreter','latex');
xlabel('Filling Time $t_f$ ($\mu$s)','Interpreter','latex');
ylabel('Energy Efficiency $\eta_{fo}$','Interpreter','latex');
grid on;
ylim([0 1]);
text(tf_range(100)*1e6 + 50, eta_fo(100)+0.05, sprintf('f_{1/2} = %d Hz', f_half_bandwidth),'Interpreter','latex', 'FontSize',9);
disp('PLOT 1: Shows the optimal energy efficiency SEL aims to achieve as a function of filling time.');

% --- 2. SEL Gain Factor G_SEL(t) & Resonance Tracking Explanation ---
G_SEL_t = zeros(size(t));
for k = 1:length(t)
    if t(k) < small_epsilon
        G_SEL_t(k) = NaN;
    else
        G_SEL_t(k) = (2 * omega_half_bandwidth) / (1 - exp(-2 * omega_half_bandwidth * t(k)));
    end
end

figure('Name', 'SEL Gain Factor & Implication for Resonance Tracking');
plot(t * 1e6, G_SEL_t, 'LineWidth', 1.5);
title({'SEL Gain Factor $G_{SEL}(t)$ vs. Time'; '$u(t) = G_{SEL}(t)v(t)$'},'Interpreter','latex');
xlabel('Time ($\mu$s)','Interpreter','latex');
ylabel('$G_{SEL}(t)$ (1/s)','Interpreter','latex');
grid on;
first_valid_gain_idx = find(~isnan(G_SEL_t), 1, 'first');
current_ylim_gain = [0, 1e5]; % Default y-limit
if ~isempty(first_valid_gain_idx)
    val_at_first_valid = G_SEL_t(first_valid_gain_idx);
    if ~isinf(val_at_first_valid) && val_at_first_valid > 0
      current_ylim_gain = [0, min(val_at_first_valid * 1.5, 5e5)];
    end
end
ylim(current_ylim_gain);
disp('PLOT 2: Shows the SEL gain G_SEL(t).');
disp('  - Quote Insight: "The optimal control realized as a self-exciting loop... forces cavity itself to track its resonance frequency."');
disp('  - DEMONSTRATION: G_SEL(t) is real and positive for t > 0. From u(t) = G_SEL(t)v(t), complex u(t) and v(t) are co-phasal.');
disp('    If v(t) phase is correct for resonance, u(t) from SEL inherently has this resonant phase, thus "tracking".');
disp(' ');


% --- 3. Cavity Variables under SEL (Optimal Path) & Klystron Power Constraint ---
v_opt_t_abs = Vc_target * (sinh(omega_half_bandwidth * t) / sinh(omega_half_bandwidth * tf_specific));
v_opt_t_abs(1) = 0;

u_SEL_t_abs = omega_half_bandwidth * Vc_target * (exp(omega_half_bandwidth * t) / sinh(omega_half_bandwidth * tf_specific));
u_SEL_t_abs(1) = (omega_half_bandwidth * Vc_target) / sinh(omega_half_bandwidth * tf_specific);

i_g_t_SEL = u_SEL_t_abs / omega0_rho_val;
P_f_t_SEL = 0.5 * R_L * (i_g_t_SEL.^2);
P_r_t_SEL = (abs(v_opt_t_abs - R_L * i_g_t_SEL).^2) / (2 * R_L);

% Calculate Stored Energy and Efficiency for annotation
W_tf_SEL = (Vc_target^2) / (2 * omega0_rho_val); % Stored energy at t_f
expended_energy_SEL = trapz(t, P_f_t_SEL);
eta_calculated_SEL = W_tf_SEL / expended_energy_SEL;


% Plotting for specific case - FIGURE 3 (TRAJECTORY AND POWER VISUALIZATION)
fig_sel_details = figure('Name', 'SEL: Trajectories & Power (with Energy Display)');
fig_font_size = 8; % Reduced font size for this figure's elements

% Subplot 1: Cavity Voltage and Current
subplot(2,1,1);
yyaxis left;
plot(t * 1e6, v_opt_t_abs / 1e6, 'b-', 'LineWidth', 1);
hold on;
legend_entries_volt = {'$|v(t)|$'};
if Vc_target >= 20e6
    plot([t(1) t(end)]*1e6, [20 20], 'k--', 'LineWidth', 0.7);
    legend_entries_volt = {'$|v(t)|$', '20 MV ref.'};
end
hold off;
ylabel('$|v(t)|$ (MV)','Interpreter','latex', 'FontSize', fig_font_size);
ylim([0, Vc_target/1e6 * 1.05]);
set(gca, 'FontSize', fig_font_size);

yyaxis right;
plot(t * 1e6, i_g_t_SEL * 1e3, 'r-', 'LineWidth', 1);
ylabel('$i_g(t)$ (mA)','Interpreter','latex', 'FontSize', fig_font_size);
min_ig_plot = min(0, floor(min(i_g_t_SEL*1e3)*0.9));
ylim([min_ig_plot, ceil(max(i_g_t_SEL*1e3)*1.05)]);
legend_entries_volt{end+1} = '$i_g(t)$';
set(gca, 'FontSize', fig_font_size);


title(sprintf('SEL: $V_{cav}$ & $I_{gen}$ ($V_{target}=%.0f$MV)', Vc_target/1e6), 'Interpreter','latex', 'FontSize', fig_font_size+1);
xlabel('Time ($\mu$s)','Interpreter','latex', 'FontSize', fig_font_size);
lgd1 = legend(legend_entries_volt, 'Location', 'East', 'Interpreter','latex');
lgd1.FontSize = fig_font_size -1;
grid on;
ax_volt_curr = gca; ax_volt_curr.YAxis(1).Color = 'b'; ax_volt_curr.YAxis(2).Color = 'r';

% Subplot 2: Forward and Reflected Power
subplot(2,1,2);
plot(t * 1e6, P_f_t_SEL / 1e3, 'b-', 'LineWidth', 1);
hold on;
plot(t * 1e6, P_r_t_SEL / 1e3, 'r-', 'LineWidth', 1);
hold off;
title(sprintf('SEL: Forward & Reflected Power ($V_{target}=%.0f$MV)', Vc_target/1e6),'Interpreter','latex', 'FontSize', fig_font_size+1);
xlabel('Time ($\mu$s)','Interpreter','latex', 'FontSize', fig_font_size);
ylabel('Power (kW)','Interpreter','latex', 'FontSize', fig_font_size);
lgd2 = legend('$P_f(t)$', '$P_r(t)$', 'Location', 'West','Interpreter','latex');
lgd2.FontSize = fig_font_size - 1;
grid on;
ax_power = gca;
set(ax_power, 'FontSize', fig_font_size);

% Add Stored Energy and Efficiency to the power plot
text_y_pos1 = max(P_f_t_SEL/1e3) * 0.8;
text_y_pos2 = max(P_f_t_SEL/1e3) * 0.7;
if isempty(text_y_pos1) || isnan(text_y_pos1) || text_y_pos1 <=0, text_y_pos1 = 100; text_y_pos2 = 80; end % Fallback positions

text(tf_specific*1e6 * 0.05, text_y_pos1, sprintf('Stored Energy $W(t_f) = %.1f$ J', W_tf_SEL), ...
     'Interpreter','latex', 'FontSize', fig_font_size-1, 'BackgroundColor', [0.95 0.95 0.95]);
text(tf_specific*1e6 * 0.05, text_y_pos2, sprintf('Energy Efficiency $\\eta_f = %.3f$', eta_calculated_SEL), ...
     'Interpreter','latex', 'FontSize', fig_font_size-1, 'BackgroundColor', [0.95 0.95 0.95]);


sgtitle('SEL Control: Trajectories & Power (with Energy Display)', 'Interpreter', 'latex', 'FontSize', fig_font_size+2);

disp(' ');
disp('PLOT 3 (Cavity/Current & Power Visualization - with Energy Display):');
disp('  - Font sizes reduced for clarity.');
disp('  - Stored energy at t_f and overall filling energy efficiency are now displayed on the power subplot.');
disp('  - Top subplot: Cavity voltage |v(t)|. If Vc_target >= 20MV, a 20MV reference line is shown.');
disp('  - Bottom subplot: Forward Power Pf(t). Observe the power rise at end of fill for');
disp('    high gradients (e.g., >=20MV), illustrating the klystron power constraint.');
disp(' ');
disp('SEL Analysis Script Finished.');