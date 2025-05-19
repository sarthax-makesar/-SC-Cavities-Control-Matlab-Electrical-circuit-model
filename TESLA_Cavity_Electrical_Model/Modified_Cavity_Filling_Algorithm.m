% MATLAB Script for Modified Cavity Filling Algorithms 


clearvars;
close all;
clc;

% --- Parameters ---
f_hb = 217; % Cavity half-bandwidth (Hz)
w_hb = 2 * pi * f_hb; % (rad/s)

Vc = 25e6; % Target V_cav amplitude (V)
tf_main = 500e-6; % Main filling time for detailed plots (s)

w0_rho = 4.247e12; % Cavity parameter omega0*rho (Ohm/s or V/(A*s))
R_L_cav = w0_rho / (2 * w_hb); % Loaded shunt impedance (Ohms)

% --- Global setting for r_detuning for Plots 2 and 3 (for non-optimal modified trajectories) ---
r_for_plots_2_3 = 0; % Default to resonant case (r=0)

disp('--- Modified Cavity Filling Algorithms (parameter b, varying detuning r, optimal modified trajectory) ---');
disp(['Vc_target = ' num2str(Vc/1e6) ' MV, Default tf = ' num2str(tf_main*1e6) ' us.']);
disp(['Plots 2 & 3 for standard modified trajectories will use r = ' num2str(r_for_plots_2_3) '.']);
disp(' ');

% --- Helper function for metrics given 'b', detuning 'r', and trajectory type ---

function [v_abs, u_abs, ig, Pf, Pr, W_stored, eta_f] = calc_metrics_b_algo(b, r_detuning, trajectory_type, t_vec, Vc_in, tf_in, w_1_2_in, w0_rho_in, R_L_in)
    v_abs = zeros(size(t_vec));
    u_complex = zeros(size(t_vec)); % u(t) can be complex

    % 1. Calculate |v(t)| based on trajectory_type
    if strcmp(trajectory_type, 'standard')
        if abs(b) < 1e-6 % b = 0 (Eq. 4.1-14)
            for i = 1:length(t_vec)
                if tf_in == 0, v_abs(i)=0; else, v_abs(i) = Vc_in * (t_vec(i) / tf_in); end
            end
        else % b ~= 0 (Eq. 4.1-13)
            den_v = 1 - exp(-b * w_1_2_in * tf_in);
            if abs(den_v) < 1e-9, v_abs(:)=0; else
                for i = 1:length(t_vec)
                    v_abs(i) = Vc_in * (1 - exp(-b * w_1_2_in * t_vec(i))) / den_v;
                end
            end
        end
    elseif strcmp(trajectory_type, 'optimal_modified_sinh') % Eq. 4.1-15
         if abs(b) < 1e-6 % b = 0 (same as t/tf)
            for i = 1:length(t_vec)
                if tf_in == 0, v_abs(i)=0; else, v_abs(i) = Vc_in * (t_vec(i) / tf_in); end
            end
        else % b ~= 0
            den_v_sinh = sinh(b * w_1_2_in * tf_in);
            if abs(den_v_sinh) < 1e-9, v_abs(:)=0; else
                for i = 1:length(t_vec)
                    v_abs(i) = Vc_in * sinh(b * w_1_2_in * t_vec(i)) / den_v_sinh;
                end
            end
        end
    else
        error('Unknown trajectory_type specified.');
    end
    if ~isempty(v_abs), v_abs(1) = 0; end

    % 2. Calculate u(t) (complex,  take abs value)
    if strcmp(trajectory_type, 'optimal_modified_sinh') % For this trajectory, u(t) is for r=0
        if abs(b) < 1e-6 % b = 0
            if tf_in == 0, u_complex(:)=0; else
                for i = 1:length(t_vec)
                    u_complex(i) = Vc_in * (1 + w_1_2_in * t_vec(i)) / tf_in; % Resonant u(t)
                end
            end
        else % b ~= 0, resonant condition: u = dv/dt + w_1/2*v
            den_v_sinh = sinh(b * w_1_2_in * tf_in);
            if abs(den_v_sinh) < 1e-9, u_complex(:)=0; else
                for i = 1:length(t_vec)
                    dvdt_term = Vc_in * b * w_1_2_in * cosh(b * w_1_2_in * t_vec(i)) / den_v_sinh;
                    w12v_term = w_1_2_in * v_abs(i); % v_abs is already calculated for this trajectory
                    u_complex(i) = dvdt_term + w12v_term;
                end
            end
        end
         if ~isempty(u_complex) % Set u(0) specifically for optimal_modified_sinh
            if abs(b) < 1e-6
                if tf_in~=0, u_complex(1) = Vc_in/tf_in; else u_complex(1)=0; end
            else
                den_v_sinh0 = sinh(b * w_1_2_in * tf_in);
                if abs(den_v_sinh0) > 1e-9
                    u_complex(1) = Vc_in * b * w_1_2_in * cosh(0) / den_v_sinh0 + w_1_2_in * v_abs(1);
                else
                    u_complex(1) = 0;
                end
            end
        end
    else % 'standard' trajectory, u(t) depends on r_detuning
        if abs(r_detuning) < 1e-6 % Resonant case u(t) (from Eq. 4.1-13, 4.1-14 magnitudes)
            if abs(b) < 1e-6 % b = 0
                if tf_in == 0, u_complex(:)=0; else
                    for i = 1:length(t_vec), u_complex(i) = Vc_in * (1 + w_1_2_in * t_vec(i)) / tf_in; end
                end
                if ~isempty(u_complex), if tf_in ~=0, u_complex(1) = Vc_in/tf_in; else u_complex(1)=0; end; end
            else % b ~= 0
                den_u = 1 - exp(-b * w_1_2_in * tf_in);
                if abs(den_u) < 1e-9, u_complex(:)=0; else
                    for i = 1:length(t_vec), u_complex(i) = w_1_2_in * Vc_in * (1 - (1 - b) * exp(-b * w_1_2_in * t_vec(i))) / den_u; end
                    if ~isempty(u_complex), if abs(den_u) > 1e-9, u_complex(1) = w_1_2_in * Vc_in * b / den_u; else u_complex(1)=0; end; end
                end
            end
        else % Non-resonant case u(t) (Eq. 4.1-17, 4.1-18)
            if abs(b) < 1e-6 % b = 0 (Eq. 4.1-18)
                if tf_in == 0, u_complex(:)=0; else
                    for i = 1:length(t_vec), u_complex(i) = Vc_in * (1 + w_1_2_in * (1 - 1j * r_detuning) * t_vec(i)) / tf_in; end
                end
                if ~isempty(u_complex), if tf_in ~=0, u_complex(1) = Vc_in/tf_in; else u_complex(1)=0; end; end
            else % b ~= 0 (Eq. 4.1-17)
                den_u = 1 - exp(-b * w_1_2_in * tf_in);
                if abs(den_u) < 1e-9, u_complex(:)=0; else
                    for i = 1:length(t_vec)
                        term_exp = exp(-b * w_1_2_in * t_vec(i));
                        num_u_cplx = (1 - 1j * r_detuning) - (1 - b - 1j * r_detuning) * term_exp;
                        u_complex(i) = w_1_2_in * Vc_in * num_u_cplx / den_u;
                    end
                    if ~isempty(u_complex)
                        if abs(den_u) > 1e-9
                             num_u0_cplx = (1 - 1j * r_detuning) - (1 - b - 1j * r_detuning);
                             u_complex(1) = w_1_2_in * Vc_in * num_u0_cplx / den_u;
                        else u_complex(1) = 0; end
                    end
                end
            end
        end
    end
    u_abs = abs(u_complex);

    % 3. Calculate derived metrics
    ig = u_abs / w0_rho_in;
    Pf = 0.5 * R_L_in * (ig.^2);
    Pr = (abs(v_abs - R_L_in * ig).^2) / (2 * R_L_in); % Simplified Pr
    W_stored = (Vc_in^2) / (2 * w0_rho_in);
    
    if isempty(t_vec) || length(t_vec) < 2, E_exp = inf; else E_exp = trapz(t_vec, Pf); end
    if E_exp == 0 || isinf(E_exp) || isnan(E_exp), eta_f = 0; else eta_f = W_stored / E_exp; end
end

% --- Plot: Efficiency vs. 'b' for different 'r' and optimal modified (Fig. 4.1-5 style) ---
b_vals_plot1 = linspace(-2, 2, 101);
t_for_plot1 = linspace(0, tf_main, 200); 
r_values_plot1 = [0, 0.5, 1]; % For 'standard' trajectory
r_labels_plot1 = {'$r=0$ (Std. Traj.)', '$r=0.5$ (Std. Traj.)', '$r=1$ (Std. Traj.)'};
colors_plot1 = lines(length(r_values_plot1) + 2); % +2 for optimal lines

figure('Name', 'Efficiency vs. b - Fig 4.1-5 Replication');
hold on;

% Theoretical Optimal Efficiency (Eq. 4.1-6)
eta_theory_opt_val = 1 - exp(-2*w_hb*tf_main);
plot(b_vals_plot1, ones(size(b_vals_plot1))*eta_theory_opt_val, 'k--', 'LineWidth', 1);
legend_entries_plot1 = {'Theor. Max $\eta_{fo}$ (Eq.4.1-6)'};
all_etas_plot1_curves = [];

% Modified Algorithm with v(t) from Eq. 4.1-13/14, for different 'r'
for j = 1:length(r_values_plot1)
    current_r = r_values_plot1(j);
    eta_f_b_plot1_std = zeros(size(b_vals_plot1));
    for i = 1:length(b_vals_plot1)
        [~,~,~,~,~,~,eta_f_b_plot1_std(i)] = calc_metrics_b_algo(b_vals_plot1(i), current_r, 'standard', t_for_plot1, Vc, tf_main, w_hb, w0_rho, R_L_cav);
    end
    plot(b_vals_plot1, eta_f_b_plot1_std, 'LineWidth', 1.5, 'Color', colors_plot1(j,:));
    legend_entries_plot1{end+1} = r_labels_plot1{j};
    all_etas_plot1_curves = [all_etas_plot1_curves; eta_f_b_plot1_std]; %#ok<AGROW>
    
    if abs(current_r) < 1e-6 % Mark max for resonant standard trajectory (r=0)
        [max_eta_b_r0_std, idx_max_eta_b_r0_std] = max(eta_f_b_plot1_std);
        plot(b_vals_plot1(idx_max_eta_b_r0_std), max_eta_b_r0_std, 'o', 'MarkerSize', 5, 'MarkerFaceColor',colors_plot1(j,:), 'MarkerEdgeColor',colors_plot1(j,:));
    end
end

% "Optimal Modified Solution" (Eq. 4.1-15 for v(t), resonant, red line in Fig 4.1-5)
eta_f_b_opt_mod = zeros(size(b_vals_plot1));
for i = 1:length(b_vals_plot1)
    % This trajectory is always considered under resonant conditions (r=0 for its u(t) calc)
    [~,~,~,~,~,~,eta_f_b_opt_mod(i)] = calc_metrics_b_algo(b_vals_plot1(i), 0, 'optimal_modified_sinh', t_for_plot1, Vc, tf_main, w_hb, w0_rho, R_L_cav);
end
plot(b_vals_plot1, eta_f_b_opt_mod, 'LineWidth', 1.5, 'Color', colors_plot1(length(r_values_plot1)+1,:));
legend_entries_plot1{end+1} = 'Opt. Mod. Traj. (Eq.4.1-15, r=0)';
all_etas_plot1_curves = [all_etas_plot1_curves; eta_f_b_opt_mod];

hold off;
title({'$\eta_f$ vs. $b$ for different trajectories & detuning $r$';sprintf('($t_f = %d \\mu s, V_c = %.0f$MV, Fig. 4.1-5 style)', round(tf_main*1e6), Vc/1e6)}, 'Interpreter','latex');
xlabel('Control Parameter $b$','Interpreter','latex'); ylabel('Energy Efficiency $\eta_f$','Interpreter','latex');
legend(legend_entries_plot1, 'Location','SouthOutside','Orientation','horizontal','Interpreter','latex', 'FontSize', 7); % Smaller font for legend
grid on; 
ylim_min_plot1 = 0.55;
current_max_plot1_data = 0.75;
if ~isempty(all_etas_plot1_curves), current_max_plot1_data = max(all_etas_plot1_curves(:)); end
ylim([ylim_min_plot1, max(0.76, current_max_plot1_data * 1.02)]);
disp('PLOT 1: Efficiency vs. parameter b, replicating Fig 4.1-5 lines.');


% --- Determine optimal 'b' for r=0 STANDARD trajectory for use in subsequent plots ---
temp_eta_for_opt_b_std_r0 = zeros(size(b_vals_plot1));
for k_opt_b = 1:length(b_vals_plot1)
    [~,~,~,~,~,~,temp_eta_for_opt_b_std_r0(k_opt_b)] = calc_metrics_b_algo(b_vals_plot1(k_opt_b), 0, 'standard', t_for_plot1, Vc, tf_main, w_hb, w0_rho, R_L_cav);
end
[~, idx_opt_b_r0_for_plots] = max(temp_eta_for_opt_b_std_r0);
optimal_b_for_std_r0_val = b_vals_plot1(idx_opt_b_r0_for_plots(1));


% --- Plot: Efficiency vs. Filling Time for selected 'b' (using r_for_plots_2_3, 'standard' trajectory) ---
tf_vals_plot2 = linspace(200e-6, 1200e-6, 100);
b_sel_plot2 = {0, 0.5, 1, optimal_b_for_std_r0_val};
b_lbl_plot2 = {'$b=0$', '$b=0.5$', '$b=1$', sprintf('$b\\approx%.2f$ (Near Opt. $\\eta_f$ for Std. r=0)',optimal_b_for_std_r0_val)};
colors_plot2 = lines(length(b_sel_plot2)+1);

figure('Name', ['Efficiency vs. tf (Std. Traj., r=' num2str(r_for_plots_2_3) ')']);
hold on;
eta_opt_tf_plot2 = 1 - exp(-2*w_hb*tf_vals_plot2);
plot(tf_vals_plot2*1e6, eta_opt_tf_plot2, '--', 'Color', colors_plot2(1,:), 'LineWidth', 1.5);
leg_plot2 = {'Theor. Max $\eta_{fo}$ (r=0)'};

t_pts_calc_plot2 = linspace(0,1,200);
for i = 1:length(b_sel_plot2)
    eta_f_tf_plot2 = zeros(size(tf_vals_plot2));
    for j = 1:length(tf_vals_plot2)
        tf_curr = tf_vals_plot2(j); if tf_curr < 1e-7, eta_f_tf_plot2(j)=0; continue; end
        t_act = t_pts_calc_plot2 * tf_curr;
        [~,~,~,~,~,~,eta_f_tf_plot2(j)] = calc_metrics_b_algo(b_sel_plot2{i}, r_for_plots_2_3, 'standard', t_act, Vc, tf_curr, w_hb, w0_rho, R_L_cav);
    end
    plot(tf_vals_plot2*1e6, eta_f_tf_plot2, 'Color', colors_plot2(i+1,:), 'LineWidth', 1.5);
    leg_plot2{end+1} = b_lbl_plot2{i};
end
hold off;
title({sprintf('$\\eta_f$ vs. $t_f$ for Std. Traj. & sel. $b$ ($r=%.1f$)',r_for_plots_2_3);'(Fig 4.1-4 style)'},'Interpreter','latex');
xlabel('$t_f$ ($\mu$s)','Interpreter','latex'); ylabel('$\eta_f$','Interpreter','latex');
legend(leg_plot2, 'Location','SouthEast','Interpreter','latex', 'FontSize', 8);
grid on; ylim([0.3 1]);
disp(['PLOT 2: Efficiency vs. tf for selected b (Standard Traj., r = ' num2str(r_for_plots_2_3) ').']);

% --- Detailed plots for specific 'b' (Std. Traj., using r_for_plots_2_3) ---
b_sel_plot3 = {optimal_b_for_std_r0_val, 0, 0.5, 1};
b_lbl_plot3 = {sprintf('$b\\approx%.2f$',optimal_b_for_std_r0_val), '$b=0$', '$b=0.5$', '$b=1$'};
t_plot3 = linspace(0, tf_main, 500);
ls_plot3 = {'-', '--', ':', '-.'}; 
c_plot3 = lines(length(b_sel_plot3));
plot_font_size = 8;

figure('Name', ['Detailed Trajectories (Std. Traj., r=' num2str(r_for_plots_2_3) ')']);
fprintf('\n--- Stored Energy & Efficiency for Detailed ''b'' Cases (Std. Traj., tf=%.0fus, r=%.1f) ---\n', tf_main*1e6, r_for_plots_2_3);

% Voltage subplot
subplot(3,1,1); hold on;
for i=1:length(b_sel_plot3)
    [v,~,~,~,~,~,~] = calc_metrics_b_algo(b_sel_plot3{i}, r_for_plots_2_3, 'standard', t_plot3, Vc, tf_main, w_hb, w0_rho, R_L_cav);
    plot(t_plot3*1e6, v/1e6, 'LineWidth',1, 'LineStyle', ls_plot3{i}, 'Color', c_plot3(i,:));
end
hold off; grid on;
title(sprintf('$|v(t)|$ ($t_f=%d\\mu s, V_c=%.0f$MV, Std.Traj, $r=%.1f$)', round(tf_main*1e6), Vc/1e6, r_for_plots_2_3),'Interpreter','latex','FontSize',plot_font_size+1);
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',plot_font_size); ylabel('$|v(t)|$ (MV)','Interpreter','latex','FontSize',plot_font_size);
legend(b_lbl_plot3, 'Location','NorthWest','Interpreter','latex','FontSize',plot_font_size-1); set(gca,'FontSize',plot_font_size);

% Current subplot
subplot(3,1,2); hold on;
for i=1:length(b_sel_plot3)
    [~,~,ig,~,~,~,~] = calc_metrics_b_algo(b_sel_plot3{i}, r_for_plots_2_3, 'standard', t_plot3, Vc, tf_main, w_hb, w0_rho, R_L_cav);
    plot(t_plot3*1e6, ig*1e3, 'LineWidth',1, 'LineStyle', ls_plot3{i}, 'Color', c_plot3(i,:));
end
hold off; grid on;
title(sprintf('$i_g(t)$ (Std. Traj., $r=%.1f$)',r_for_plots_2_3),'Interpreter','latex','FontSize',plot_font_size+1);
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',plot_font_size); ylabel('$i_g(t)$ (mA)','Interpreter','latex','FontSize',plot_font_size);
legend(b_lbl_plot3, 'Location','NorthWest','Interpreter','latex','FontSize',plot_font_size-1); set(gca,'FontSize',plot_font_size);

% Power subplot
subplot(3,1,3); hold on;
leg_power_handles = [];
for i=1:length(b_sel_plot3)
    [~,~,~,Pf,Pr,Wst,eta] = calc_metrics_b_algo(b_sel_plot3{i}, r_for_plots_2_3, 'standard', t_plot3, Vc, tf_main, w_hb, w0_rho, R_L_cav);
    p_main = plot(t_plot3*1e6, Pf/1e3, 'LineWidth',1, 'LineStyle', ls_plot3{i}, 'Color', c_plot3(i,:));
    leg_power_handles(end+1) = p_main; %#ok<AGROW>
    plot(t_plot3*1e6, Pr/1e3, 'LineWidth',0.5, 'LineStyle', ls_plot3{i}, 'Color', [c_plot3(i,:)*0.7, 0.5]);
    fprintf('b=%s (Std.Traj, r=%.1f): Wst=%.1f J, eta_f=%.4f\n', strrep(b_lbl_plot3{i},'$',''), r_for_plots_2_3, Wst, eta);
end
hold off; grid on;
title(sprintf('$P_f(t)$ (solid) & $P_r(t)$ (dimmed) (Std. Traj, $r=%.1f$)',r_for_plots_2_3),'Interpreter','latex','FontSize',plot_font_size+1);
xlabel('Time ($\mu$s)','Interpreter','latex','FontSize',plot_font_size); ylabel('Power (kW)','Interpreter','latex','FontSize',plot_font_size);
legend(leg_power_handles, b_lbl_plot3, 'Location','NorthWest','Interpreter','latex','FontSize',plot_font_size-1); set(gca,'FontSize',plot_font_size);

sgtitle(sprintf('Modified Algorithms (Std. Traj. for Plots 2&3 use $r=%.1f$)',r_for_plots_2_3), 'Interpreter','latex','FontSize',plot_font_size+2);
disp(['PLOT 3: Detailed trajectories (Standard Traj., using r = ' num2str(r_for_plots_2_3) ').']);
disp('  - Stored energy & efficiency for these cases are printed above.');