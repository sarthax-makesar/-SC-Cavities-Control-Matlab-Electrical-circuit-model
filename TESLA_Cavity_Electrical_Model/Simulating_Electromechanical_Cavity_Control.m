% MATLAB Script for Electromechanical Cavity Control (Fig. 4.2-2 style - Reworked)
% Focus on one dominant inertial mode + non-inertial for mechanical model

clearvars; close all; clc;

disp('--- Electromechanical Model Control Simulation (Reworked) ---');

% --- Simulation Parameters ---
T_s = 1e-6;         % Sampling time (s)
t_fill = 500e-6;    % Filling duration (s)
t_flattop = 800e-6; % Flattop duration (s)
t_decay = 700e-6;   % Decay duration (s)
t_pulse_end = t_fill + t_flattop + t_decay;
time_vec_sim = (0:T_s:t_pulse_end-T_s)';
n_pts = length(time_vec_sim);

n_fill_end = round(t_fill / T_s);
n_flattop_end = n_fill_end + round(t_flattop / T_s);

% --- Electrical Parameters ---
f_hb = 217;           % Half-bandwidth (Hz)
w_hb = 2 * pi * f_hb;   % omega_1/2 (rad/s)
w0_rho = 4.247e12;    % omega0*rho (Ohm/s or V/(A*s))

% --- Simulation Cases (Vc, Predetuning, Beam Current) ---
Vc_lvls_MV = [15, 20, 25];       % MV
df0_lvls_Hz = [130, 230, 370];   % Approx. optimal predetuning (Hz)
Ib_lvls_mA = [5.0, 6.5, 8.0];   % Approx. optimal beam current (mA)
num_sim_cases = length(Vc_lvls_MV);
sim_results = cell(num_sim_cases, 1);

% --- Mechanical Model Constants ---
L_cav_sq = (1.037)^2; 
K_L_Hz_per_MVm2_all = [0.1, 0.1, 0.1, 0.5]; % K factor Hz/(MV/m)^2
km_Hz_per_MV2_all = K_L_Hz_per_MVm2_all / L_cav_sq; % k factor Hz/MV^2

% Using ONE dominant inertial mode + non-inertial
fm_mech_Hz_inertial = 235; % Hz (first mode from Table 3.1-1)
wm_mech_rs_inertial = 2 * pi * fm_mech_Hz_inertial;
Qm_mech_inertial = 100;
km_Hz_per_MV2_inertial = km_Hz_per_MV2_all(1);

tau_noni_s = 0.1e-3; 
k_noni_Hz_MV2 = km_Hz_per_MV2_all(4);

% --- Main Simulation Loop ---
for i_case = 1:num_sim_cases
    Vc_mag_target = Vc_lvls_MV(i_case) * 1e6; 
    df0_Hz = df0_lvls_Hz(i_case);             
    Ib_mag_flattop = Ib_lvls_mA(i_case) * 1e-3; 
    
    Vc_phase_target = 0; 
    Vc_cplx_target = Vc_mag_target * exp(1j * Vc_phase_target);
    dw0_rs = 2 * pi * df0_Hz; 
    phi_b_rel = 0; 
    Ub_phase = Vc_phase_target + phi_b_rel;
    Ub_cplx_flattop = w0_rho * Ib_mag_flattop * exp(1j * Ub_phase);

    fprintf('Simulating Case %d: Vc=%.0fMV, df0=%.0fHz, Ib=%.1fmA\n', i_case, Vc_lvls_MV(i_case), df0_Hz, Ib_lvls_mA(i_case));

    % Initialize states per case
    W_m_inert_single = [0, 0]; % [delta_w_m, dot_delta_w_m] for the single inertial mode
    dw_noni_rs = 0;
    v_cav_cplx = zeros(n_pts, 1);
    u_drive_Vps = zeros(n_pts, 1);
    i_total_abs_A = zeros(n_pts, 1);
    phi_cav_target_fill = Vc_phase_target; % Initial target phase for filling logic
                                       % This needs to be chosen to make phase at t_f = Vc_phase_target
                                       % For simplicity, start with Vc_phase_target and let it evolve via dphi/dt = DeltaOmega
    
    dw_LFD_rs = zeros(n_pts, 1);
    dw_total_rs = zeros(n_pts, 1);

    for n = 1:(n_pts-1)
        v_n_V_abs = abs(v_cav_cplx(n));
        v_n_sq_MV2 = (v_n_V_abs / 1e6)^2;

        % Simplified Mechanical Model Update (1 inertial + 1 non-inertial)
        % Inertial Mode
        km_Hz_V2_i = km_Hz_per_MV2_inertial / (1e6)^2;
        Am_m = [0, 1; -wm_mech_rs_inertial^2, -wm_mech_rs_inertial/Qm_mech_inertial];
        Bm_m_term = -2 * pi * wm_mech_rs_inertial^2 * km_Hz_V2_i;
        Bm_m = [0; Bm_m_term];
        W_dot = Am_m * W_m_inert_single' + Bm_m * v_n_V_abs^2;
        W_m_inert_single = W_m_inert_single + (W_dot * T_s)';
        dw_inert_LFD_rs = W_m_inert_single(1);
        
        % Non-Inertial Mode
        km_noni_Hz_V2_val = k_noni_Hz_MV2 / (1e6)^2;
        drv_term_noni = -2 * pi * km_noni_Hz_V2_val * v_n_V_abs^2;
        dot_dw_noni = (drv_term_noni - dw_noni_rs) / tau_noni_s;
        dw_noni_rs = dw_noni_rs + dot_dw_noni * T_s;
        
        dw_LFD_rs(n) = dw_inert_LFD_rs + dw_noni_rs;
        dw_total_rs(n) = dw0_rs + dw_LFD_rs(n);

        % Electrical Control Input
        if n <= n_fill_end % FILLING
            % Target phase for u_drive should be such that d(phi_cav)/dt = delta_w_total
            % Let the phase of u_drive be phi_u
            % v_cav_phase_current_step = angle(v_cav_cplx(n)); % This can be noisy if v_cav is small
            % phi_u = v_cav_phase_current_step + delta_w_total_rs(n) * T_s; % Predict next phase for resonance
            % Alternative: Maintain a target phase accumulator for u_drive
            phi_cav_target_fill = phi_cav_target_fill + dw_total_rs(n) * T_s;

            t_fill_curr = (n-1)*T_s;
            sinh_tf = sinh(w_hb * t_fill);
            u_mag = 0;
            if abs(sinh_tf) > 1e-9
                u_mag = w_hb * Vc_mag_target * exp(w_hb * t_fill_curr) / sinh_tf;
            end
            u_drive_Vps(n) = u_mag * exp(1j * phi_cav_target_fill);
            i_total_abs_A(n) = abs(u_drive_Vps(n)) / w0_rho;
        elseif n > n_fill_end && n <= n_flattop_end % FLATTOP
            Ac_n = -w_hb + 1j * dw_total_rs(n);
            u_drive_Vps(n) = -Ac_n * Vc_cplx_target; % This is omega0*rho*(ig-ib)
            i_total_abs_A(n) = abs(u_drive_Vps(n)) / w0_rho; % This is |ig-ib|
        else % DECAY
            u_drive_Vps(n) = 0;
            i_total_abs_A(n) = 0;
        end
        
        % Electrical Model Update
        An_elec = (1 - w_hb * T_s) + 1j * dw_total_rs(n) * T_s;
        un_elec_drive = T_s * u_drive_Vps(n);
        if n < n_pts
            v_cav_cplx(n+1) = An_elec * v_cav_cplx(n) + un_elec_drive;
        end
    end
    dw_LFD_rs(n_pts) = dw_LFD_rs(n_pts-1); 
    dw_total_rs(n_pts) = dw_total_rs(n_pts-1);
    i_total_abs_A(n_pts) = i_total_abs_A(n_pts-1);
    u_drive_Vps(n_pts) = u_drive_Vps(n_pts-1);

    sim_results{i_case}.v_cav_abs = abs(v_cav_cplx);
    sim_results{i_case}.v_cav_phase = unwrap(angle(v_cav_cplx));
    sim_results{i_case}.i_total_abs_A = i_total_abs_A;
    u_drive_phase_temp = unwrap(angle(u_drive_Vps));
    u_drive_phase_temp(abs(u_drive_Vps) < 1e-9) = NaN;
    sim_results{i_case}.u_drive_phase = u_drive_phase_temp;
    sim_results{i_case}.dw_total_rs = dw_total_rs;
    sim_results{i_case}.Vc_level_MV = Vc_lvls_MV(i_case);
end
disp('All simulations finished.');

% --- Plotting (Overlaying results) ---
plot_fs = 8; 
line_styles = {'-', '--', ':'};
plot_colors = lines(num_sim_cases);
Vc_legend = arrayfun(@(x) sprintf('%dMV',x), Vc_lvls_MV, 'UniformOutput', false);

figure('Name', 'Electromechanical Sim (Multi Vc - Reworked)');

subplot(2,2,1); hold on;
for i=1:num_sim_cases, plot(time_vec_sim*1e6, sim_results{i}.v_cav_abs/1e6, 'LineWidth',1, 'Color',plot_colors(i,:), 'LineStyle',line_styles{i}); end
hold off; title('$|V_{cav}(t)|$','Int','latex','FS',plot_fs+1); xlabel('Time ($\mu$s)','Int','latex','FS',plot_fs); ylabel('$|V_{cav}|$ (MV)','Int','latex','FS',plot_fs);
grid on; set(gca,'FS',plot_fs); xlim([0 t_pulse_end*1e6]); legend(Vc_legend,'Loc','SE','FS',plot_fs-1);

subplot(2,2,2); hold on;
for i=1:num_sim_cases, plot(time_vec_sim*1e6, sim_results{i}.i_total_abs_A*1e3, 'LineWidth',1, 'Color',plot_colors(i,:), 'LineStyle',line_styles{i}); end
hold off; title('$|i_{total}(t)|$','Int','latex','FS',plot_fs+1); xlabel('Time ($\mu$s)','Int','latex','FS',plot_fs); ylabel('Current (mA)','Int','latex','FS',plot_fs);
grid on; set(gca,'FS',plot_fs); xlim([0 t_pulse_end*1e6]); legend(Vc_legend,'Loc','NE','FS',plot_fs-1);

subplot(2,2,3); hold on;
for i=1:num_sim_cases
    plot(time_vec_sim*1e6, sim_results{i}.v_cav_phase, 'LineWidth',1, 'Color',plot_colors(i,:), 'LineStyle',line_styles{i});
end
hold off; title('Phase of $V_{cav}$','Int','latex','FS',plot_fs+1); xlabel('Time ($\mu$s)','Int','latex','FS',plot_fs); ylabel('Phase (rad)','Int','latex','FS',plot_fs);
grid on; set(gca,'FS',plot_fs); xlim([0 t_pulse_end*1e6]); legend(Vc_legend,'Loc','SE','FS',plot_fs-1);

subplot(2,2,4); hold on;
for i=1:num_sim_cases, plot(time_vec_sim*1e6, sim_results{i}.dw_total_rs/(2*pi), 'LineWidth',1, 'Color',plot_colors(i,:), 'LineStyle',line_styles{i}); end
hold off; title('Total Cavity Detuning $\Delta f(t)$','Int','latex','FS',plot_fs+1); xlabel('Time ($\mu$s)','Int','latex','FS',plot_fs); ylabel('$\Delta f$ (Hz)','Int','latex','FS',plot_fs);
grid on; set(gca,'FS',plot_fs); xlim([0 t_pulse_end*1e6]); legend(Vc_legend,'Loc','NE','FS',plot_fs-1);

sgtitle('Electromechanical Model ($V_{flattop}$ varied, Simplified Mech. Model)', 'Int','latex', 'FS', plot_fs+2);
disp('Plots generated using a simplified mechanical model (1 inertial + 1 non-inertial mode).');