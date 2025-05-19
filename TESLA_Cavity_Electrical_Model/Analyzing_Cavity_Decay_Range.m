% MATLAB Script for Analyzing Cavity Decay Range

clearvars;
close all;
clc;

% --- Parameters ---
f_hb_default = 217; % Default Cavity half-bandwidth (Hz)
w_hb_default = 2 * pi * f_hb_default; % Default omega_1/2 (rad/s)

V_start_mag = 25e6; % Initial cavity voltage magnitude at start of decay (V)
phi_start_rad = 0;  % Initial cavity voltage phase at start of decay (rad)

decay_duration = 2e-3; % Duration of decay to observe (s), e.g., 2ms
num_points = 1000;    % Number of points for time vector
t_decay = linspace(0, decay_duration, num_points);

disp(' Cavity Decay Range Analysis');
disp(['Initial |V_start| = ' num2str(V_start_mag/1e6) ' MV, Initial Phase = ' num2str(phi_start_rad) ' rad.']);
disp(['Observing decay for ' num2str(decay_duration*1e3) ' ms.']);
disp(' ');

%  Case 1: Effect of omega_1/2 on Amplitude Decay (DeltaOmega = 0)
delta_w_case1_rad_s = 0; % Zero detuning for this case
f_hb_values_Hz = [100, 217, 400]; % Different half-bandwidths to test
colors_case1 = lines(length(f_hb_values_Hz));

figure('Name', 'Decay: Amplitude vs. Time for different omega_1/2');
subplot(2,1,1); hold on;
title_text_case1_mag = {'$|v(t)|$ during Decay for different $\omega_{1/2}$';'($\Delta\omega = 0$)'};
title(title_text_case1_mag, 'Interpreter','latex');
xlabel('Time (ms)','Interpreter','latex'); 
ylabel('$|v(t)|$ (MV)','Interpreter','latex');
grid on;
legend_entries_case1 = {};

subplot(2,1,2); hold on;
title_text_case1_log = {'$\ln|v(t)|$ during Decay (Slope $\approx -\omega_{1/2}$)';'($\Delta\omega = 0$)'};
title(title_text_case1_log, 'Interpreter','latex');
xlabel('Time (ms)','Interpreter','latex'); 
ylabel('$\ln|v(t)|$','Interpreter','latex');
grid on;

for i = 1:length(f_hb_values_Hz)
    current_f_hb = f_hb_values_Hz(i);
    current_w_hb = 2 * pi * current_f_hb;
    
    v_mag_t = V_start_mag * exp(-current_w_hb * t_decay);
    ln_v_mag_t = log(v_mag_t); % log is natural log in MATLAB
    
    subplot(2,1,1);
    plot(t_decay * 1e3, v_mag_t / 1e6, 'LineWidth', 1.5, 'Color', colors_case1(i,:));
    
    subplot(2,1,2);
    plot(t_decay * 1e3, ln_v_mag_t, 'LineWidth', 1.5, 'Color', colors_case1(i,:));
    
    legend_entries_case1{end+1} = sprintf('$f_{1/2} = %d$ Hz', current_f_hb);
end
subplot(2,1,1); legend(legend_entries_case1,'Location','NorthEast','Interpreter','latex', 'FontSize', 8);
subplot(2,1,2); legend(legend_entries_case1,'Location','NorthEast','Interpreter','latex', 'FontSize', 8);
disp('PLOT 1: Shows amplitude decay for different half-bandwidths (f_1/2).');
disp('  - Top: |v(t)|. Bottom: ln|v(t)|, where slope illustrates -omega_1/2.');

%  Case 2: Effect of constant DeltaOmega on Phase & I/Q Plot 
delta_w_values_Hz_case2 = [-200, 0, 200]; % Hz; (0 for reference, +/- for detuning)
colors_case2 = lines(length(delta_w_values_Hz_case2));
w_hb_case2 = w_hb_default; % Use default omega_1/2

figure('Name', 'Decay: Phase Evolution & I/Q Plot for constant DeltaOmega');
sgtitle(sprintf('Decay Characteristics for Constant $\\Delta f$ ($f_{1/2}=%.0f$Hz)', w_hb_case2/(2*pi)), 'Interpreter','latex');

subplot(1,2,1); hold on; % Phase vs Time
title('Phase Evolution $\angle v(t)$','Interpreter','latex');
xlabel('Time (ms)','Interpreter','latex'); 
ylabel('Phase (rad)','Interpreter','latex');
grid on;
legend_entries_case2_phase = {};

subplot(1,2,2); hold on; % I/Q plot
title('I/Q Trajectory of $v(t)$','Interpreter','latex');
xlabel('In-Phase (I) [MV]','Interpreter','latex'); 
ylabel('Quadrature (Q) [MV]','Interpreter','latex');
grid on; axis equal;
legend_entries_case2_iq = {};

v_mag_t_case2 = V_start_mag * exp(-w_hb_case2 * t_decay); % Magnitude decays same way

for i = 1:length(delta_w_values_Hz_case2)
    current_delta_f = delta_w_values_Hz_case2(i);
    current_delta_w = 2 * pi * current_delta_f;
    
    % Phase evolution: phi(t) = phi_start + delta_w * t (since delta_w is constant)
    v_phase_t = phi_start_rad + current_delta_w * t_decay;
    
    v_complex_t = v_mag_t_case2 .* exp(1j * v_phase_t);
    
    subplot(1,2,1); % Phase vs Time
    plot(t_decay * 1e3, unwrap(v_phase_t), 'LineWidth', 1.5, 'Color', colors_case2(i,:));
    legend_entries_case2_phase{end+1} = sprintf('$\\Delta f = %d$ Hz', current_delta_f);
    
    subplot(1,2,2); % I/Q plot
    plot(real(v_complex_t) / 1e6, imag(v_complex_t) / 1e6, 'LineWidth', 1.5, 'Color', colors_case2(i,:));
    legend_entries_case2_iq{end+1} = sprintf('$\\Delta f = %d$ Hz', current_delta_f);
end
subplot(1,2,1); legend(legend_entries_case2_phase,'Location','NorthWest','Interpreter','latex', 'FontSize', 8);
subplot(1,2,2); legend(legend_entries_case2_iq,'Location','NorthEast','Interpreter','latex', 'FontSize', 8);
disp('PLOT 2: Shows phase evolution and I/Q trajectory for different constant detunings (Delta_f).');

%  Case 3: Effect of Time-Varying DeltaOmega(t) on Phase 
w_hb_case3 = w_hb_default;
v_mag_t_case3 = V_start_mag * exp(-w_hb_case3 * t_decay);

% Example time-varying detuning profile: ramp during decay
delta_w_start_Hz_profile = -100; % Hz
delta_w_end_Hz_profile   = 100;  % Hz
delta_w_varying_profile_rad_s = 2 * pi * (delta_w_start_Hz_profile + ...
    (delta_w_end_Hz_profile - delta_w_start_Hz_profile) * (t_decay / decay_duration));

% Integrate delta_w(t) to get phase contribution
phi_from_delta_w_t = zeros(size(t_decay));
if length(t_decay) > 1
    for i = 2:length(t_decay)
        phi_from_delta_w_t(i) = trapz(t_decay(1:i), delta_w_varying_profile_rad_s(1:i));
    end
end
v_phase_t_varying = phi_start_rad + phi_from_delta_w_t;
v_complex_t_varying = v_mag_t_case3 .* exp(1j * v_phase_t_varying);

figure('Name', 'Decay: Phase for Time-Varying DeltaOmega(t)');
sgtitle(sprintf('Decay with Time-Varying $\\Delta f(t)$ ($f_{1/2}=%.0f$Hz)',w_hb_case3/(2*pi)), 'Interpreter','latex');

subplot(2,1,1);
plot(t_decay * 1e3, delta_w_varying_profile_rad_s / (2*pi), 'k-', 'LineWidth', 1.5);
title('Example Time-Varying Detuning Profile $\Delta f(t)$','Interpreter','latex');
xlabel('Time (ms)','Interpreter','latex'); 
ylabel('$\Delta f(t)$ (Hz)','Interpreter','latex');
grid on;

subplot(2,1,2);
plot(t_decay * 1e3, unwrap(v_phase_t_varying), 'm-', 'LineWidth', 1.5); hold on;
% For comparison, plot phase if detuning were average of the ramp
avg_delta_w = mean(delta_w_varying_profile_rad_s);
v_phase_t_avg_delta_w = phi_start_rad + avg_delta_w * t_decay;
plot(t_decay * 1e3, unwrap(v_phase_t_avg_delta_w), 'c--', 'LineWidth', 1);
title('Resulting Phase Evolution $\angle v(t)$','Interpreter','latex');
xlabel('Time (ms)','Interpreter','latex'); 
ylabel('Phase (rad)','Interpreter','latex');
legend({'Actual $\angle v(t)$', 'Phase for Avg. $\Delta f$'},'Location','NorthWest','Interpreter','latex', 'FontSize', 8);
grid on;
disp('PLOT 3: Shows phase evolution for an example time-varying detuning profile (ramp).');
disp('  - Compares with phase evolution if detuning were constant at the average of the ramp.');
disp(' ');
disp('Decay Range Analysis Script Finished.');
