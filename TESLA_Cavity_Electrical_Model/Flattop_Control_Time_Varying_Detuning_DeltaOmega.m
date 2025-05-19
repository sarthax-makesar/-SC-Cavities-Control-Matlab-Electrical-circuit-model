% MATLAB Script for Flattop Control with Time-Varying Detuning DeltaOmega(t)


clearvars;
close all;
clc;

%  Parameters 
f_hb = 217; % Cavity half-bandwidth (Hz)
w_hb = 2 * pi * f_hb; % omega_1/2 (rad/s)

Vc_mag = 25e6; % Target cavity voltage magnitude |Vc| (V)
Vc_phase_rad = 0; % Target cavity voltage phase (rad)
Vc_complex = Vc_mag * exp(1j * Vc_phase_rad); % Complex Vc

Ib_mag = 8e-3; % Beam current magnitude |Ib| (A)
phi_b_rel_Vc = 0; % Relative beam phase to Vc (rad, 0 for on-crest)
Ub_phase_rad = Vc_phase_rad + phi_b_rel_Vc; % Absolute phase of beam current
w0_rho = 4.247e12; % Cavity parameter omega0*rho (Ohm/s or V/(A*s))
Ub_complex = w0_rho * Ib_mag * exp(1j * Ub_phase_rad); % Complex unified beam loading (constant)

R_L_cav = w0_rho / (2 * w_hb); % Loaded shunt impedance (Ohms)

flattop_duration = 0.8e-3; % Flattop duration (s)
num_points = 800; % Number of points for time vector
t_vec = linspace(0, flattop_duration, num_points);
dt = t_vec(2)-t_vec(1);

disp(' Flattop Control with Time-Varying Detuning DeltaOmega(t) ');
disp(['Target Vc = ' num2str(Vc_mag/1e6) ' MV at phase ' num2str(Vc_phase_rad) ' rad.']);
disp(['Beam Ib = ' num2str(Ib_mag*1e3) ' mA at relative phase ' num2str(phi_b_rel_Vc) ' rad to Vc.']);
disp(['Flattop duration = ' num2str(flattop_duration*1e3) ' ms.']);
disp(' ');

%  Define DeltaOmega(t) Profiles 
% Profile 1: Constant Detuning
delta_w_const_Hz = 100; % Example: 100 Hz detuning
delta_w_profile1 = 2 * pi * delta_w_const_Hz * ones(size(t_vec));
profile1_name = sprintf('Constant Detuning (%.0f Hz)', delta_w_const_Hz);

% Profile 2: Linear Ramp Detuning
delta_w_start_Hz = -200; % Hz
delta_w_end_Hz   = 200;  % Hz
delta_w_profile2 = 2 * pi * (delta_w_start_Hz + (delta_w_end_Hz - delta_w_start_Hz) * (t_vec / flattop_duration));
profile2_name = sprintf('Ramp Detuning (%.0f to %.0f Hz)', delta_w_start_Hz, delta_w_end_Hz);

% Profile 3: Sinusoidal (Microphonic-like) Detuning
delta_w_offset_Hz = 50;   % Hz
Am_micro_Hz       = 150;  % Amplitude of microphonic detuning in Hz
fm_micro          = 250;  % Frequency of microphonics in Hz 
phi_m_micro       = 0;    % Initial phase of microphonics
delta_w_profile3 = 2 * pi * (delta_w_offset_Hz + Am_micro_Hz * sin(2 * pi * fm_micro * t_vec + phi_m_micro));
profile3_name = sprintf('Sinusoidal Detuning (Offset %.0fHz, Amp %.0fHz, Freq %.0fHz)', delta_w_offset_Hz, Am_micro_Hz, fm_micro);

detuning_profiles = {delta_w_profile1, delta_w_profile2, delta_w_profile3};
profile_names     = {profile1_name, profile2_name, profile3_name};

%  Loop through profiles and plot 
for k_profile = 1:length(detuning_profiles)
    current_delta_w_profile = detuning_profiles{k_profile};
    current_profile_name    = profile_names{k_profile};
    
    % Arrays to store results for current profile
    Ac_complex_t  = zeros(1, num_points);
    u_g_complex_t = zeros(1, num_points);
    u_g_abs_t     = zeros(1, num_points);
    u_g_phase_t   = zeros(1, num_points);
    i_g_abs_t     = zeros(1, num_points);
    Pf_t          = zeros(1, num_points);
    Pr_t          = zeros(1, num_points); % Optional: Reflected Power
    
    for i = 1:num_points
        delta_w_inst = current_delta_w_profile(i); % Instantaneous detuning
        
        Ac_complex_t(i)  = -w_hb + 1j * delta_w_inst;
        u_g_complex_t(i) = Ub_complex - Ac_complex_t(i) * Vc_complex;
        
        u_g_abs_t(i)   = abs(u_g_complex_t(i));
        u_g_phase_t(i) = angle(u_g_complex_t(i));
        
        i_g_abs_t(i)   = u_g_abs_t(i) / w0_rho;
        Pf_t(i)        = 0.5 * R_L_cav * (i_g_abs_t(i)^2);
        
        % Optional: Reflected Power
        % Vg_complex_t(i) = R_L_cav * (u_g_complex_t(i) / w0_rho); % This is not Vg, this is R_L * I_g
        % If u_g is drive voltage onto a matched load representing cavity input:
        % Pr_t(i)        = (abs(Vc_complex - Vg_complex_t(i)).^2) / (2 * R_L_cav); % Simplified
        % A more common representation for Pr involves forward and cavity voltage:
        % V_forward_equivalent = u_g_complex_t(i) / (2*w_hb); % Assuming u_g is from definition u_g = 2*w_1/2*V_fwd_at_cavity_transformed
      
        
        Pr_t(i) = (abs(Vc_mag - R_L_cav * i_g_abs_t(i)).^2) / (2 * R_L_cav); % Highly simplified
    end
    
    % Calculation of average Pf and efficiency
    Pb_const = Vc_mag * Ib_mag * cos(phi_b_rel_Vc);
    Pf_avg = mean(Pf_t);
    eta_b_avg = 0;
    if Pf_avg > 1e-9 % Avoid division by zero if no power
        eta_b_avg = Pb_const / Pf_avg;
    end
    if eta_b_avg < 0, eta_b_avg=0; end
    if eta_b_avg > 1, eta_b_avg=1; end


    % Plotting for the current profile
    figure('Name', ['Flattop Analysis: ' current_profile_name]);
    sgtitle({'Flattop Control with Time-Varying Detuning $\Delta\omega(t)$'; current_profile_name}, 'Interpreter', 'latex');
    
    % Subplot 1: Detuning profile
    subplot(2,2,1);
    plot(t_vec * 1e3, current_delta_w_profile / (2*pi), 'k-', 'LineWidth', 1.5);
    title('Detuning Profile $\Delta f(t)$','Interpreter','latex');
    xlabel('Time (ms)','Interpreter','latex');
    ylabel('$\Delta f(t)$ (Hz)','Interpreter','latex');
    grid on;

    % Subplot 2: Required |u_g(t)|
    subplot(2,2,2);
    plot(t_vec * 1e3, u_g_abs_t / 1e9, 'b-', 'LineWidth', 1.5); % Show in GV/s for scale
    title('Magnitude of Req. Generator Input $|u_g(t)|$','Interpreter','latex');
    xlabel('Time (ms)','Interpreter','latex');
    ylabel('$|u_g(t)|$ (GV/s)','Interpreter','latex');
    grid on;

    % Subplot 3: Required Phase of u_g(t)
    subplot(2,2,3);
    plot(t_vec * 1e3, unwrap(u_g_phase_t), 'm-', 'LineWidth', 1.5);
    title('Phase of Req. Generator Input $\angle u_g(t)$','Interpreter','latex');
    xlabel('Time (ms)','Interpreter','latex');
    ylabel('$\angle u_g(t)$ (rad)','Interpreter','latex');
    grid on;

    % Subplot 4: Forward Power P_f(t)
    subplot(2,2,4);
    plot(t_vec * 1e3, Pf_t / 1e3, 'r-', 'LineWidth', 1.5);
    title(sprintf('Forward Power $P_f(t)$ (Avg: %.1fkW, Peak: %.1fkW)', Pf_avg/1e3, max(Pf_t)/1e3),'Interpreter','latex');
    xlabel('Time (ms)','Interpreter','latex');
    ylabel('$P_f(t)$ (kW)','Interpreter','latex');
    grid on;
    text(flattop_duration*1e3*0.1, max(Pf_t)*0.9/1e3, sprintf('Avg. $\\eta_b \\approx %.3f$', eta_b_avg), 'Interpreter','latex','FontSize',8);

    disp([' Results for Profile: ' current_profile_name ' ']);
    disp(['  Average Forward Power Pf_avg: ' num2str(Pf_avg/1e3, '%.1f') ' kW']);
    disp(['  Peak Forward Power Pf_peak: ' num2str(max(Pf_t)/1e3, '%.1f') ' kW']);
    disp(['  Beam Power Pb: ' num2str(Pb_const/1e3, '%.1f') ' kW']);
    disp(['  Average Flattop Efficiency eta_b_avg: ' num2str(eta_b_avg, '%.3f')]);
    disp(' ');
end
