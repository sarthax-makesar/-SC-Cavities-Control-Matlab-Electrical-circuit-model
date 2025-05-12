
% Electromechanical Cavity Simulation
% Aims for more stable demonstration of Lorentz detuning effects
% Using 16mA drive and adjusting mechanical Q-factors 

clear; clc; close all;

% 1) Simulation Parameters
T = 1e-6;         % Time step [s]
tpe = 10e-3;      % Pulse end time [s]
tsm = 100e-3;     % Total simulation time [s]
tvec = 0:T:tsm;
N = numel(tvec);

% 2) Electrical Model Parameters 
fhb = 217;        % Half‐bandwidth [Hz]
omhb = 2*pi*fhb;  % Cavity half-bandwidth [rad/s]
C_cav = 0.235e-12;% Cavity capacitance [F]
idr = 16e-3;      % Drive current [A] (original from Fig 3.1-3 caption)

fpd = 417;        % Predetuning [Hz]
ompd = 2*pi*fpd;  % Predetuning [rad/s]

vc = zeros(1,N);
detHz = zeros(1,N);
detHz(1) = fpd;

%% 3) Mechanical Model Parameters
fmc = [235, 290, 450];
% --- ADJUSTED PARAMETER for more damping ---
% Qmc = [30, 30, 30];  % Mechanical quality factors (originally [100,100,100])
Qmc = [100,100,100]; % Original 

Km = [0.1, 0.1, 0.1];  % Lorentz constants for resonant modes [Hz/(MV)^2]
Kinl = 0.5;            % Lorentz constant for inertialess mode [Hz/(MV)^2]
tauinl = 0.1e-3;       % Time constant for inertialess mode [s]

nrm = numel(fmc);
omm = 2*pi*fmc;

wrs = zeros(2*nrm, N);
wis = zeros(1,N);

%% 4) Time-loop (Euler integration)
for n = 1:N-1
    sum_om_res_n = sum(wrs(1:2:end, n));
    omDyn_n = sum_om_res_n + wis(n);
    omTot_n = omDyn_n + ompd;
    detHz(n) = omTot_n / (2*pi);

    ic = 0;
    if tvec(n) < tpe
        ic = idr;
    end
    dvcdt = (-omhb + 1j*omTot_n) * vc(n) + (1/C_cav)*ic;
    vc(n+1) = vc(n) + T*dvcdt;

    vMVSq_n = (abs(vc(n))*1e-6)^2;

    for m_idx = 1:nrm
        idx = 2*m_idx-1;
        om_m = omm(m_idx);
        Qm_val = Qmc(m_idx); % Using potentially modified Q
        K_val = Km(m_idx);

        force_res = -2*pi*K_val * (om_m^2) * vMVSq_n;
        d_om_dt = wrs(idx+1,n);
        d2_om_dt2 = force_res - (om_m/Qm_val)*wrs(idx+1,n) - (om_m^2)*wrs(idx,n);
        wrs(idx,  n+1) = wrs(idx,  n) + T*d_om_dt;
        wrs(idx+1,n+1) = wrs(idx+1,n) + T*d2_om_dt2;
    end

    force_inl = -2*pi*Kinl * vMVSq_n;
    d_wis_dt = (force_inl - wis(n)) / tauinl;
    wis(n+1) = wis(n) + T*d_wis_dt;
end

sum_om_res_N = sum(wrs(1:2:end, N));
omDyn_N = sum_om_res_N + wis(N);
detHz(N) = (omDyn_N + ompd) / (2*pi);

%% 5) Plots
figure('Position',[100 100 700 950]);
plot_lim_ms = 15;

h_ax1 = subplot(4,1,1);
plot(tvec*1e3, abs(vc)/1e6,'LineWidth',1.5)
xlabel('Time [ms]'); ylabel('Amp. [MV]')
title(['Amplitude (Drive Current = ', num2str(idr*1e3), ' mA, Modified Qm)'])
grid on; ylim('auto'); xlim([0 plot_lim_ms]);

h_ax2 = subplot(4,1,2);
plot(tvec*1e3, unwrap(angle(vc)),'LineWidth',1.5)
xlabel('Time [ms]'); ylabel('Phase [rad]')
title('Phase evolution')
grid on; ylim('auto'); xlim([0 plot_lim_ms]);

h_ax3 = subplot(4,1,3);
plot(tvec*1e3, detHz,'LineWidth',1.5)
xlabel('Time [ms]'); ylabel('Detuning [Hz]')
title('Total Detuning (during and after pulse)')
grid on; ylim('auto'); xlim([0 plot_lim_ms]);

h_ax4 = subplot(4,1,4);
idx_start_after = find(tvec >= tpe, 1);
plot(tvec(idx_start_after:end)*1e3, detHz(idx_start_after:end),'LineWidth',1.5)
xlabel('Time [ms]'); ylabel('Detuning [Hz]')
title('Detuning after pulse (Zoomed, extended time)')
grid on; ylim('auto'); xlim([tpe*1e3 tsm*1e3]);