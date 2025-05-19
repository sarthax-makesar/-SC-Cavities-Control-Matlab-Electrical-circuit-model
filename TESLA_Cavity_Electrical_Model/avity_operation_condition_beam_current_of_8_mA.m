
% MATLAB: Complex‐Envelope Cavity Simulation 


clear; close all; clc;

%% 1) TIME GRID
dt     = 1e-6;             % 1 µs
t_end  = 3.2e-3;           % 3.2 ms
t      = (0:dt:t_end)';    % column vector
N      = numel(t);

% define fill + flattop windows
t_fill   = 500e-6;         % 500 µs
t_flat   = 800e-6;         % 800 µs
idx_on   = t <= (t_fill + t_flat);
idx_beam = (t > t_fill) & (t <= t_fill + t_flat);

%% 2) CAVITY & LLRF CONSTANTS
f0       = 1.3e9;             % Hz
omega0   = 2*pi*f0;           % rad/s
rho      = 520;               % Ω (shunt impedance)
f_half   = 217;               % Hz
omega_h  = 2*pi*f_half;       % rad/s
D        = rho * omega0 / 2;  % drive coefficient

%% 3) LORENTZ DETUNING CONSTANT (single‐mode inertialess)
% want Δf = 217Hz at |V|=25MV → K = 217/(25e6)^2
K_det    = 217 / (25e6)^2;     % [Hz/V^2]

%% 4) DRIVE & BEAM PHASORS
Ig       = 16e-3;    phi_g = -0.30;   % 16 mA, phase = −0.30 rad 
Ib       =  8e-3;    phi_b =  0.11;   %  8 mA, phase = +0.11 rad

i_g      = zeros(N,1);
i_b      = zeros(N,1);
i_g(idx_on)   = Ig * exp(1j*phi_g);
i_b(idx_beam) = Ib * exp(1j*phi_b);

%% 5) STATE VECTORS
V_c      = zeros(N,1);   % complex envelope [V]
phi_c    = zeros(N,1);   % cavity phase [rad]
detuneHz = zeros(N,1);   % instantaneous detuning [Hz]

%% 6) INTEGRATION (Euler)
for n = 2:N
  % 6.1) instantaneous Lorentz detuning, inertialess
  detuneHz(n) = K_det * abs(V_c(n-1))^2;
  
  % 6.2) complex‐envelope ODE
  %    dV/dt = −(ω½ + j·2π·Δf) V + D·[i_g − i_b]
  A = -(omega_h + 1j*2*pi*detuneHz(n));
  drive = i_g(n) - i_b(n);
  dV   = A * V_c(n-1) + D * drive;
  
  V_c(n)   = V_c(n-1) + dt * dV;
  phi_c(n) = angle(V_c(n));
end

%% 7) PLOTTING – exactly the four panels of Fig 3.3-1

figure('Position',[100 50 600 900]);

% Panel 1: Cavity I,Q,|Abs|
subplot(4,1,1);
plot(t*1e6, real(V_c)/1e6, 'b', ...
     t*1e6, imag(V_c)/1e6, 'r', ...
     t*1e6, abs(V_c)/1e6, 'k','LineWidth',1.2);
ylabel('Cavity Envelope [MV]');
title('Cavity Output Envelope');
legend('I','Q','|Abs|','Location','NorthWest');
grid on;

% Panel 2: Klystron I,Q,|Abs|
subplot(4,1,2);
plot(t*1e6, real(i_g)*1e3, 'b', ...
     t*1e6, imag(i_g)*1e3, 'r', ...
     t*1e6, abs(i_g)*1e3, 'k','LineWidth',1.2);
ylabel('Klystron Output [mA]');
legend('I','Q','|Abs|','Location','NorthWest');
grid on;

% Panel 3: Phase of both cavity & klystron
subplot(4,1,3);
plot(t*1e6, phi_c, 'b', 'LineWidth',1.2); hold on;
plot(t*1e6, angle(i_g), 'r--','LineWidth',1);
ylabel('Phase [rad]');
title('Phase of Cavity Field & Klystron');
legend('Cavity','Klystron','Location','NorthEast');
grid on;

% Panel 4: Detuning
subplot(4,1,4);
plot(t*1e6, detuneHz, 'LineWidth',1.2);
ylabel('Detuning [Hz]');
xlabel('Time [\mus]');
grid on;
