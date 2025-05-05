% Electromechanical Cavity Model Simulation 
clear; clc; close all;

%  Physical & model constants from Table 
rho      = 520;                      %  shunt impedance
QL       = 3e6;                      % Loaded Q
Rl       = rho / (4 * QL);           % Correct loaded shunt impedance 
vc       = 25;                       % MV 
k_vec    = [0.1, 0.1, 0.1, 0.5];      % Hz/(MV)^2 (4 modes)
df0_Hz   = sum(k_vec .* vc^2);       % Predetuning in Hz to cancel Lorentz pull
Dw0      = 2 * pi * df0_Hz;          % rad/s

% Mechanical mode properties
f_mech   = [235, 290, 450];          % Hz
Q_mech   = [100, 100, 100];
tau4     = 0.1e-3;                   % s (1st-order mode)

% Derived current for steady-state
ig       = vc * 1e6 / (2 * Rl);       % in A

%  Discretization 
T        = 1e-6;                     % s
T_ms     = T * 1e3;
t_total  = 20e-3;                   % 20 ms
N        = round(t_total / T);
n        = (0:N)';
t_ms     = n * T * 1e3;

% Input current pulse (10 ms on, then off)
ig_n     = ig * ones(size(n));
ig_n(t_ms > 10) = 0;

%  Preallocate 
v        = zeros(size(n));           % complex voltage [V]
ampl     = zeros(size(n));           % [MV]
phase    = zeros(size(n));           % [rad]
detHz    = zeros(size(n));           % [Hz]

% Mechanical states: 3 second-order + 1 first-order
dw4      = zeros(size(n));           % mode 4 (1st order)
W2       = zeros(2, 3, N+1);          % modes 1-3: [pos; vel] per mode

% Build system matrices for each mode
Am = zeros(2,2,3);
Bm = zeros(2,3);
for m = 1:3
    wm = 2 * pi * f_mech(m);
    qm = Q_mech(m);
    Am(:,:,m) = [0 1; -wm^2, -wm/qm];
    Bm(:,m)   = [0; -2*pi * k_vec(m)];
end
b4 = 2 * pi * k_vec(4) / tau4;

%  Simulation Loop 
for k = 1:N
    % Detuning computation
    Dw = Dw0 + sum(W2(1,:,k)) + dw4(k);
    detHz(k) = Dw / (2*pi);

    % Electrical forward-Euler update
    Ae = 2 * pi * 217 - 1i * Dw;  % omega_h - j*Dw
    dv = -Ae * v(k) + 2 * Rl * 2 * pi * 217 * ig_n(k);
    v(k+1) = v(k) + T * dv;

    % Record amplitude & phase
    ampl(k)  = abs(v(k)) / 1e6;   % MV
    phase(k) = angle(v(k));

    % Voltage squared (in MV^2)
    v2 = ampl(k)^2;

    % Mechanical updates
    for m = 1:3
        x = W2(:,m,k);
        W2(:,m,k+1) = x + T * (Am(:,:,m) * x + Bm(:,m) * v2);
    end
    dw4(k+1) = dw4(k) + T * (-dw4(k) / tau4 + b4 * v2);
end

% Trim trailing point
ampl(end) = []; phase(end) = []; detHz(end) = []; t_ms(end) = [];

% Plotting 
figure('Position',[100 100 1000 600]);

subplot(2,2,1);
plot(t_ms, ampl, 'b', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('Amplitude [MV]');
title('Amplitude during pulse'); xlim([0 10]); grid on;

subplot(2,2,2);
plot(t_ms, phase, 'r', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('Phase [rad]');
title('Phase during pulse'); xlim([0 10]); grid on;

subplot(2,2,3);
plot(t_ms, detHz, 'k', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('Detuning [Hz]');
title('Detuning during pulse'); xlim([0 10]); grid on;

subplot(2,2,4);
plot(t_ms, detHz, 'k', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('Detuning [Hz]');
title('Detuning after pulse'); xlim([10 20]); grid on;

sgtitle('Electromechanical Cavity Model Response (Corrected)', 'FontSize', 14);
