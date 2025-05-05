% Discrete vs Continuous Step-Response Ratio |r[n]| 
clear; clc; close all;

%   Parameters from Table 3.1-1 & Sec 2.3 
f_half   = 217;                 % Hz half-bandwidth
omega_h  = 2*pi*f_half;         % rad/s
C        = 0.235e-12;           % F
Rl       = 1/(2*omega_h*C);     % Ω
ig       = 16e-3;               % A (step input)

%   Time discretization 
T     = 1e-6;                   % 1 µs sampling interval
Nmax  = 5000;                   % number of samples to plot
n     = (0:Nmax)';              % index vector
t     = n * T;                  % continuous time vector

%  Relative detunings Δω/ω_half 
rel = [0, 0.5, 1, 1.5, 2];       % Δω/ω_half values
Delta_w = rel * omega_h;        % actual detunings [rad/s]
colors  = lines(numel(rel));

%   Preallocate ratio array 
R = zeros(numel(rel), Nmax+1);

%   Compute discrete and continuous step responses 
for k = 1:numel(rel)
    Dw = Delta_w(k);
    A  = omega_h - 1i * Dw;
    
    %  Continuous‐time step response: v_c(t) 
    % steady-state: v_inf = 2*Rl*omega_h*ig / A
    v_inf = 2*Rl*omega_h*ig / A;
    v_cont = v_inf * (1 - exp(-A * t));  
    
    %  Discrete‐time response via forward Euler:
    v_disc = zeros(size(t));
    for ii = 1:Nmax
        dv = -A * v_disc(ii) + 2*Rl*omega_h*ig;
        v_disc(ii+1) = v_disc(ii) + T * dv;
    end
    
    %  Ratio sequence r[n] = v_disc[n] / v_cont[n], then abs:
    R(k,:) = abs(v_disc ./ v_cont);
end

%  Plot |r[n]| for each detuning 
figure('Position',[100 100 700 450]); hold on; grid on;
for k = 1:numel(rel)
    plot(n, R(k,:), 'LineWidth',1.5, 'Color', colors(k,:));
end
xlabel('discrete time index n','FontSize',12);
ylabel('|v_{disc}[n] / v_{cont}(nT)|','FontSize',12);
title('Ratio of Discrete to Continuous Step Response (Eq. 3.1-5)','FontSize',14);
xlim([0 Nmax]);
ylim([0.9993 1.0011]);

%  On-curve labels at different n to avoid overlap
label_ns = [500, 1500, 2500, 3500, 4500];
for k = 1:numel(rel)
    nk = label_ns(k);
    text(nk, R(k,nk), sprintf('Δω/ω_{1/2}=%.1f', rel(k)), ...
         'Color', colors(k,:), 'FontWeight','bold', ...
         'BackgroundColor','w', 'Margin',1, ...
         'HorizontalAlignment','center');
end

hold off;
