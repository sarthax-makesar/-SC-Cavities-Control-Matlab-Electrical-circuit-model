%% Power Curves 
clear; clc; close all;

%  Constants 
f_half  = 217;                     
omega_h = 2*pi*f_half;            
C       = 0.235e-12;              
Rl      = 1 / (2 * omega_h * C);  
ig      = 16e-3;                   % A

Pf = 0.5 * Rl * ig^2 / 1e3;        % Forward power 

% Detunings & time span 
detune_list = [0, 150, 250, 400];
t_end  = 5 / omega_h;
tspan  = linspace(0, t_end, 3000);
opts   = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Storage for results 
N = numel(detune_list);
Pr = zeros(N, numel(tspan));
Pc = zeros(N, numel(tspan));

%  Loop through detunings and solve ODE 
for k = 1:N
    df      = detune_list(k);
    Delta_w = 2 * pi * df;
    A       = omega_h - 1i * Delta_w;
    rhs     = @(t,v) -A*v + 2 * Rl * omega_h * ig;

    % Solve for v(t)
    sol = ode45(@(t,y)[ real(rhs(t,y(1)+1i*y(2)));
                        imag(rhs(t,y(1)+1i*y(2))) ], ...
                tspan, [0; 0], opts);
    v_t = deval(sol, tspan, 1) + 1i * deval(sol, tspan, 2);

    % Compute powers
    i_r = ig - v_t ./ Rl;
    Pr(k,:) = 0.5 * Rl * abs(i_r).^2 / 1e3;     % Reflected power [kW]
    Pc(k,:) = Pf - Pr(k,:);                    % Cavity power [kW]
end

% Plot 
figure('Position',[100 100 800 500]); hold on; grid on;

% Plot Pf (black solid line)
plot([0 t_end]*1e6, [Pf Pf], 'k-', 'LineWidth', 2);

% Colors for Pr and Pc
col_pr = [0.2 0.4 0.8];   % blue
col_pc = [0.9 0.4 0.0];   % orange

% Plot all Pr and Pc (thinner lines)
for k = 1:N
    plot(tspan*1e6, Pr(k,:), '-', 'Color', col_pr, 'LineWidth', 1.2);
    plot(tspan*1e6, Pc(k,:), '-', 'Color', col_pc, 'LineWidth', 1.2);
end

xlabel('time [\mus]', 'FontSize',12);
ylabel('Power [kW]', 'FontSize',12);
title('Pf (black), Pr (blue), Pc (orange) vs Time', 'FontSize',14);
xlim([0 t_end*1e6]);
label_times = [0.2, 0.3, 0.4, 0.55, 0.65, 0.75, 0.83];  % as fraction of total time
for k = 1:N
    idx = round(label_times(k) * numel(tspan));
    x = tspan(idx)*1e6;
    ypr = Pr(k,idx);
    ypc = Pc(k,idx);

    text(x, ypr + 8, sprintf('Pr, df=%dHz', detune_list(k)), ...
        'Color', col_pr, 'FontWeight','bold', ...
        'BackgroundColor','w','Margin',1, ...
        'HorizontalAlignment','center');

    text(x, ypc - 8, sprintf('Pc, df=%dHz', detune_list(k)), ...
        'Color', col_pc, 'FontWeight','bold', ...
        'BackgroundColor','w','Margin',1, ...
        'HorizontalAlignment','center');
end

hold off;
