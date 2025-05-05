%% Cavity Envelope Amplitude vs Time 
clear; clc; close all;

% Parameters 
f_half   = 217;                   
omega_h  = 2*pi*f_half;           
C        = 0.235e-12;             
Rl       = 1 / (2 * omega_h * C); 

ig_mag   = 16e-3;                 
ig_phase = 0;                     
ig       = ig_mag * exp(1i*ig_phase);

detune_list = [0, 50, 150, 250, 400, 600, 1000];
t_end = 5/omega_h;                
tspan = linspace(0, t_end, 3000);
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Pre‐allocate & integrate as before
Vamp = zeros(numel(detune_list), numel(tspan));
for k = 1:numel(detune_list)
    df      = detune_list(k);
    Delta_w = 2*pi*df;
    A       = omega_h - 1i*Delta_w;
    rhs     = @(t,v) -A*v + 2*Rl*omega_h * ig;

    sol = ode45(@(t,y)[ real(rhs(t,y(1)+1i*y(2)));
                        imag(rhs(t,y(1)+1i*y(2))) ], ...
                tspan, [0;0], opts);
    v_t = deval(sol, tspan, 1) + 1i*deval(sol, tspan, 2);
    Vamp(k,:) = abs(v_t) * 1e-6;  
end

% Compute overall amplitude range for offset scaling
Vmin = min(Vamp(:));
Vmax = max(Vamp(:));
Vrange = Vmax - Vmin;

% --- Plot & label ---
figure('Position',[100 100 700 450]);
hold on; grid on;
colors = lines(numel(detune_list));

for k = 1:numel(detune_list)
    % plot curve
    plot(tspan*1e6, Vamp(k,:), 'LineWidth',1.8, 'Color', colors(k,:));

    % staggered time‐fraction between 0.6 and 0.9
    frac = 0.6 + 0.3*(k-1)/(numel(detune_list)-1);
    idx  = round(frac * numel(tspan));
    
    % base label position
    x_lbl = tspan(idx)*1e6;       % µs
    y_lbl = Vamp(k, idx);
    
    % offset up/down by up to ±5% of full range
    y_offset = Vrange * 0.05 * ((-1)^(k));
    
    text(x_lbl, y_lbl + y_offset, ...
         sprintf('\\Deltaf = %d Hz', detune_list(k)), ...
         'Color', colors(k,:), ...
         'FontWeight','bold', ...
         'HorizontalAlignment','center', ...
         'BackgroundColor','white', ...
         'Margin',1);
end

xlabel('time [\mus]', 'FontSize',12);
ylabel('|\!v(t)\!| [MV]', 'FontSize',12);
title('Cavity Envelope Amplitude vs Time', 'FontSize',14);
xlim([0 t_end*1e6]);
hold off;
