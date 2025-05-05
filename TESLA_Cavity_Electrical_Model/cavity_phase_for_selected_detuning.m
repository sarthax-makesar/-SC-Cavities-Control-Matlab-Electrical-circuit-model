%% Step Response of Cavity Envelope Amplitude & Phase
clear; clc; close all;

% Cavity & Generator Parameters 
f_half   = 217;                   % Hz (half-bandwidth)
omega_h  = 2*pi*f_half;           % rad/s
C        = 0.235e-12;             % F
Rl       = 1 / (2 * omega_h * C); % Ω

ig_mag   = 16e-3;                 % A
ig_phase = 0;                     % align with real axis
ig       = ig_mag * exp(1i * ig_phase);

% Detuning List in Hz 
detune_list = [0, 50, 150, 250, 400, 600, 1000];

%  Time Span & ODE Options 
t_end = 5 / omega_h;              % approx. 3.7 ms
tspan = linspace(0, t_end, 3000);
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Pre‐allocate
V  = zeros(numel(detune_list), numel(tspan));
PH = zeros(size(V));

%  Integrate ODE for detuning 
for k = 1:numel(detune_list)
    df      = detune_list(k);
    Delta_w = 2*pi*df;
    A       = omega_h - 1i*Delta_w;
    rhs     = @(t,v) -A*v + 2*Rl*omega_h * ig;
    
    sol = ode45(@(t,y)[ real(rhs(t,y(1)+1i*y(2)));
                        imag(rhs(t,y(1)+1i*y(2))) ], ...
                tspan, [0;0], opts);
    v_t = deval(sol, tspan, 1) + 1i*deval(sol, tspan, 2);
    
    V(k,:)  = abs(v_t) * 1e-6;      % MV
    PH(k,:) = unwrap(angle(v_t));
end

% compute overall ranges for offsets
Vmin = min(V(:));   Vmax = max(V(:));   Vrange = Vmax - Vmin;
Pmin = min(PH(:));  Pmax = max(PH(:));  Prange = Pmax - Pmin;

% choose a colormap
colors = lines(numel(detune_list));

%% Plot Amplitude with On‐Curve Labels 
figure('Position',[100 100 600 450]);
hold on; grid on;
for k = 1:numel(detune_list)
    plot(tspan*1e6, V(k,:), 'LineWidth',1.6, 'Color', colors(k,:));
    
    % stagger label at frac∈[0.6,0.9]
    frac = 0.6 + 0.3*(k-1)/(numel(detune_list)-1);
    idx  = round(frac * numel(tspan));
    x_lbl = tspan(idx)*1e6;           % µs
    y_lbl = V(k,idx) + Vrange*0.04*((-1)^k);  % ±4% offset
    
    text(x_lbl, y_lbl, sprintf('\\Deltaf = %d Hz', detune_list(k)), ...
         'Color', colors(k,:), 'FontWeight','bold', ...
         'HorizontalAlignment','center', ...
         'BackgroundColor','white', 'Margin',1);
end
xlabel('time [\mus]','FontSize',12);
ylabel('|\!v(t)\!| [MV]','FontSize',12);
title('Cavity Envelope Amplitude vs Time','FontSize',14);
xlim([0 t_end*1e6]);
hold off;

Plot Phase with On‐Curve Labels 
figure('Position',[100 100 600 450]);
hold on; grid on;
for k = 1:numel(detune_list)
    plot(tspan*1e6, PH(k,:), 'LineWidth',1.6, 'Color', colors(k,:));
    
    % stagger label at frac∈[0.6,0.9]
    frac = 0.6 + 0.3*(k-1)/(numel(detune_list)-1);
    idx  = round(frac * numel(tspan));
    x_lbl = tspan(idx)*1e6;           % µs
    y_lbl = PH(k,idx) + Prange*0.04*((-1)^(k+1)); % alternate offset
    
    text(x_lbl, y_lbl, sprintf('\\Deltaf = %d Hz', detune_list(k)), ...
         'Color', colors(k,:), 'FontWeight','bold', ...
         'HorizontalAlignment','center', ...
         'BackgroundColor','white', 'Margin',1);
end
xlabel('time [\mus]','FontSize',12);
ylabel('\angle v(t) [rad]','FontSize',12);
title('Cavity Envelope Phase vs Time','FontSize',14);
xlim([0 t_end*1e6]);
hold off;
