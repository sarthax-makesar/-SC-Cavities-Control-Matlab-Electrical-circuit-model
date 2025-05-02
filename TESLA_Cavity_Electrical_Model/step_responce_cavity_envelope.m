%% Fig. 2.3-4 Right Panel: Phasor Trajectories with θ ≥ 0 for Δf = +150 Hz
clear; clc; close all;

% --- 1. Cavity and Generator Parameters ---
f_half   = 217;                      % Hz (half-bandwidth)
omega_h  = 2*pi*f_half;              % rad/s
C        = 0.235e-12;                % F
Rl       = 1 / (2 * omega_h * C);    % Ω

ig_mag   = 16e-3;                    % A
ig_phase = 0;                        % align with real axis
ig       = ig_mag * exp(1i * ig_phase);

% --- 2. Detuning and Styles ---
detune_list = [-400, -150, -50, 0, 50, 150, 400];  % Hz
styles      = {'-','--',':','-.','-','--',':'};

% --- 3. Time Span and ODE Options ---
t_end = 5 / omega_h;
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);

% --- 4. Figure Setup ---
figure('Position',[200 200 720 600]);
hold on; grid on; axis equal;
xlabel('Re\{v(t)\} (MV)', 'FontSize',12);
ylabel('Im\{v(t)\} (MV)', 'FontSize',12);
title('Cavity Phasor Trajectories with θ≥0 at Δf = +150 Hz', 'FontSize',14);

% --- 5. Loop Over Detunings ---
for k = 1:numel(detune_list)
    df      = detune_list(k);
    Delta_w = 2*pi*df;
    
    % Use PDF sign convention: dv/dt + jΔω v + ω_h v = ...
    A   = omega_h - 1i * Delta_w;
    rhs = @(t,v) -A*v + 2*Rl*omega_h * ig;

    % Integrate
    [~, Y] = ode45(@(t,y)[real(rhs(t,y(1)+1i*y(2))); ...
                          imag(rhs(t,y(1)+1i*y(2)))], ...
                  [0 t_end], [0;0], opts);
    v_t  = Y(:,1) + 1i*Y(:,2);
    Pv_t = v_t * 1e-6;  % MV

    % Steady-state phasor
    K     = (2*Rl) / (1 + 1i * (Delta_w / omega_h));
    v_ss  = K * ig;
    Pv_ss = v_ss * 1e-6;

    % Plot trajectory and endpoint
    plot(real(Pv_t), imag(Pv_t), 'k', 'LineStyle', styles{k}, 'LineWidth',1.8);
    plot(real(Pv_ss), imag(Pv_ss), 'ko', 'MarkerFaceColor','k');
    text(real(Pv_t(end))*1.02, imag(Pv_t(end))*1.02, ...
         sprintf('\\Deltaf = %+d Hz', df), 'Color','k','FontSize',9);

    % Angle for +150 Hz, made positive
    if df == 150
        ang_rad = abs(atan2(imag(Pv_ss), real(Pv_ss)));  % absolute value
        fprintf('→ Positive angle for Δf = +150 Hz: θ = %.4f rad\n', ang_rad);
        text(real(Pv_ss)*0.5, imag(Pv_ss)*0.5, ...
             sprintf('\\psi = %.3f rad', ang_rad), ...
             'Color','r','FontSize',12,'FontWeight','bold', ...
             'HorizontalAlignment','center');
    end
end

% --- 6. Reference Generator Phasor 2R_L i_g ---
Pg = 2 * Rl * ig * 1e-6;
quiver(0, 0, real(Pg), imag(Pg), 0, 'b', 'LineWidth',2, 'MaxHeadSize',0.6);
text(real(Pg)*1.02, imag(Pg)*1.02, '2R_L i_g', 'Color','b','FontSize',11);

hold off;
