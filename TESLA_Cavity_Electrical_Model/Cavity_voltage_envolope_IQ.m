%% Phasor Trajectories + Angle Calculation for Δf = +150 Hz
clear; clc; close all;

% 1) Cavity + generator params
f_half   = 217;               
omega_h  = 2*pi*f_half;       
C        = 0.235e-12;         
Rl       = 1/(2*omega_h*C);   

ig_mag   = 16e-3;             
ig_phase = 0;                 
ig       = ig_mag * exp(1i*ig_phase);

% 2) Detuning sweep (Hz) and line styles
detune_list = [-400, -150, -50, 0, 50, 150, 400];
styles      = {'-','--',':','-.','-','--',':'};

% 3) Time span and ODE options
t_end  = 5/omega_h;
opts   = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Prepare figure
figure('Position',[200 200 700 600]);
hold on; grid on; axis equal;
xlabel('Re\{v\} (MV)'); ylabel('Im\{v\} (MV)');
title('Cavity Phasor Trajectories & Angle @ \Deltaf=+150Hz','FontSize',14);

% 4) Loop over detunings
for k = 1:numel(detune_list)
    df      = detune_list(k);
    Delta_w = 2*pi*df;
    A       = omega_h + 1i*Delta_w;
    rhs     = @(t,v) -A*v + 2*Rl*omega_h * ig;
    
    [~, Y] = ode45(@(t,y)[real(rhs(t,y(1)+1i*y(2))); imag(rhs(t,y(1)+1i*y(2)))], ...
                  [0 t_end], [0;0], opts);
    v_t = Y(:,1) + 1i*Y(:,2);
    Pv_t = v_t * 1e-6;
    
    % plot trajectory
    plot(real(Pv_t), imag(Pv_t), 'k', 'LineStyle', styles{k}, 'LineWidth',1.8);
    
    % steady‐state phasor
    K     = (2*Rl)/(1 + 1i*(Delta_w/omega_h));
    v_ss  = K * ig;
    Pv_ss = v_ss * 1e-6;
    
    % mark final phasor
    plot(real(Pv_ss), imag(Pv_ss), 'ko', 'MarkerFaceColor','k');
    
    % inline label
    text(real(Pv_t(end))*1.02, imag(Pv_t(end))*1.02, ...
         sprintf('\\Deltaf=%+dHz', df), 'Color','k','FontSize',9);
    
    % if this is the +150 Hz case, compute & annotate angle
    if df == 150
        ang_rad = atan2(imag(Pv_ss), real(Pv_ss));
        ang_deg = ang_rad * (180/pi);
        fprintf('Angle for Δf = +150 Hz endpoint: %.2f°\n', ang_deg);
        
        % annotate on plot
        txt = sprintf('θ = %.1f°', ang_deg);
        text(real(Pv_ss)*0.5, imag(Pv_ss)*0.5, txt, ...
             'Color','r','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
    end
end

% common generator arrow
Pg = 2*Rl*ig * 1e-6;
quiver(0,0, real(Pg), imag(Pg), 0, 'b', 'LineWidth',2, 'MaxHeadSize',0.6);
text(real(Pg)*1.02, imag(Pg)*1.02, '2R_L i_g', 'Color','b','FontSize',11);

hold off;
