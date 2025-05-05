Required Forward Power vs Relative Detuning


% 1) Physical constants
f_half   = 217;               % Hz (half-bandwidth)
omega_h  = 2*pi*f_half;       % rad/s
C        = 0.235e-12;         % F
Rl       = 1/(2*omega_h*C);   % Ohm (shunt impedance)
vc       = 25e6;              % V (cavity voltage)

% 2) Beam currents [A]
ib_list = [0, 2e-3, 4e-3, 6e-3, 8e-3];

% 3) Normalized detuning axis: x = Δω/ω_half = tan(ψ)
x = linspace(-1, 1, 500);

% 4) Precompute cavity current term
ic = vc / (2 * Rl);  % [A]

% 5) Plot setup
figure('Position',[100 100 720 500]); hold on; grid on;
cols = lines(numel(ib_list));

for k = 1:numel(ib_list)
    ib = ib_list(k);
    
    % 5.1) Steady-state generator current i_g(ψ)
    %     = i_b + (v_c / 2R_L)*(1 + j*tanψ)
    ig = ib + ic*(1 + 1i * x);
    
    % 5.2) Forward power Pf(ψ) = 1/2 R_L |i_g|^2  [W] → [kW]
    Pf = 0.5 * Rl * abs(ig).^2 / 1e3;
    
    % 5.3) Plot
    plot(x, Pf, 'LineWidth', 1.6, 'Color', cols(k,:));
end

% 6) Labels & styling
xlabel('Relative detuning \Delta\omega / \omega_{1/2} = \tan\psi', 'FontSize',12);
ylabel('Forward power P_f [kW]', 'FontSize',12);
title('Required Forward Power vs Normalized Detuning (v_c = 25 MV)', 'FontSize',14);
xlim([-1 1]); ylim([0 max(Pf)*1.05]);

legend("i_b = 0 mA", "2 mA", "4 mA", "6 mA", "8 mA", ...
       'Location','northwest', 'Box','off');

hold off;
