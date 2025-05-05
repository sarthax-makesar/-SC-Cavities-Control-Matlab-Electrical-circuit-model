
clear; clc; close all;

%  Cavity Parameters 
omega_0 = 2 * pi * 1.3e9;
omega_half = 2 * pi * 217;
C_value = 0.235e-12;
Rl = 1 / (2 * omega_half * C_value);

% Currents 
ig = 16e-3 * exp(1i * (-0.3));
ib = 8e-3 * exp(1i * (-0.11));
i0 = ig - ib;

% Detuning 
tan_psi = 0.5;
psi = atan(tan_psi);

% Scaling to MV 
scale_factor = 2 * Rl;
scale_to_MV = 1e-6;

ig_plot = scale_factor * ig * scale_to_MV;
ib_plot = scale_factor * ib * scale_to_MV;
i0_plot = scale_factor * i0 * scale_to_MV;
v_magnitude = 2 * Rl * abs(i0) * cos(psi);
v_direction = (i0/abs(i0)) * exp(1i * psi);
v_plot = v_magnitude * v_direction * scale_to_MV;

% Setup Figure 
figure;
hold on;
axis equal;
grid on;
box on;
set(gca, 'GridLineStyle', '--');

% Set plot limits
xlim([0, 60]);
ylim([-20, 20]);

xlabel('Real Axis (MV)', 'FontSize', 12);
ylabel('Imaginary Axis (MV)', 'FontSize', 12);
title('Phasor Diagram with Annotation Arrows', 'FontSize', 14);

% Save axis limits
ax = gca;
ax_pos = ax.Position;    % normalized position of the axis
ax_xlim = xlim;
ax_ylim = ylim;

% Function to map (x,y) to normalized coordinates 
map_to_norm = @(x, y) [(x-ax_xlim(1))/(ax_xlim(2)-ax_xlim(1))*ax_pos(3) + ax_pos(1), ...
                       (y-ax_ylim(1))/(ax_ylim(2)-ax_ylim(1))*ax_pos(4) + ax_pos(2)];


% ig 
p1 = map_to_norm(0, 0);
p2 = map_to_norm(real(ig_plot), imag(ig_plot));
annotation('arrow', [p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'b', 'LineWidth', 2);

% ib 
p1 = map_to_norm(0, 0);
p2 = map_to_norm(real(ib_plot), imag(ib_plot));
annotation('arrow', [p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'g', 'LineWidth', 2);

% i0 
p1 = map_to_norm(0, 0);
p2 = map_to_norm(real(i0_plot), imag(i0_plot));
annotation('arrow', [p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'm', 'LineWidth', 2);

% v 
p1 = map_to_norm(0, 0);
p2 = map_to_norm(real(v_plot), imag(v_plot));
annotation('arrow', [p1(1) p2(1)], [p1(2) p2(2)], 'Color', 'r', 'LineWidth', 2);

text(real(ig_plot)*1.05, imag(ig_plot)*1.05, '2Rl \times i_g', 'Color', 'b', 'FontSize', 12);
text(real(ib_plot)*1.05, imag(ib_plot)*1.05, '2Rl \times i_b', 'Color', 'g', 'FontSize', 12);
text(real(i0_plot)*1.05, imag(i0_plot)*1.05, '2Rl \times i_0', 'Color', 'm', 'FontSize', 12);
text(real(v_plot)*1.05, imag(v_plot)*1.05, 'Cavity Voltage v', 'Color', 'r', 'FontSize', 12);

hold off;
