% Example: compute L for a ferrite rod inductor

N       = 283;        % turns
R_core  = 0.0325;       % 1 cm radius  [m]
L_coil  = 0.29715;       % cm length [m]
mu_r    = 1;       % guess effective mu_r 
Nz      = 40000000;       % more points -> more accurate

[L, K_N, z, Bz] = ferrite_solenoid_L(N, R_core, L_coil, mu_r, Nz);

fprintf('Inductance L  = %.6f H (%.3f mH)\n', L, L*1e3);
fprintf('Nagaoka-like K_N = %.6f\n', K_N);

% Plot Bz along axis (just to see the field shape)
figure;
plot(z, Bz);
xlabel('z [m]');
ylabel('B_z on axis [T]');
title('On-axis B field inside ferrite-core solenoid (I = 1 A)');
grid on;
