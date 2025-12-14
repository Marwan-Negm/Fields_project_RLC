function [L, K_N, z, Bz] = ferrite_solenoid_L(N, R_core, L_coil, mu_r, Nz)
% FERRITE_SOLENOID_L  Numerically compute inductance of a finite solenoid
% on a ferrite rod using on-axis field averaging.
%
% Inputs:
%   N       - total number of turns
%   R_core  - core radius   [m]
%   L_coil  - coil (winding) length  [m]
%   mu_r    - relative permeability of ferrite (effective)
%   Nz      - number of z-samples for numerical integration (e.g. 2000)
%
% Outputs:
%   L       - inductance [H]
%   K_N     - Nagaoka-like coefficient (L / L_ideal)
%   z       - z positions along the axis (for plotting) [m]
%   Bz      - B_z(z) on axis inside the coil [T] for I = 1 A

    if nargin < 5
        Nz = 2000; % default resolution if not given
    end

    mu0 = 4*pi*1e-7;         % vacuum permeability [H/m]
    I   = 1;                 % 1 A excitation
    a   = R_core;            % coil radius ≈ core radius
    n   = N / L_coil;        % turns per unit length

    % z samples along the axis inside the coil
    z = linspace(-L_coil/2, +L_coil/2, Nz);

    % Precompute some constants
    pref = mu0 * mu_r * n * I / 2;

    % On-axis field Bz(z) for a finite solenoid
    z_plus  = z + L_coil/2;
    z_minus = z - L_coil/2;

    Bz = pref * ( z_plus ./ sqrt(a^2 + z_plus.^2) ...
                - z_minus ./ sqrt(a^2 + z_minus.^2) );

    % Average B over the coil length (numerical integration)
    B_avg = trapz(z, Bz) / L_coil;

    % Core cross-sectional area
    A_core = pi * R_core^2;

    % Flux per turn (approx) and inductance
    phi = B_avg * A_core;        % Wb (for I = 1 A)
    L   = N * phi / I;           % H

    % Ideal long-solenoid inductance with same mu_r
    L_ideal = mu0 * mu_r * (N^2 * A_core) / L_coil;

    % Nagaoka-like correction factor
    K_N = L / L_ideal;
end
