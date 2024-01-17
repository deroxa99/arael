function [Cnm_mod, Snm_mod] = gravCoeff(filename,n)
% -------------------------------------------------------------------------
% GRAVCOEFF - Read Gravitational Coefficients from file
% -------------------------------------------------------------------------
% Compute modified Gravitational Coefficients Cnm_mod and Snm_mod by 
% reading normalized coefficients data from GRAIL mission data to implement
% a spherical harmonic expansion (see gravPot).
% The coeffieicnts from the file are fully normalized and the file also 
% contains the reference values for the mean radius and the
% gravity constant
% -------------------------------------------------------------------------
% INPUTS:
% filename     : [string] - file path containing all the coefficients
% n            : [1,1] - Truncation degree for the legendre polynomials [-]
% -------------------------------------------------------------------------
% OUTPUTS:
% Cnm_mod      : [(n+2)*(n+3)/2,3] - Cnm modified gravity coefficients 
%                                    along x,y,z [-]
% Snm_mod      : [(n+2)*(n+3)/2,3] - Snm modified gravity coefficients 
%                                    along x,y,z [-]
% -------------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% -------------------------------------------------------------------------
% CHANGELOG:   2023/09/26 - V1: First draft
%              2023/10/16 - V2: Added computation of modified coefficients
%              2024/01/17 - V2: Matlab release
% -------------------------------------------------------------------------

% Define the vector length
g = (n + 1) * (n + 2) / 2;
g = int32(g);

% Get the coefficients
coeff_data = readmatrix(filename);
Cnm_norm = coeff_data(:, 3);
Snm_norm = coeff_data(:, 4);

% Pre-compute modified coefficients
Cnm_mod = zeros(g, 3);
Snm_mod = zeros(g, 3);

for i = 3:n+1
    % Define pointers
    g = int32((i + 1) * (i + 2) / 2 - i);
    g_p = int32((i) * (i + 1) / 2 - i + 1);

    % Initialize sum quantities with j = [0,1]
    an0 = (1 / sqrt(2)) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt((i - 1) * i);
    an1 = -(1 / sqrt(2)) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt(i * (i + 1));
    cn1 = (1 / 2) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt((i - 2) * (i - 1));
    dn0 = -(sqrt(2 * i - 1) / sqrt(2 * i + 1)) * (i);
    dn1 = -(sqrt(2 * i - 1) / sqrt(2 * i + 1)) * (sqrt(i^2 - 1));

    % x direction
    Cnm_mod(g, 1) = an0 * Cnm_norm(g_p + 1); % Cn0_x
    Cnm_mod(g + 1, 1) = an1 * Cnm_norm(g_p) + cn1 * Cnm_norm(g_p + 2); % Cn1_x
    Snm_mod(g + 1, 1) = cn1 * Snm_norm(g_p + 2); % Sn1_x

    % y direction
    Cnm_mod(g, 2) = an0 * Snm_norm(g_p + 1); % Cn0_y
    Cnm_mod(g + 1, 2) = cn1 * Snm_norm(g_p + 2); % Cn1_y
    Snm_mod(g + 1, 2) = an1 * Cnm_norm(g_p) - cn1 * Cnm_norm(g_p + 2); % Sn1_y

    % z direction
    Cnm_mod(g, 3) = dn0 * Cnm_norm(g_p); % Cn0_z
    Cnm_mod(g + 1, 3) = dn1 * Cnm_norm(g_p + 1); % Cn1_z
    Snm_mod(g + 1, 3) = dn1 * Snm_norm(g_p + 1); % Sn1_z

    if i > 1
        % Compute sum quantities for j = i

        % Update pointer
        g = g + i;
        g_p = g_p + i;

        % Update multipliers
        bnn = -(1 / 2) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt((i + i - 1) * (i + i));

        % x direction
        Cnm_mod(g, 1) = bnn * Cnm_norm(g_p - 1); % Cnn_x
        Snm_mod(g, 1) = bnn * Snm_norm(g_p - 1); % Snn_x

        % y direction
        Cnm_mod(g, 2) = -bnn * Snm_norm(g_p - 1); % Cnn_y
        Snm_mod(g, 2) = bnn * Cnm_norm(g_p - 1); % Snn_y

        % z direction (no update)
        if i > 2
            % Compute sum quantities for j = i-1
            % Update pointer
            g = g - 1;
            g_p = g_p - 1;

            % Update multipliers
            bnnp = -(1 / 2) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt((i + (i - 1) - 1) * (i + (i - 1)));
            dnnp = -(sqrt(2 * i - 1) / sqrt(2 * i + 1)) * (sqrt(i^2 - (i - 1)^2));

            % x direction
            Cnm_mod(g, 1) = bnnp * Cnm_norm(g_p - 1); % Cnnp_x
            Snm_mod(g, 1) = bnnp * Snm_norm(g_p - 1); % Snnp_x

            % y direction
            Cnm_mod(g, 2) = -bnnp * Snm_norm(g_p - 1); % Cnnp_y
            Snm_mod(g, 2) = bnnp * Cnm_norm(g_p - 1); % Snnp_y

            % z direction
            Cnm_mod(g, 3) = dnnp * Cnm_norm(g_p); % Cnnp_z
            Snm_mod(g, 3) = dnnp * Snm_norm(g_p); % Snnp_z

            if i > 3
                for j = 2:i-2
                    % Compute sum quantities for 1 < j < i-1
                    % Update pointer
                    g_j = g + j - i + 1;
                    g_p_j = g_p + j - i + 1;

                    % Update multipliers
                    bnm = -(1 / 2) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt((i + j - 1) * (i + j));
                    cnm = (1 / 2) * (sqrt(2 * i - 1) / sqrt(2 * i + 1)) * sqrt((i - j - 1) * (i - j));
                    dnm = -(sqrt(2 * i - 1) / sqrt(2 * i + 1)) * (sqrt(i^2 - j^2));

                    % x direction
                    Cnm_mod(g_j, 1) = bnm * Cnm_norm(g_p_j - 1) + cnm * Cnm_norm(g_p_j + 1); % Cnm_x
                    Snm_mod(g_j, 1) = bnm * Snm_norm(g_p_j - 1) + cnm * Snm_norm(g_p_j + 1); % Snm_x

                    % y direction
                    Cnm_mod(g_j, 2) = -bnm * Snm_norm(g_p_j - 1) + cnm * Snm_norm(g_p_j + 1); % Cnm_y
                    Snm_mod(g_j, 2) = bnm * Cnm_norm(g_p_j - 1) - cnm * Cnm_norm(g_p_j + 1); % Snm_y

                    % z direction
                    Cnm_mod(g_j, 3) = dnm * Cnm_norm(g_p_j); % Cnm_z
                    Snm_mod(g_j, 3) = dnm * Snm_norm(g_p_j); % Snm_z
                end
            end
        end
    end
end
end