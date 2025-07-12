function [A, u_sc] = generate_system_matrix(kb, rho_s, x_range, y_range, h, rx, ry, chi)
%GENERATE_SYSTEM_MATRIX Generates the system matrix A and scattered field vector u_sc.
%   Inputs:
%       kb          - Background wavenumber
%       rho_s       - Source position as a 2x1 vector [x; y]
%       x_range     - x range for object domain [x_min, x_max]
%       y_range     - y range for object domain [y_min, y_max]
%       h           - Grid spacing
%       rx, ry      - Receiver positions (1xM vectors)
%       x_chi       - True object contrast vector (Nx1)
%   Outputs:
%       A           - System matrix of size MxN
%       u_sc        - Scattered field (Mx1)



    % Grid setup for D_obj
    x = x_range(1):h:x_range(2);
    y = y_range(1):h:y_range(2);
    [X, Y] = meshgrid(x, y);
    [n_y, n_x] = size(X);
    N = n_x * n_y;                  % Number of elements in discrete D_obj
    M = length(rx);                 % Number of receivers

    % Distance from grid points to source
    R = sqrt((X - rho_s(1)).^2 + (Y - rho_s(2)).^2);
    u_inc = -1i/4 * besselh(0, 2, kb * R);  % Incident field at D_obj

    x_chi = reshape(chi, [N, 1]);  % Vectorized sampled contrast

    % Flatten grid
    x_grid = reshape(X, [N, 1]);
    y_grid = reshape(Y, [N, 1]);
    u_inc_vec = reshape(u_inc, [N, 1]);

    % Build system matrix A
    A = zeros(M, N);
    for m = 1:M
        r_m = [rx(m), ry(m)];
        dists = sqrt((x_grid - r_m(1)).^2 + (y_grid - r_m(2)).^2);
        G = -1i/4 * besselh(0, 2, kb * dists);
        A(m, :) = (kb^2 * h^2) * G.' .* u_inc_vec.';
    end

    % Simulate scattered field
    u_sc = A * x_chi;
end
