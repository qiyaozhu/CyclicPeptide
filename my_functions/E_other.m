function energy = E_other(coordinates, D, LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V, w_atr, w_elec, w_sol, w_intra_sol)
atr = fa_atr(coordinates, D, LJ_radius, LJ_well);
elec = fa_elec(coordinates, D, Coulumb);
[sol, intra_sol] = fa_sol(coordinates, D, LJ_radius, LK_DG, LK_Lambda, LK_V);
energy = w_atr*atr + w_elec*elec + w_sol*sol + w_intra_sol*intra_sol;
end


% Calculate Lennard-Jones attractive
function fa_atr = fa_atr(coordinates, D, LJ_radius, LJ_well)

fa_atr = 0;
N = size(coordinates,1); % number of atoms

% Calculate energy for all atom pairs (i,j)
for i = 1 : N-1
    for j = i+1 : N

        % only consider atoms from different residues
        if ceil(i/7) ~= ceil(j/7)

            % only consider atoms separated by at least 4 bonds
            k = D(i,j);
            if k >= 4

                % set the connectivity weight
                if k == 4
                    w = 0.2;
                else
                    w = 1;
                end

                % compute distance d_ij and fa_atr
                d = norm(coordinates(i,:) - coordinates(j,:));
                s_i = LJ_radius(i);
                e_i = LJ_well(i);
                s_j = LJ_radius(j);
                e_j = LJ_well(j);

                epsilon = sqrt(e_i*e_j);
                sigma = s_i + s_j;
                fa_atr = fa_atr + w*E_fa_atr(d, epsilon, sigma);
            end
        end
    end
end
end


function fa_atr = E_fa_atr(d, epsilon, sigma)

if d <= sigma
    fa_atr = -epsilon;
elseif d <= 4.5
    fa_atr = epsilon*((sigma/d)^12 - 2*(sigma/d)^6);
elseif d <= 6
    [C0, C1, C2, C3] = Cubic_inter(4.5, 6, epsilon*((sigma/4.5)^12 - 2*(sigma/4.5)^6), ...
        0, epsilon/4.5*(-12*(sigma/4.5)^12 + 12*(sigma/4.5)^6), 0);
    fa_atr = C0 + C1*d + C2*d^2 + C3*d^3;
else
    fa_atr = 0;
end
end


% Calculate Coulombic electrostatic potential
function fa_elec = fa_elec(coordinates, D, Coulumb)

fa_elec = 0;
N = size(coordinates,1); % number of atoms

% Calculate energy for all atom pairs (i,j)
for i = 1 : N-1
    for j = i+1 : N
        % only consider atoms separated by at least 4 bonds
        k = D(i,j);
        if k >= 4

            % set the connectivity weight
            if k == 4
                w = 0.2;
            else
                w = 1;
            end

            % compute distance d_ij and fa_elec
            d = norm(coordinates(i,:) - coordinates(j,:));
            q_i = Coulumb(i);
            q_j = Coulumb(j);
            fa_elec = fa_elec + w*E_fa_elec(d, q_i, q_j);
        end
    end
end
end


function fa_elec = E_fa_elec(d, q_i, q_j)

D = 80;
D0 = 4;
S = 0.4;
Co = 322.0637;
d_min = 1.5;
d_max = 5.5;

elf = @(x) D - (D-D0)/2 * ((x*S)^2+2*x*S+2) * exp(-x*S); % sigmoid function between D0 and D
C1 = Co*q_i*q_j / (d_max*elf(d_max));
E = @(x) Co*q_i*q_j / (x*elf(x)) - C1;
E_deriv = @(x) -Co*q_i*q_j/(x^2*elf(x)) - Co*q_i*q_j/(x*elf(x)^2) * ...
    (-(D-D0)/2*(2*x*S^2+2*S)*exp(-x*S) + (D-D0)/2*((x*S)^2+2*x*S+2)*exp(-x*S)*S);

if d < d_min-0.25
    fa_elec = E(d_min);
elseif d < d_min+0.25
    [C0, C1, C2, C3] = Cubic_inter(d_min-0.25, d_min+0.25, E(d_min), E(d_min+0.25), 0, E_deriv(d_min+0.25));
    fa_elec = C0 + C1*d + C2*d^2 + C3*d^3;
elseif d < d_max-1
    fa_elec = E(d);
elseif d < d_max
    [C0, C1, C2, C3] = Cubic_inter(d_max-1, d_max, E(d_max-1), 0, E_deriv(d_max-1), 0);
    fa_elec = C0 + C1*d + C2*d^2 + C3*d^3;
else
    fa_elec = 0;
end
end


% Calculate isotropic solvation energy
function [fa_sol, fa_intra_sol] = fa_sol(coordinates, D, LJ_radius, LK_DG, LK_Lambda, LK_V)

fa_sol = 0;
fa_intra_sol = 0;
N = size(coordinates,1); % number of atoms

% Calculate energy for all atom pairs (i,j)
for i = 1 : N-1
    for j = i+1 : N
        % only consider atoms separated by at least 4 bonds
        k = D(i,j);
        if k >= 4

            % set the connectivity weight
            if k == 4
                w = 0.2;
            else
                w = 1;
            end

            % compute distance d_ij and fa_sol
            d = norm(coordinates(i,:) - coordinates(j,:));
            s_i = LJ_radius(i);
            s_j = LJ_radius(j);
            DG_i = LK_DG(i);
            DG_j = LK_DG(j);
            lambda_i = LK_Lambda(i);
            lambda_j = LK_Lambda(j);
            V_i = LK_V(i);
            V_j = LK_V(j);

            % check if atoms from different residues
            if ceil(i/7) ~= ceil(j/7)
                fa_sol = fa_sol + w * (E_fa_sol(d, s_i, s_j, DG_i, lambda_i, V_j) + E_fa_sol(d, s_j, s_i, DG_j, lambda_j, V_i));
            else
                fa_intra_sol = fa_intra_sol + w * (E_fa_sol(d, s_i, s_j, DG_i, lambda_i, V_j) + E_fa_sol(d, s_j, s_i, DG_j, lambda_j, V_i));
            end            
        end
    end
end
end


function fa_sol = E_fa_sol(d, sigma_i, sigma_j, DG, lambda, V)

c0 = 0.3;
c1 = 0.2;
sigma = sigma_i + sigma_j;

E = @(x) -V*DG / (2*pi^1.5*lambda*x^2) * exp(-((x-sigma)/lambda)^2);
E_deriv = @(x) 2*V*DG / (2*pi^1.5*lambda*x^3) * exp(-((x-sigma)/lambda)^2) - ...
    V*DG / (2*pi^1.5*lambda*x^2) * exp(-((x-sigma)/lambda)^2) * (-2*(x-sigma)/lambda^2);

if d <= sigma-c0
    fa_sol = E(sigma);
elseif d <= sigma+c1
    [C0, C1, C2, C3] = Cubic_inter(sigma-c0, sigma+c1, E(sigma), E(sigma+c1), 0, E_deriv(sigma+c1));
    fa_sol = C0 + C1*d + C2*d^2 + C3*d^3;
elseif d <= 4.5
    fa_sol = E(d);
elseif d <= 6
    [C0, C1, C2, C3] = Cubic_inter(4.5, 6, E(4.5), 0, E_deriv(4.5), 0);
    fa_sol = C0 + C1*d + C2*d^2 + C3*d^3;
else
    fa_sol = 0;
end
end


% Cubic interpolation with two endpoints and their derivatives
% f(x) = C0 + C1*x + C2*x^2 + C3*x^3
% f(a) = A, f(b) = B, f'(a) = alpha, f'(b) = beta
function [C0, C1, C2, C3] = Cubic_inter(a, b, A, B, alpha, beta)
C0 = 1/(a-b)^3 * (3*a*A*b^2 - a^2*alpha*b^2 - A*b^3 + a*alpha*b^3 + ...
    a^3*B - 3*a^2*b*B - a^3*b*beta + a^2*b^2*beta);
C1 = 1/(a-b)^3 * (-6*a*A*b + 2*a^2*alpha*b - a*alpha*b^2 - alpha*b^3 + ...
    6*a*b*B + a^3*beta + a^2*b*beta - 2*a*b^2*beta);
C2 = 1/(a-b)^3 * (3*a*A - a^2*alpha + 3*A*b - a*alpha*b + 2*alpha*b^2 - ...
    3*a*B - 3*b*B - 2*a^2*beta + a*b*beta + b^2*beta);
C3 = 1/(a-b)^3 * (-2*A + a*alpha - alpha*b + 2*B + a*beta - b*beta);
end