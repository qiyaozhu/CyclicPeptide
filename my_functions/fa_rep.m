% Calculate Lennard-Jones repulsive
function [fa_rep, fa_intra_rep] = fa_rep(coordinates, D, LJ_radius, LJ_well)

fa_rep = 0;
fa_intra_rep = 0;
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

            % compute distance d_ij and fa_rep
            d = norm(coordinates(i,:) - coordinates(j,:));
            s_i = LJ_radius(i);
            e_i = LJ_well(i);
            s_j = LJ_radius(j);
            e_j = LJ_well(j);

            epsilon = sqrt(e_i*e_j);
            sigma = s_i + s_j;

            if d <= 0.6*sigma
                m = 20*epsilon/sigma * (-(5/3)^12 + (5/3)^6);
                b = epsilon * (13*(5/3)^12 - 14*(5/3)^6 + 1);
                E_fa_rep = m*d+b;
            elseif d <= sigma
                E_fa_rep = epsilon*((sigma/d)^12 - 2*(sigma/d)^6 + 1);
            else
                E_fa_rep = 0;
            end

            % check if atoms from different residues
            if ceil(i/7) ~= ceil(j/7)
                fa_rep = fa_rep + w*E_fa_rep;
            else
                fa_intra_rep = fa_intra_rep + w*E_fa_rep;
            end
        end
    end
end
end