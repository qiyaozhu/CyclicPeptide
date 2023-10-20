% Sample the rama space to minimize energy while keeping cyclic
function Energy_24residue_parameter

addpath("/scratch/qz886/my_functions");

% rotation matrices and translation vector
T_NR = [0.5255, -0.8508, 0; 0.8508, 0.5255, 0; 0, 0, 1];
T_CaR = [0.3616, -0.9323, 0; 0.9323, 0.3616, 0; 0, 0, 1];
T_CR = [0.4415, -0.8973, 0; 0.8973, 0.4415, 0; 0, 0, 1];
R_omega = [1, 0, 0; 0, -1, 0; 0, 0, -1];
q = [3.5620; 1.3322; 0];

% Ramachandran plot for symmetric glycine
Boundaries = load('ramabin_glycine.mat').Boundaries.';
A_rama = load('rama_bicubic_interpolation.mat').A;

% energy weights
w_rep = 0.55;
w_intra_rep = 0.005;

% atom properties
n = 24; % number of residues
[LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V] = FA_parameter(n); % parameters
D = n_bonds(n); % number of bonds between any atom pair

% Simulated annealing parameters
M = 10000;
t0_rama = 30;
t0_rep = 100;
t0_cyc = 6;
t0_hbond = 6;
Repeat = 10;

% thresholds
rama_threshold = 8*n;
rep_threshold = 20;
cyc_threshold = 1;
count_threshold = ceil(n/3);

% Simulated annealing to find cyclic peptides with at least n/3 hydrogen
% bonds and low energies, espeically the repulsive VDW energy, with
% different parameters
Parameters = load('parameters_combdesign24.mat').Parameters;
GoodAngleCount = zeros(Repeat, size(Parameters,2));
AcceptCounts = zeros(Repeat, size(Parameters,2));
InitAngles = [-1.7, 1.7; 0.17, -0.17];

parfor parind = 1 : size(Parameters,2)

    fprintf('index = %d\n', parind);

    k0 = Parameters(1,parind);
    b = Parameters(2,parind);
    c_rama = Parameters(3,parind);
    c_rep = Parameters(4,parind);
    c_cyc = Parameters(5,parind);
    c_hbond = Parameters(6,parind);

    % record good candidates
    for repeat = 1 : Repeat
        % initial angles
%         bin = [1,2,1,2,1,1,1,1,1,1,2,1,1,2,1,1,2,2,2,1,2,1,2,2];
        bin = ones(1,n);
        angles = [InitAngles(:,bin(1)); InitAngles(:,bin(2)); InitAngles(:,bin(3)); ...
        InitAngles(:,bin(4)); InitAngles(:,bin(5)); InitAngles(:,bin(6)); InitAngles(:,bin(7)); ...
        InitAngles(:,bin(8)); InitAngles(:,bin(9)); InitAngles(:,bin(10)); InitAngles(:,bin(11)); ...
        InitAngles(:,bin(12)); InitAngles(:,bin(13)); InitAngles(:,bin(14)); InitAngles(:,bin(15)); ...
        InitAngles(:,bin(16)); InitAngles(:,bin(17)); InitAngles(:,bin(18)); ...
        InitAngles(:,bin(19)); InitAngles(:,bin(20)); InitAngles(:,bin(21)); InitAngles(:,bin(22)); ...
        InitAngles(:,bin(23)); InitAngles(:,bin(24))];
        coordinates = peptide(angles);
        rama = E_rama(angles, A_rama);
        [out_rep, intra_rep] = fa_rep(coordinates, D, LJ_radius, LJ_well);
        rep = out_rep*w_rep + intra_rep*w_intra_rep;
        cyc = f_cyc24(angles, R_omega, T_CR, T_NR, T_CaR, q);
        [hbond, ~] = E_hbond(coordinates);

        i = 1;
        goodcount = 0;
        acceptcount = 0;

        % stop when angles outside of the Rama space or reached maximum number of iterations
        while i <= M
            % generate new random move
            k = k0/(1+b*i/M);
            p = random_move24(k);
            [~, angles_new] = feasible24(p, angles, Boundaries);

            % first layer: check rama plot
            rama_new = E_rama(angles_new, A_rama);
            rama_explore = false;
            temp_rama = t0_rama/(1+c_rama*i/M);
            if rama_new <= rama || rama_new <= rama_threshold
                rama_explore = true;
            else
                prob = exp(1)^((rama-rama_new)/temp_rama);
                if rand <= prob
                    rama_explore = true;
                end
            end

            % second layer: check repulsive energy
            if rama_explore
                coordinates = peptide(angles_new);
                [out_rep, intra_rep] = fa_rep(coordinates, D, LJ_radius, LJ_well);
                rep_new = out_rep*w_rep + intra_rep*w_intra_rep;
                rep_explore = false;
                temp_rep = t0_rep/(1+c_rep*i/M);
                if rep_new <= rep || rep_new <= rep_threshold
                    rep_explore = true;
                else
                    prob = exp(1)^((rep-rep_new)/temp_rep);
                    if rand <= prob
                        rep_explore = true;
                    end
                end

                % third layer: check cyclic requirement
                if rep_explore
                    cyc_new = f_cyc24(angles_new, R_omega, T_CR, T_NR, T_CaR, q);
                    cyc_explore = false;
                    temp_cyc = t0_cyc/(1+c_cyc*i/M);
                    if cyc_new <= cyc || cyc_new <= cyc_threshold
                        cyc_explore = true;
                    else
                        prob = exp(1)^((cyc-cyc_new)/temp_cyc);
                        if rand <= prob
                            cyc_explore = true;
                        end
                    end

                    % fourth layer: check hydrogen bonds
                    if cyc_explore
                        [hbond_new, count_new] = E_hbond(coordinates);
                        hbond_explore = false;
                        temp_hbond = t0_hbond/(1+c_hbond*i/M);
                        if hbond_new <= hbond || count_new >= count_threshold
                            hbond_explore = true;
                        else
                            prob = exp(1)^((hbond-hbond_new)/temp_hbond);
                            if rand <= prob
                                hbond_explore = true;
                            end
                        end

                        % pass all layers, this random move is accepted
                        if hbond_explore
                            angles = angles_new;
                            rama = rama_new;
                            rep = rep_new;
                            cyc = cyc_new;
                            hbond = hbond_new;
                            count = count_new;
                            acceptcount = acceptcount + 1;

                            if rep <= 15 && cyc <= cyc_threshold && count >= ceil(n/3)
                                goodcount = goodcount + 1;
                            end
                        end
                    end
                end
            end
            i = i+1;
        end
        GoodAngleCount(repeat,parind) = goodcount;
        AcceptCounts(repeat,parind) = acceptcount;
    end
end

save('GoodCount24.mat', 'GoodAngleCount', 'AcceptCounts');
end
