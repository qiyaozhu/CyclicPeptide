function Energy_15residue

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

% atom properties
n = 15; % number of residues
[LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V] = FA_parameter(n); % parameters
D = n_bonds(n); % number of bonds between any atom pair

% energy weights
w_atr = 1;
w_elec = 1;
w_sol = 1;
w_intra_sol = 1;
w_rep = 0.55;
w_intra_rep = 0.005;

% Simulated annealing parameters
M = 10000;
t0_rama = 20;
t0_rep = 60;
t0_cyc = 4;
t0_hbond = 4;
t0_other = 40;

k0 = 0.7;
b = 16;
c_rama = 2;
c_rep = 14;
c_cyc = 16;
c_hbond = 18;
c_other = 4;
Repeat = 3;

% thresholds
rama_threshold = 8*n;
rep_threshold = 15;
cyc_threshold = 0.3;
count_threshold = ceil(n/3);
energy_other_threshold = 20;

% Different initial angles
InitAngles = [-2, -1.7, -1.9, 1.9, 1.7, 2; ...
    -2.4, 0.17, 2.6, -2.6, -0.17, 2.4];
bin_class = load('InitialPoints_15residue.mat').bin_class;

% Simulated annealing to find cyclic peptides with at least n/3 hydrogen
% bonds and low energies, espeically the repulsive VDW energy
batches = floor(size(bin_class, 1)/1000);
for bat = 0 : batches
    GoodAngles = [];
    GoodRama = [];
    GoodRep = [];
    GoodCyc = [];
    GoodHbond = [];
    GoodCount = [];
%     GoodEnergy = [];

    if bat < batches
        uplim = 1000*bat+1000;
    else
        uplim = size(bin_class, 1);
    end

    parfor pos = 1000*bat+1 : uplim


        fprintf('index = %d\n', pos);

        % record good candidates
        repeat = 1;
        goodcount = 0;

        while repeat <= Repeat && goodcount == 0
            % initial angles
            bin = bin_class(pos,:);
            angles = [InitAngles(:,bin(1)); InitAngles(:,bin(2)); InitAngles(:,bin(3)); ...
                InitAngles(:,bin(4)); InitAngles(:,bin(5)); InitAngles(:,bin(6)); InitAngles(:,bin(7)); ...
                InitAngles(:,bin(8)); InitAngles(:,bin(9)); InitAngles(:,bin(10)); InitAngles(:,bin(11)); ...
                InitAngles(:,bin(12)); InitAngles(:,bin(13)); InitAngles(:,bin(14)); InitAngles(:,bin(15))];
            coordinates = peptide(angles);
            rama = E_rama(angles, A_rama);
            [out_rep, intra_rep] = fa_rep(coordinates, D, LJ_radius, LJ_well);
            rep = out_rep*w_rep + intra_rep*w_intra_rep;
            cyc = f_cyc15(angles, R_omega, T_CR, T_NR, T_CaR, q);
            [hbond, ~] = E_hbond(coordinates);
            %         energy_other = E_other(coordinates, D, LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V, w_atr, w_elec, w_sol, w_intra_sol);

            i = 1;
            % stop when angles outside of the Rama space or reached maximum number of iterations
            while i <= M
                % generate new random move
                k = k0/(1+b*i/M);
                p = random_move15(k);
                [~, angles_new] = feasible15(p, angles, Boundaries);

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
                        cyc_new = f_cyc15(angles_new, R_omega, T_CR, T_NR, T_CaR, q);
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

                            % fifth layer: check other energies including fa_atr, f_elec, fa_sol
                            if hbond_explore
                                accept = true;
                                %                             energy_other_new = E_other(coordinates, D, LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V, w_atr, w_elec, w_sol, w_intra_sol);
                                %                             accept = false;
                                %                             temp_other = t0_other/(1+c_other*i/M);
                                %                             if energy_other_new <= energy_other || energy_other_new <= energy_other_threshold
                                %                                 accept = true;
                                %                             else
                                %                                 prob = exp(1)^((energy_other-energy_other_new)/temp_other);
                                %                                 if rand <= prob
                                %                                     accept = true;
                                %                                 end
                                %                             end

                                % pass all layers, this random move is accepted
                                if accept
                                    angles = angles_new;
                                    rama = rama_new;
                                    rep = rep_new;
                                    cyc = cyc_new;
                                    hbond = hbond_new;
                                    count = count_new;
                                    %                                 energy_other = energy_other_new;

                                    if rep <= 10 && cyc <= 0.3 && count >= ceil(n/3)
                                        GoodAngles = [GoodAngles, angles];
                                        GoodRama = [GoodRama, rama];
                                        GoodRep = [GoodRep, rep];
                                        GoodCyc = [GoodCyc, cyc];
                                        GoodHbond = [GoodHbond, hbond];
                                        GoodCount = [GoodCount, count];
                                        %                                     GoodEnergy = [GoodEnergy, energy_other];
                                        goodcount = goodcount + 1;
                                    end
                                end
                            end
                        end
                    end
                end
                i = i+1;
            end
            repeat = repeat+1;
        end
    end
    filename = sprintf('GoodAngles15res_%d.mat', bat+1);
    save(filename, 'GoodAngles', 'GoodRama', 'GoodRep', 'GoodCyc', 'GoodHbond', 'GoodCount');
end
end