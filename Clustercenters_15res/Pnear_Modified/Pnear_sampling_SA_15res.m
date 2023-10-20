function Pnear_sampling_SA_15res(seq, name)

addpath("/scratch/qz886/my_functions");
n = 15;
SampleSize = 1000;
seq = strsplit(seq, '-');

% Get corresponding Rama information for the given sequence
RamaMap = zeros(37,37,n);
RamaEnergy = zeros(37,37,n);
for res = 1 : n
    RamaMap(:,:,res) = load("/scratch/qz886/RamaMap/ramamap_"+seq(res)+".mat").Map;
    amiaci = char(seq(res));
    if amiaci(1) == 'd'
        RamaEnergy(:,:,res) = load("/scratch/qz886/RamaMap/ramamap_"+seq(res)+".mat").E;
    else
        RamaEnergy(:,:,res) = load("/scratch/qz886/RamaMap/ramamap_"+seq(res)+".mat").F;
    end
end

% Get native designed structure coordinates
file = "/scratch/qz886/Clustercenters_15res/Pnear_Candidates/output_"+name+".pdb";
backbone = pdb2mat(file);
coordinate = [backbone.X.', backbone.Y.', backbone.Z.'];
atomName = backbone.atomName;
native_backbone = zeros(4*n,3);

N_coor = coordinate(atomName=="N",:);
Ca_coor = coordinate(atomName=="CA",:);
C_coor = coordinate(atomName=="C",:);
O_coor = coordinate(atomName=="O",:);

for res = 1 : n
    native_backbone(4*res-3,:) = N_coor(res,:);
    native_backbone(4*res-2,:) = Ca_coor(res,:);
    native_backbone(4*res-1,:) = C_coor(res,:);
    native_backbone(4*res,:) = O_coor(res,:);
end

% Simulated annealing on the filtered backbones
fprintf('Running simulated annealing ...\n');
candAngles = [];
candScores = [];

% Initial points for minimizing energies and RMSD
points = InitialPoints(seq, SampleSize);
parfor cand = 1 : size(points,2)
    angles = points(:,cand);
    [Angles, Scores] = SA(angles, RamaMap, RamaEnergy);
    candAngles = [candAngles, Angles];
    candScores = [candScores, Scores];
    fprintf('Energy minimizing SA center %d, number of candidates obtained %d\n', cand, size(Angles,2));
end
parfor cand = 1 : size(points,2)
    angles = points(:,cand);
    [Angles, RMSD] = SA_RMSD(angles, RamaMap, native_backbone);
    candAngles = [candAngles, Angles];
    candScores = [candScores, 40*ones(1,length(RMSD))];
    fprintf('RMSD minimizing SA center %d, number of candidates obtained %d\n', cand, size(Angles,2));
end

% Write into file the good candidates
candPhi = zeros(n,length(candScores));
candPsi = zeros(n,length(candScores));
candOmega = 180*ones(n,length(candScores));
candCoordinates_x = zeros(n,length(candScores));
candCoordinates_y = zeros(n,length(candScores));
candCoordinates_z = zeros(n,length(candScores));
for cand = 1 : length(candScores)
    angles = candAngles(:,cand);
    coordinates = peptide(angles);
    candPhi(:,cand) = angles(1:2:2*n-1);
    candPsi(:,cand) = angles(2:2:2*n);
    candCoordinates_x(:,cand) = coordinates(3:7:7*n,1);
    candCoordinates_y(:,cand) = coordinates(3:7:7*n,2);
    candCoordinates_z(:,cand) = coordinates(3:7:7*n,3);
end

candPhi = rad2deg(candPhi);
candPsi = rad2deg(candPsi);
filename = sprintf('After_SA_%s_15res.mat', name);
save(filename, 'candPhi', 'candPsi', 'candOmega', 'candCoordinates_x', 'candCoordinates_y', 'candCoordinates_z', 'candScores');
fprintf('Number of candidates obtained %d\n', size(candAngles,2));
end


function [candAngles, candScores] = SA(angles, RamaMap, RamaEnergy)

% parameters
n = 15;
[LJ_radius, LJ_well, ~, ~, ~, ~] = FA_parameter(n); 
D = n_bonds(n);
w_rep = 0.55;
w_intra_rep = 0.005;

M = 10000;
t0_rama = 20;
t0_rep = 60;
t0_cyc = 4;
t0_hbond = 4;

k0 = 0.7;
b = 16;
c_rama = 2;
c_rep = 14;
c_cyc = 16;
c_hbond = 18;

% thresholds
rama_threshold = 8*n;
rep_threshold = 15;
cyc_threshold = 0.3;
count_threshold = ceil(n/3);

% rotation matrices and translation vector
T_NR = [0.5255, -0.8508, 0; 0.8508, 0.5255, 0; 0, 0, 1];
T_CaR = [0.3616, -0.9323, 0; 0.9323, 0.3616, 0; 0, 0, 1];
T_CR = [0.4415, -0.8973, 0; 0.8973, 0.4415, 0; 0, 0, 1];
R_omega = [1, 0, 0; 0, -1, 0; 0, 0, -1];
q = [3.5620; 1.3322; 0];

% calculate initial energies
coordinates = peptide(angles);
rama = rama_calc(angles, RamaEnergy);
[out_rep, intra_rep] = fa_rep(coordinates, D, LJ_radius, LJ_well);
rep = out_rep*w_rep + intra_rep*w_intra_rep;
cyc = f_cyc15(angles, R_omega, T_CR, T_NR, T_CaR, q);
[hbond, count] = E_hbond(coordinates);

if rep <= 10 && cyc <= 0.3 && count >= ceil(n/3)
    score = 0.45*rama + rep + hbond;
    candAngles = angles;
    candScores = score;
else
    candAngles = [];
    candScores = [];
end

i = 1;
% stop when angles outside of the Rama space or reached maximum number of iterations
while i <= M
    % generate new random move
    k = k0/(1+b*i/M);
    p = random_move15(k);
    angles_new = feasible_move(p, angles, RamaMap);

    % first layer: check rama plot
    rama_new = rama_calc(angles_new, RamaEnergy);
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

                % this random move is accepted
                if hbond_explore
                    angles = angles_new;
                    rama = rama_new;
                    rep = rep_new;
                    cyc = cyc_new;
                    hbond = hbond_new;
                    count = count_new;

                    if rep <= 10 && cyc <= 0.3 && count >= ceil(n/3)
                        score = 0.45*rama + rep + hbond;
                        candAngles = [candAngles, angles];
                        candScores = [candScores, score];
                    end
                end
            end
        end
    end
    i = i+1;
end
end


function [candAngles, candRMSD] = SA_RMSD(angles, RamaMap, native_backbone)
% parameters
M = 10000;
t0_cyc = 4;
t0_rmsd = 1;

k0 = 0.7;
b = 16;
c_cyc = 16;
c_rmsd = 50;

cyc_threshold = 0.3;
rmsd_threshold = 1.5;

% rotation matrices and translation vector
T_NR = [0.5255, -0.8508, 0; 0.8508, 0.5255, 0; 0, 0, 1];
T_CaR = [0.3616, -0.9323, 0; 0.9323, 0.3616, 0; 0, 0, 1];
T_CR = [0.4415, -0.8973, 0; 0.8973, 0.4415, 0; 0, 0, 1];
R_omega = [1, 0, 0; 0, -1, 0; 0, 0, -1];
q = [3.5620; 1.3322; 0];

% calculate initial energies
coordinates = get_coordinates(angles);
cyc = f_cyc15(angles, R_omega, T_CR, T_NR, T_CaR, q);
[~, ~, rmsd] = kabsch_algorithm(coordinates, native_backbone);

if cyc <= cyc_threshold && rmsd <= rmsd_threshold
    candAngles = angles;
    candRMSD = rmsd;
else
    candAngles = [];
    candRMSD = [];
end

i = 1;
% stop when angles outside of the Rama space or reached maximum number of iterations
while i <= M
    % generate new random move
    k = k0/(1+b*i/M);
    p = random_move15(k);
    angles_new = feasible_move(p, angles, RamaMap);

    % check cyclic requirement
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

    % check rmsd
    if cyc_explore
        coordinates = get_coordinates(angles_new);
        [~, ~, rmsd_new] = kabsch_algorithm(coordinates, native_backbone);
        rmsd_explore = false;
        temp_rmsd = t0_rmsd/(1+c_rmsd*i/M);
        if rmsd_new <= rmsd || rmsd_new <= rmsd_threshold
            rmsd_explore = true;
        else
            prob = exp(1)^((rmsd-rmsd_new)/temp_rmsd);
            if rand <= prob
                rmsd_explore = true;
            end
        end

        % this random move is accepted
        if rmsd_explore
            angles = angles_new;
            cyc = cyc_new;
            rmsd = rmsd_new;

            if cyc <= cyc_threshold && rmsd <= rmsd_threshold
                candAngles = [candAngles, angles];
                candRMSD = [candRMSD, rmsd];
            end
        end
    end
    i = i+1;
end
end


function coordinates = get_coordinates(angles)
% bond angles
NR = deg2rad(121.7);
CaR = deg2rad(111.2);
CR = deg2rad(116.2);
COR = deg2rad(123.0);

% bond lengths
NCa = 1.458;
CaC = 1.524;
CN = 1.329;
CO = 1.231;

n = length(angles)/2;
coordinates = zeros(4*n, 3);
N = [0;0;0];
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n
    rotation = rotation*T(NR)*R(angles(2*i-1));
    Ca = N + rotation*[NCa; 0; 0];
    rotation = rotation*T(CaR)*R(angles(2*i));
    C = Ca + rotation*[CaC; 0; 0];
    rotation = rotation*T(CR)*R(pi);
    O = C + rotation*[-CO*cos(pi-COR); CO*sin(pi-COR); 0];
    N = C + rotation*[CN; 0; 0];

    coordinates(4*i-2,:) = Ca.';
    coordinates(4*i-1,:) = C.';
    coordinates(4*i,:) = O.';
    if i < n
        coordinates(4*i+1,:) = N.';
    end
end
end


function [R, t, rmsd] = kabsch_algorithm(P, Q)

% Calculate the centroids of the two sets of points
centroid_P = mean(P, 1);
centroid_Q = mean(Q, 1);

% Center the points by subtracting their centroids
P_centered = P - centroid_P;
Q_centered = Q - centroid_Q;

% Calculate the covariance matrix of the centered points
covariance_matrix = P_centered' * Q_centered;

% Calculate the optimal rotation matrix using singular value decomposition (SVD)
[U, ~, V] = svd(covariance_matrix);
rotation_matrix = V * U';

% If the determinant of the rotation matrix is negative, we need to flip one axis
if det(rotation_matrix) < 0
    V(:, 3) = -V(:, 3);
    rotation_matrix = V * U';
end

% Calculate the translation vector
translation_vector = centroid_Q' - rotation_matrix * centroid_P';

% Apply the rotation and translation to the original set of points
P_aligned = (rotation_matrix * P')' + translation_vector';

% Calculate the root-mean-square deviation (RMSD) between the aligned points
rmsd = sqrt(sum(sum((Q - P_aligned).^2)) / size(P, 1));

% Output the rotation matrix, translation vector, and RMSD
R = rotation_matrix;
t = translation_vector;
end


function chi = torsion(p1, p2, p3, p4)
b1 = p2-p1;
b2 = p3-p2;
b3 = p4-p3;
n1 = cross(b1, b2)/norm(cross(b1, b2));
n2 = cross(b2, b3)/norm(cross(b2, b3));
m = cross(n1, b2/norm(b2));
x = dot(n1, n2);
y = dot(m, n2);
chi = atan2(y, x);
end


function rama = rama_calc(angles, RamaEnergy)
locations = floor(rad2deg(angles)/10) + 19;
n = length(angles)/2;
rama = 0;
for res = 1 : n
    rama = rama + RamaEnergy(locations(2*res-1), locations(2*res), res);
end
end


function angles_new = feasible_move(p, angles, RamaMap)
angles_new = angles + p;
angles_new(angles_new>pi) = angles_new(angles_new>pi)-2*pi;
angles_new(angles_new<-pi) = angles_new(angles_new<-pi)+2*pi;

locations = floor(rad2deg(angles_new)/10) + 19;
n = length(angles_new)/2;
for res = 1 : n
    if RamaMap(locations(2*res-1), locations(2*res), res) == 0
        angles_new(2*res-1:2*res) = angles(2*res-1:2*res);
    end
end
end


function points = InitialPoints(seq, SampleSize)
n = length(seq);
centers = deg2rad([-110, 60, -80, 110, -60, 80; ...
    90, 30, 60, -90, -30, -60]);

for res = 1 : n
    amiaci = char(seq(res));
    if amiaci(1) == 'd'
        if seq(res) == "dpro"
            IP{res} = 6;
        elseif seq(res) == "dile" || seq(res) == "dthr" || seq(res) == "dval"
            IP{res} = 4;
        else
            IP{res} = [4, 5];
        end
    else
        if seq(res) == "pro"
            IP{res} = 3;
        elseif seq(res) == "ile" || seq(res) == "thr" || seq(res) == "val"
            IP{res} = 1;
        else
            IP{res} = [1, 2];
        end
    end
end

[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15]=ndgrid(IP{1},IP{2},IP{3},IP{4},IP{5},IP{6},IP{7},IP{8},IP{9},IP{10},IP{11},IP{12},IP{13},IP{14},IP{15});
bins = [x1(:),x2(:),x3(:),x4(:),x5(:),x6(:),x7(:),x8(:),x9(:),x10(:),x11(:),x12(:),x13(:),x14(:),x15(:)];

if size(bins, 1) < SampleSize
    points = zeros(n*2,size(bins, 1));
    for pos = 1 : size(bins, 1)
        bin = bins(pos,:);
        points(:,pos) = [centers(:,bin(1)); centers(:,bin(2)); centers(:,bin(3)); ...
            centers(:,bin(4)); centers(:,bin(5)); centers(:,bin(6)); centers(:,bin(7)); ...
            centers(:,bin(8)); centers(:,bin(9)); centers(:,bin(10)); centers(:,bin(11)); ...
            centers(:,bin(12)); centers(:,bin(13)); centers(:,bin(14)); centers(:,bin(15))];
    end
else
    points = zeros(n*2,SampleSize);
    samples = randperm(size(bins, 1), SampleSize);
    for pos = 1 : SampleSize
        bin = bins(samples(pos),:);
        points(:,pos) = [centers(:,bin(1)); centers(:,bin(2)); centers(:,bin(3)); ...
            centers(:,bin(4)); centers(:,bin(5)); centers(:,bin(6)); centers(:,bin(7)); ...
            centers(:,bin(8)); centers(:,bin(9)); centers(:,bin(10)); centers(:,bin(11)); ...
            centers(:,bin(12)); centers(:,bin(13)); centers(:,bin(14)); centers(:,bin(15))];
    end
end
end


function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end