function Pnear_Clustering_15res(Names)

addpath("/scratch/qz886/my_functions");
n = 15;
score_threshold = 0;

results = cell(1, length(Names));
parfor i = 1 : length(Names)
    name = Names{i};
    seconds = 0;
    tic;

    filename = sprintf('After_SA_%s_Angles_15res.txt', name);
    Peptides = readmatrix(filename).';
    Scores = Peptides(3*n+1,:);
    RMSD = Peptides(end,:);
    Angles = Peptides(1:3*n,:);
    Angles = deg2rad(Angles);

    N = size(Angles, 2);
    Coordinates = zeros(4*n, 3*N);
    for cand = 1 : N
        angle = Angles(:,cand);
        Coordinates(:,3*(cand-1)+1:3*cand) = get_coordinates(angle);
    end
    Coordinates = reshape(Coordinates, [4*n*3, N]);
    disp("A total population of "+N+" after SA.");

    % Clustering to select lowest energy centers
    [Scores, bestIndices] = mink(Scores, length(Scores));
    Angles = Angles(:, bestIndices);
    RMSD = RMSD(bestIndices);
    Coordinates = Coordinates(:, bestIndices);

    indices = clustering(Coordinates);
    Angles = Angles(:, indices);
    Scores = Scores(indices);
    RMSD = RMSD(indices);
    Coordinates = Coordinates(:, indices);
    disp("A total of "+size(Angles,2)+" cluster centers.");

    candAngles = Angles(:,Scores<score_threshold);
    candScores = Scores(Scores<score_threshold);
    candRMSD = RMSD(Scores<score_threshold);
    candCoordinates = Coordinates(:,Scores<score_threshold);

    results{i} = struct('Angles', Angles, 'Coordinates', Coordinates, 'RMSD', RMSD, 'Scores', Scores, ...
                        'candAngles', candAngles, 'candScores', candScores, 'candRMSD', candRMSD, 'candCoordinates', candCoordinates);

    seconds = seconds + toc;
    disp("Clustering time for "+name+" is "+seconds+" seconds.")
end

for i = 1 : length(Names)
    name = Names{i};
    filename = sprintf('After_SA_clustering_%s_15res.mat', name);
    result = results{i};
    save(filename, '-struct', 'result');
end
end


% KABSCH_ALGORITHM calculates the optimal rigid body transformation
% that aligns two sets of 3D points (P and Q) using the Kabsch algorithm.
%
% Inputs:
%   P: Nx3 array of points to be aligned
%   Q: Nx3 array of target points
%
% Outputs:
%   R: 3x3 rotation matrix
%   t: 3x1 translation vector
%   rmsd: root-mean-square deviation between the aligned points
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


function indices = clustering(Coordinates)
n = 15;
N = size(Coordinates,2);
lib = 1 : N;
indices = [];

while ~isempty(lib)
    members = [];
    center = lib(1);
    centercoor = Coordinates(:,center);
    for i = 1 : length(lib)
        coor = Coordinates(:,lib(i));
        [~, ~, rmsd] = kabsch_algorithm(reshape(coor, [4*n,3]), reshape(centercoor, [4*n,3]));
        if rmsd < 0.5
            members = [members, lib(i)];
        end
    end
    lib = setdiff(lib, members);
    indices = [indices, center];
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

n = length(angles)/3;
coordinates = zeros(4*n, 3);
N = [0;0;0];
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n
    rotation = rotation*T(NR)*R(angles(3*i-2));
    Ca = N + rotation*[NCa; 0; 0];
    rotation = rotation*T(CaR)*R(angles(3*i-1));
    C = Ca + rotation*[CaC; 0; 0];
    rotation = rotation*T(CR)*R(angles(3*i));
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


function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end
