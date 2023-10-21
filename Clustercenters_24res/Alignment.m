clc;
clear all;

n = 24;
Names = ["S4-1", "S4_inside_out_crystal"];

filename = "After_GA_clustercenterS4-Design_0001_0001_24res.mat";
data = load(filename);
Score = data.candScores;
Coordinates = data.candCoordinates;

for ros = 1 : length(Names)
    % Get Rosetta designed structure coordinates
    name = Names(ros);
    file = name+".pdb";
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

    RMSD = zeros(1,size(Coordinates,2));
    for i = 1 : size(Coordinates,2)
        coordinates = Coordinates(:,i);
        coordinates = reshape(coordinates, [4*n, 3]);
        [Rotation, translation, rmsd] = kabsch_algorithm(coordinates, native_backbone);
        RMSD(i) = rmsd;
    end
    [m,i] = min(RMSD)
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