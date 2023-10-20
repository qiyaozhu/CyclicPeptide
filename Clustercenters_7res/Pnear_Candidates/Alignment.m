function Alignment

addpath("/scratch/qz886/my_functions");
n = 7;
% Names = ["design7.1", "design7.2", "design7.3", "c.11.18", "c.2.8", "c.3.45", ...
%    "c.4.35", "c.4.59", "c.4.78", "c.5.4", "c.8.1", "c.9.2"];
Names = ["6be9", "6bew", "6bf3", "6bf5"];
FILE = fopen("Pnear_list.txt", 'r');
data = textscan(FILE, '%s\t%f\t%s\n');
fclose(FILE);
DesignNames = data{1};

DDNames = {};
DNRMSD = zeros(1, length(Names));
parfor ros = 1 : length(Names)
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

    % Perform cyclic permutation to align best with the designed structure
    dnRMSD = 100;
    candRMSD = 0;
    align_peptide = backbone;
    for cand = 1 : length(DesignNames)
        file = "output_"+DesignNames(cand)+".pdb";
        backbone = pdb2mat(file);
        coordinate1 = [backbone.X.', backbone.Y.', backbone.Z.'];
        atomName = backbone.atomName;
        design_backbone = zeros(4*n,3);

        N_coor = coordinate1(atomName=="N",:);
        Ca_coor = coordinate1(atomName=="CA",:);
        C_coor = coordinate1(atomName=="C",:);
        O_coor = coordinate1(atomName=="O",:);

        for res = 1 : n
            design_backbone(4*res-3,:) = N_coor(res,:);
            design_backbone(4*res-2,:) = Ca_coor(res,:);
            design_backbone(4*res-1,:) = C_coor(res,:);
            design_backbone(4*res,:) = O_coor(res,:);
        end

        [R, t, rmsd] = kabsch_algorithm(design_backbone, native_backbone);
        for res = 1 : n-1
            design_backbone1 = [design_backbone(4*res+1:end,:); design_backbone(1:4*res,:)];
            [R1, t1, rmsd1] = kabsch_algorithm(design_backbone1, native_backbone);
            if rmsd1 < rmsd
                rmsd = rmsd1;
                R = R1;
                t = t1;
            end
        end

        if rmsd < dnRMSD
            dnRMSD = rmsd;
            candRMSD = cand;
            coordinate_align = (R * coordinate1')' + t';
            align_peptide = backbone;
            align_peptide.X = coordinate_align(:,1);
            align_peptide.Y = coordinate_align(:,2);
            align_peptide.Z = coordinate_align(:,3);
            align_peptide.outfile = name+"_align_"+DesignNames(candRMSD)+".pdb";
        end
    end

    FILE = fopen(align_peptide.outfile, 'w');
    fprintf(FILE, 'SEQRES   1 A  7\n');
    fclose(FILE);
    mat2pdb(align_peptide);
    DDNames{ros} = DesignNames(candRMSD);
    DNRMSD(ros) = dnRMSD;
end

fileID = fopen("Design_align.txt",'w');
for ros = 1 : length(Names)
    name = Names(ros);
    ddname = DDNames{ros};
    dnRMSD = DNRMSD(ros);
    fprintf(fileID,'%s\t%s\t%.3f\n', name, ddname{1}, dnRMSD);
end
fclose(fileID);
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
