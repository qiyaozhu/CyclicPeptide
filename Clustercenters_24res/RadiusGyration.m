clc;
clear all;

n = 24;

% FILE = fopen("Pnear_list.txt", 'r');
% data = textscan(FILE, '%s\t%f\t%s\n');
% fclose(FILE);
% CenterNames = data{1};
CenterNames = ["LowEnergy21698_0001", "LowEnergy10052_0001", "LowEnergy759_0001"];
Rg = zeros(1,length(CenterNames));

for ros = 1 : length(CenterNames)
    % Get Rosetta designed structure coordinates
    name = CenterNames(ros)
    % file = "output_"+name+".pdb";
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

    % Calculate the center of mass
    center_mass = mean(native_backbone);

    % Calculate the distance of each atom from the center of mass
    displacement_vectors = native_backbone - center_mass;
    squared_distances = sum(displacement_vectors.^2, 2);

    % Calculate the radius of gyration
    radius_of_gyration = sqrt(mean(squared_distances))
    Rg(ros) = radius_of_gyration;
end

% save("Rg.mat", "Rg");