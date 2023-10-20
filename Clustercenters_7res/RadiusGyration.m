clc;
clear all;

n = 7;
% Names = ["design7.1", "design7.2", "design7.3", "c.11.18", "c.2.8", "c.3.45", ...
%     "c.4.35", "c.4.59", "c.4.78", "c.5.4", "c.8.1", "c.9.2"];
Names = ["144", "5191", "9104", "10632", "372", "9541", ...
    "1058", "3860", "5063", "9191", "7437", "9772", "1828", "7997"];

for ros = 1 : length(Names)
    % Get Rosetta designed structure coordinates
    name = Names(ros);
%     file = name+".pdb";
    file = "output_clustercenter"+name+"_0001_0001.pdb";
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
    name
    radius_of_gyration = sqrt(mean(squared_distances))
end
