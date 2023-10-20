function Clustering

addpath("/scratch/qz886/my_functions");

bestAngles = load("/scratch/qz886/bestAngles_7res.mat").bestAngles;
bestScores = load("/scratch/qz886/bestAngles_7res.mat").bestScores;

N = length(bestScores);
TorsionBins = zeros(N,7);

% get torsion strings
parfor i = 1 : N
    angles = bestAngles(:,i);
    for res = 1 : 7
        angle = angles(res*2-1:res*2);
        TorsionBins(i,res) = find_bin(angle);
    end
end

[TorsionStrings, PermAngles] = equivalent_bins(TorsionBins, bestAngles);
[TS, ~, index_TS] = unique(TorsionStrings, "rows"); % TorsionStrings = TS(index_TS,:)
ClusterIndices = zeros(1,size(TS,1)); % indices of centers in PermAngles

% choose the lowest-energy peptide as the cluster center
parfor ts = 1 : size(TS,1)
    AllIndex = 1 : N;
    lib = AllIndex(index_TS==ts);
    ClusterIndices(ts) = lib(1);
end

ClusterCenters = PermAngles(:,ClusterIndices);
ClusterScores = bestScores(:,ClusterIndices);

% sort the cluster centers based on energy
[ClusterScores, bestIndices] = mink(ClusterScores, length(ClusterScores));
ClusterCenters = ClusterCenters(:, bestIndices);

candPhi = zeros(7, length(ClusterScores));
candPsi = zeros(7, length(ClusterScores));
candCoordinates_x = zeros(7, length(ClusterScores));
candCoordinates_y = zeros(7, length(ClusterScores));
candCoordinates_z = zeros(7, length(ClusterScores));
candOmega = zeros(1, length(ClusterScores));

n = 7;
parfor cand = 1 : length(ClusterScores)
    backbone = reshape(ClusterCenters(:,cand), [2,n]);
    candPhi(:,cand) = backbone(1,:).';
    candPsi(:,cand) = backbone(2,:).';
    backbone = reshape(backbone, [2*n,1]);
    [coordinates, omega] = get_coordinates(backbone);
    candCoordinates_x(:,cand) = coordinates(:,1);
    candCoordinates_y(:,cand) = coordinates(:,2);
    candCoordinates_z(:,cand) = coordinates(:,3);
    candOmega(cand) = -omega;
end

candPhi = rad2deg(candPhi);
candPsi = rad2deg(candPsi);
candOmega = rad2deg(candOmega);

save('ClusterCenters.mat', 'ClusterCenters', 'ClusterScores', 'candPhi', 'candPsi', 'candCoordinates_x', 'candCoordinates_y', 'candCoordinates_z', 'candOmega');
end


% Find which torsion bin the angles [phi, psi] belong to
% Bin 1: [phi<0, psi<-1.39626], Bin 2: [phi<0, psi<1.39626], Bin 3: [phi<0],
% Bin 4: [phi>0, psi<-1.39626], Bin 5: [phi>0, psi<1.39626], Bin 6: [phi>0]
function bin = find_bin(angle)
phi = angle(1);
psi = angle(2);
if phi < 0
    if psi < -1.39626
        bin = 1;
    elseif psi < 1.39626
        bin = 2;
    else
        bin = 3;
    end
else
    if psi < -1.39626
        bin = 4;
    elseif psi < 1.39626
        bin = 5;
    else
        bin = 6;
    end
end
end


% Identify same torsion bin strings up to some circular permutation
function [TorsionStrings, PermAngles] = equivalent_bins(TorsionBins, bestAngles)
TorsionStrings = TorsionBins;
PermAngles = bestAngles;
for i = 1 : size(TorsionBins,1)
    bin = TorsionBins(i,:);
    angles = reshape(bestAngles(:,i), [2,7]);
    indices = find(bin==min(bin));
    if length(indices) ~= length(bin)
        count = 0;
        while length(indices) ~= 1 && count < length(bin)
            indices = mod(indices+1, length(bin));
            indices(indices==0) = length(bin);
            indices = indices(bin(indices)==min(bin(indices)));
            count = count + 1;
        end
        if length(indices) ~= 1
            indices = indices(1);
        end
        TorsionStrings(i,:) = circshift(bin, length(bin)-indices+count+1);
        PermAngles(:,i) = reshape(circshift(angles, length(bin)-indices+count+1, 2), [14, 1]);
    end
end
end


% Get CA coordiantes and the linking omega angle
function [coordinates, omega] = get_coordinates(angles)
% bond angles
NR = deg2rad(121.7);
CaR = deg2rad(111.2);
CR = deg2rad(116.2);

% bond lengths
NCa = 1.458;
CaC = 1.524;
CN = 1.329;

% omega torsion angle
omega = pi;

n = length(angles)/2;
coordinates = zeros(n, 3);
N = [0;0;0];
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n
    rotation = rotation*T(NR)*R(angles(2*i-1));
    Ca = N + rotation*[NCa; 0; 0];
    rotation = rotation*T(CaR)*R(angles(2*i));
    C = Ca + rotation*[CaC; 0; 0];
    rotation = rotation*T(CR)*R(omega);
    N = C + rotation*[CN; 0; 0];
    coordinates(i,:) = Ca.';
end

omega = torsion(coordinates(n,:).', C, [0;0;0], coordinates(1,:).');
end


function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
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
