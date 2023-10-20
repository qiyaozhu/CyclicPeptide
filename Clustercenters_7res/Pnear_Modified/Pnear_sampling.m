function Pnear_sampling(fname)

file = fopen(fname,'rt');
C = textscan(file, '%s%s');
fclose(file);

Seqs = C{1};
Names = C{2};

bestAngles = load("/scratch/qz886/bestAngles_7res.mat").bestAngles;
bestScores = load("/scratch/qz886/bestAngles_7res.mat").bestScores;
bestLocations = floor(rad2deg(bestAngles)/10) + 19;

parfor dc = 1 : length(Names)
    sprintf('%d', dc)
    seq = Seqs(dc);
    seq = strsplit(seq{1}, '-');
    n = length(seq);
    RamaMap = zeros(37,37,n);

    for res = 1 : n
        RamaMap(:,:,res) = load("/scratch/qz886/RamaMap/ramamap_"+seq(res)+".mat").Map;
    end

    candAngles = [];
    candScores = [];
    LocationIndices = randperm(length(bestScores), length(bestScores));
    loc = 1;

    while loc <= length(LocationIndices) && length(candScores) <= 5000
        locations = reshape(bestLocations(:,LocationIndices(loc)), [2,n]);
        startres = [1:n];
        seqpos = 1;

        while ~isempty(startres) && seqpos<=n
            remove = [];
            for i = 1 : length(startres)
                respos = startres(i)+seqpos-1;
                if respos > 7
                    respos = respos - 7;
                end
                phipos = locations(1,respos);
                psipos = locations(2,respos);
                if RamaMap(phipos, psipos, seqpos) == 0
                    remove = [remove, startres(i)];
                end
            end
            startres = setdiff(startres, remove);
            seqpos = seqpos + 1;
        end

        angles = reshape(bestAngles(:,LocationIndices(loc)), [2,n]);
        score = bestScores(LocationIndices(loc));
        for k = 1 : length(startres)
            backbone = [angles(:,startres(k):end), angles(:,1:startres(k)-1)];
            backbone = reshape(backbone, [2*n,1]);
            candAngles = [candAngles, backbone];
            candScores = [candScores, score];
        end
        loc = loc + 1;
    end

    candPhi = zeros(7, length(candScores));
    candPsi = zeros(7, length(candScores));
    candCoordinates_x = zeros(7, length(candScores));
    candCoordinates_y = zeros(7, length(candScores));
    candCoordinates_z = zeros(7, length(candScores));
    candOmega = zeros(1, length(candScores));

    for cand = 1 : length(candScores)
        backbone = reshape(candAngles(:,cand), [2,n]);
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

    results{dc}{1} = candPhi;
    results{dc}{2} = candPsi;
    results{dc}{3} = candCoordinates_x;
    results{dc}{4} = candCoordinates_y;
    results{dc}{5} = candCoordinates_z;
    results{dc}{6} = candOmega;
    results{dc}{7} = candScores;
end

for dc = 1 : length(Names)
    name = Names(dc);
    filename = sprintf('Pnear_test_%s.mat', name{1});
    candPhi = results{dc}{1};
    candPsi = results{dc}{2};
    candCoordinates_x = results{dc}{3};
    candCoordinates_y = results{dc}{4};
    candCoordinates_z = results{dc}{5};
    candOmega = results{dc}{6};
    candScores = results{dc}{7};
    save(filename, 'candPhi', 'candPsi', 'candCoordinates_x', 'candCoordinates_y', 'candCoordinates_z', 'candOmega', 'candScores');
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
