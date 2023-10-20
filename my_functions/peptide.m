%% Glycine each residue has 7 atoms: N, H, CA, 1HA, 2HA, C, O
function coordinates = peptide(angles)
% bond angles
NR = deg2rad(121.7);
CaR = deg2rad(111.2);
CaHR = deg2rad(109.5);
CR = deg2rad(116.2);
COR = deg2rad(123.0);
NHR = deg2rad(119.2);

% bond lengths
NCa = 1.458;
CaH = 1.090;
CaC = 1.524;
CN = 1.329;
CO = 1.231;
NH = 1.010;

% omega torsion angle
omega = pi;

n = length(angles)/2;
coordinates = zeros(7*n, 3);
N = [0;0;0];
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n
    H = N + rotation*[NH*cos(pi-NHR); -NH*sin(pi-NHR); 0];
    rotation = rotation*T(NR)*R(angles(2*i-1));
    Ca = N + rotation*[NCa; 0; 0];
    rotation = rotation*T(CaR)*R(angles(2*i));
    C = Ca + rotation*[CaC; 0; 0];
    rotation = rotation*T(CR)*R(omega);
    O = C + rotation*[-CO*cos(pi-COR); CO*sin(pi-COR); 0];
    N = C + rotation*[CN; 0; 0];

    coordinates(7*i-5,:) = H.';
    coordinates(7*i-4,:) = Ca.';
    coordinates(7*i-1,:) = C.';
    coordinates(7*i,:) = O.';
    
    if i < n
        coordinates(7*i+1,:) = N.';
    end
end

% add the two hydrogen atoms at each Ca
for i = 1 : n
    Ncoord = coordinates(7*i-6,:);
    CAcoord = coordinates(7*i-4,:);
    Ccoord = coordinates(7*i-1,:);
    hor = -(Ncoord + (Ccoord-Ncoord)*1.20315/2.46077 - CAcoord);
    horcomp = hor / norm(hor) * CaH*cos(0.938693);
    pep = cross(Ccoord-CAcoord, Ncoord-CAcoord);
    pepcomp = pep / norm(pep) * CaH*sin(0.938693);
    coordinates(7*i-3,:) = CAcoord + horcomp + pepcomp;
    coordinates(7*i-2,:) = CAcoord + horcomp - pepcomp;
end
end


function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end