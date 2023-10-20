% Parameter matrix for Lennard-Jones and Coulombic electrostatic potential
% @ n is the number of residues, 7n atoms (N, H, CA, 1HA, 2HA, C, O)
% @ LJ_radius, LJ_well, Coulumb are 1x7n parameter vectors
function [LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V] = FA_parameter(n)

LJ_radius = zeros(1,7*n);
LJ_well = zeros(1,7*n);
Coulumb = zeros(1,7*n);
LK_DG = zeros(1,7*n);
LK_Lambda = zeros(1,7*n);
LK_V = zeros(1,7*n);

% LJ radius
s_N = 1.802452;
s_H = 0.901681;
s_Ca = 2.011760;
s_HA = 1.421272;
s_C = 1.916661;
s_O = 1.540580;

% LJ well depth
e_N = 0.161725;
e_H = 0.005000;
e_Ca = 0.062642;
e_HA = 0.021808;
e_C = 0.141799;
e_O = 0.142417;

% partial atomic charges
q_N = -0.6046255;
q_H = 0.3987955;
q_Ca = -0.0257287; % this is different for other residues
q_HA = 0.1157793;
q_C = 0.6884871;
q_O = -0.6884871;

% vapor-to-water transfer free energy Delta_G_free
DG_N = -9.969494;
DG_H = 0;
DG_Ca = 2.533791;
DG_HA = 0;
DG_C = 3.104248;
DG_O = -8.006829;

% correlation length lambda
lam_N = 3.5;
lam_H = 3.5;
lam_Ca = 3.5;
lam_HA = 3.5;
lam_C = 3.5;
lam_O = 3.5;

% atomic volume
V_N = 15.992;
V_H = 0;
V_Ca = 12.137;
V_HA = 0;
V_C = 13.221;
V_O = 12.196;

for i = 0 : n-1
    LJ_radius(7*i+1) = s_N;
    LJ_radius(7*i+2) = s_H;
    LJ_radius(7*i+3) = s_Ca;
    LJ_radius(7*i+4) = s_HA;
    LJ_radius(7*i+5) = s_HA;
    LJ_radius(7*i+6) = s_C;
    LJ_radius(7*i+7) = s_O;

    LJ_well(7*i+1) = e_N;
    LJ_well(7*i+2) = e_H;
    LJ_well(7*i+3) = e_Ca;
    LJ_well(7*i+4) = e_HA;
    LJ_well(7*i+5) = e_HA;
    LJ_well(7*i+6) = e_C;
    LJ_well(7*i+7) = e_O;

    Coulumb(7*i+1) = q_N;
    Coulumb(7*i+2) = q_H;
    Coulumb(7*i+3) = q_Ca;
    Coulumb(7*i+4) = q_HA;
    Coulumb(7*i+5) = q_HA;
    Coulumb(7*i+6) = q_C;
    Coulumb(7*i+7) = q_O;

    LK_DG(7*i+1) = DG_N;
    LK_DG(7*i+2) = DG_H;
    LK_DG(7*i+3) = DG_Ca;
    LK_DG(7*i+4) = DG_HA;
    LK_DG(7*i+5) = DG_HA;
    LK_DG(7*i+6) = DG_C;
    LK_DG(7*i+7) = DG_O;

    LK_Lambda(7*i+1) = lam_N;
    LK_Lambda(7*i+2) = lam_H;
    LK_Lambda(7*i+3) = lam_Ca;
    LK_Lambda(7*i+4) = lam_HA;
    LK_Lambda(7*i+5) = lam_HA;
    LK_Lambda(7*i+6) = lam_C;
    LK_Lambda(7*i+7) = lam_O;

    LK_V(7*i+1) = V_N;
    LK_V(7*i+2) = V_H;
    LK_V(7*i+3) = V_Ca;
    LK_V(7*i+4) = V_HA;
    LK_V(7*i+5) = V_HA;
    LK_V(7*i+6) = V_C;
    LK_V(7*i+7) = V_O;
end
end