% hbond_bb: Backbone-backbone hbonds
% Donor is NH (hbdon_PBA), acceptor is CO (hbacc_PBA).
% Need also sequence separation to read HBEval.csv (headers in schema),
% then know which fade intervals from HBFadeIntervals to use, and which
% polynomial coefficients from HBPoly1D to use.
% hbonds_geom.cc from Rosetta implements this calculation.
% @ E_HA = F_AHD * F_BAH * poly(d_HA)
% @ E_AHD = Fs_AH * F_BAH * poly_s(AHD)
% @ E_BAH = Fs_AH * F_AHD * poly_s(BAH)
% @ E_B2BAH = E_B2BAH_calc(chi, BAH)
% @ WH strength of donor, WA strength of acceptor
% @ energy = Sum[WH*WA*f(E_HA+E_AHD+E_BAH+E_B2BAH)], f smooth min
% @ coordinates: size 7nx3 matrix, 7n atoms (N, H, CA, 1HA, 2HA, C, O)
function [energy, count] = E_hbond(coordinates)

energy = 0;
count = 0; % number of hydrogen bonds
WH = 1.41;
WA = 1.08;
E_cutoff = -0.25;

% Fade function parameters
Fs_AH = [0, 0.1, 3.2, 3.3];
F_AHD = [-0.05, 0, 1, 1.05];
F_BAH = [-0.562949, 0, 1, 1.05];

n = size(coordinates,1) / 7; % number of residues
for i = 1 : n % acceptor residue number
    for j = 1 : n % donor residue number
        if j ~= i

            % relevant atom coordinates
            B2 = coordinates(i*7-4,:); % Ca
            B = coordinates(i*7-1,:); % C
            A = coordinates(i*7,:); % O
            H = coordinates(j*7-5,:); % H
            D = coordinates(j*7-6,:); % N

            % calculate distance and angles
            d_HA = norm(A-H);
            AHD = acos(dot((A-H)/norm(A-H), (D-H)/norm(D-H)));
            BAH = acos(dot((B-A)/norm(B-A), (H-A)/norm(H-A)));
            chi = torsion(H, A, B, B2);

            % fade functions
            fs_AH = Fade(d_HA, Fs_AH);
            f_AHD = Fade(-cos(AHD), F_AHD);
            f_BAH = Fade(-cos(BAH), F_BAH);

            % polynomial parameters for E_HA, E_AHD, and E_BAH
            xmin_HA = 1.384038127;
            xmax_HA = 2.998103943;
            fmin_HA = 1.1;
            fmax_HA = 1.1;
            C_HA = [-0.5307601, 6.47949946, -22.39522814, -55.14303544, 708.3094524, -2619.493182, 5227.88058, -6043.312116, 3806.046762, -1007.660241];

            xmin_AHD = 1.143564639;
            xmax_AHD = 3.1416;
            fmin_AHD = 1.1;
            fmax_AHD = 1.1;
            C_AHD = [0.47683259, -9.54524724, 83.62557693, -420.5586777, 1337.193549, -2786.262657, 3803.178227, -3278.628799, 1619.041162, -347.5015791];

            % calculate energies
            E_HA = f_AHD*f_BAH*polynomial1d(d_HA, xmin_HA, xmax_HA, fmin_HA, fmax_HA, C_HA);
            E_AHD = fs_AH*f_BAH*polynomial1d(AHD, xmin_AHD, xmax_AHD, fmin_AHD, fmax_AHD, C_AHD);
            E_B2BAH = E_B2BAH_calc(chi, BAH);
            if fs_AH ~= 0 && f_AHD ~= 0 && f_BAH ~= 0
                E = f(E_HA+E_AHD+E_B2BAH);
                energy = energy + WH*WA*E;
                if E < E_cutoff
                    count = count + 1;
                end
            end
        end
    end
end
end


function y = f(x)
if x < -0.1
    y = x;
elseif x < 0.1
    y = -0.025 + x/2 - 2.5*x^2;
else
    y = 0;
end
end


% Calculate energy E_B2BAH for SP2 hybrid
function energy = E_B2BAH_calc(chi, BAH)
d = 0.75;
m = 1.6;
l = 0.357;

H = (cos(2*chi)+1)/2;

if BAH > pi*2/3
    F = d/2*cos(3*(pi-BAH)) + (d-1)/2;
elseif BAH >= pi*(2/3-l)
    F = m/2*cos(pi-(2/3*pi-BAH)/l) + (m-1)/2;
else
    F = m - 1/2;
end

if BAH > pi*2/3
    G = d - 1/2;
elseif BAH >= pi*(2/3-l)
    G = (m-d)/2*cos(pi-(2/3*pi-BAH)/l) + (m+d-1)/2;
else
    G = m - 1/2;
end

energy = H*F + (1-H)*G;
end


% Fade function given I = [a,b,c,d]
% f=0 if x<=a or x>=d
% f=1 if b<=x<=c
% sigmoid function on [a,b] and [c,d] to ensure continuous derivative
function y = Fade(x, I)
a = I(1);
b = I(2);
c = I(3);
d = I(4);

if x <= a
    y = 0;
elseif x <= b
    z = (x-a) / (b-a);
    y = -2*z^3 + 3*z^2;
elseif x <= c
    y = 1;
elseif x <= d
    z = (x-c) / (d-c);
    y = 2*z^3 - 3*z^2 + 1;
else
    y = 0;
end
end


% Polynomial evaluation
% f=fmin if x<=xmin
% f=fmax if x>=xmax
% f=C[1]*x^(n-1) + C[2]*x^(n-2) + ... + C[n-1]*x + C[n] if xmin<x<xmax
function y = polynomial1d(x, xmin, xmax, fmin, fmax, C)
if x <= xmin
    y = fmin;
elseif x <= xmax
    y = C(1);
    for i = 2 : length(C)
        y = y*x + C(i);
    end
else
    y = fmax;
end
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
