% Function to find a feasible direction
function [p_f, angles_f] = feasible20(p, angles, Boundaries)

angles_periodic = periodic(angles + p);
in_or_not = inBoundary(angles_periodic, Boundaries);
in_or_not = repelem(in_or_not, 2);
p_f = p.*in_or_not.';
angles_f = periodic(angles + p_f);

end


% Function to make the angles within [-pi, pi]
function angles_periodic = periodic(angles)

angles_periodic = angles;
for i = 1 : length(angles)
    if angles(i) < -pi
        angles_periodic(i) = angles(i) + 2*pi;
    elseif angles(i) > pi
        angles_periodic(i) = angles(i) - 2*pi;
    end
end
end

% Function to check if the torsion angles are in the Ramachandran plot
function in_or_not = inBoundary(angles, Boundaries)

phi1 = angles(1);
psi1 = angles(2);
phi2 = angles(3);
psi2 = angles(4);
phi3 = angles(5);
psi3 = angles(6);
phi4 = angles(7);
psi4 = angles(8);
phi5 = angles(9);
psi5 = angles(10);
phi6 = angles(11);
psi6 = angles(12);
phi7 = angles(13);
psi7 = angles(14);
phi8 = angles(15);
psi8 = angles(16);
phi9 = angles(17);
psi9 = angles(18);
phi10 = angles(19);
psi10 = angles(20);
phi11 = angles(21);
psi11 = angles(22);
phi12 = angles(23);
psi12 = angles(24);
phi13 = angles(25);
psi13 = angles(26);
phi14 = angles(27);
psi14 = angles(28);
phi15 = angles(29);
psi15 = angles(30);
phi16 = angles(31);
psi16 = angles(32);
phi17 = angles(33);
psi17 = angles(34);
phi18 = angles(35);
psi18 = angles(36);
phi19 = angles(37);
psi19 = angles(38);
phi20 = angles(39);
psi20 = angles(40);

in_or_not = [inpolygon(phi1, psi1, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi2, psi2, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi3, psi3, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi4, psi4, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi5, psi5, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi6, psi6, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi7, psi7, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi8, psi8, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi9, psi9, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi10, psi10, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi11, psi11, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi12, psi12, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi13, psi13, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi14, psi14, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi15, psi15, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi16, psi16, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi17, psi17, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi18, psi18, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi19, psi19, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi20, psi20, Boundaries(1,:), Boundaries(2,:))];
end