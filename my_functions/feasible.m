% Function to find a feasible direction
function [p_f, angles_f] = feasible(p, angles, Boundaries)

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

in_or_not = [inpolygon(phi1, psi1, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi2, psi2, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi3, psi3, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi4, psi4, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi5, psi5, Boundaries(1,:), Boundaries(2,:)), inpolygon(phi6, psi6, Boundaries(1,:), Boundaries(2,:)), ...
    inpolygon(phi7, psi7, Boundaries(1,:), Boundaries(2,:))];
end