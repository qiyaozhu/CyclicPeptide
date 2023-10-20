% Function to compute the cyclic squared error
function error = f_cyc(angles, R_omega, T_CR, T_NR, T_CaR, q)

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

s1 = T_CR*R_omega*T_NR*R(phi1)*T_CaR*R(psi1);
s2 = T_CR*R_omega*T_NR*R(phi2)*T_CaR*R(psi2);
s3 = T_CR*R_omega*T_NR*R(phi3)*T_CaR*R(psi3);
s4 = T_CR*R_omega*T_NR*R(phi4)*T_CaR*R(psi4);
s5 = T_CR*R_omega*T_NR*R(phi5)*T_CaR*R(psi5);
s6 = T_CR*R_omega*T_NR*R(phi6)*T_CaR*R(psi6);
s7 = T_CR*R_omega*T_NR*R(phi7)*T_CaR*R(psi7);

e1 = [1;0;0];
e2 = [0;1;0];
Rotation = s1*s2*s3*s4*s5*s6*s7;
Translation = s1*s2*s3*s4*s5*s6*q + s1*s2*s3*s4*s5*q + s1*s2*s3*s4*q + s1*s2*s3*q + s1*s2*q + s1*q + q;
error = norm(Rotation*e1 - e1)^2 + norm(Rotation*e2 - e2)^2 + norm(Translation)^2;
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end