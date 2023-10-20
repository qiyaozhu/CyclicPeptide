% Function to compute the cyclic squared error
function error = f_cyc15(angles, R_omega, T_CR, T_NR, T_CaR, q)

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

s1 = Sigma(phi1, psi1, R_omega, T_CR, T_NR, T_CaR);
s2 = Sigma(phi2, psi2, R_omega, T_CR, T_NR, T_CaR);
s3 = Sigma(phi3, psi3, R_omega, T_CR, T_NR, T_CaR);
s4 = Sigma(phi4, psi4, R_omega, T_CR, T_NR, T_CaR);
s5 = Sigma(phi5, psi5, R_omega, T_CR, T_NR, T_CaR);
s6 = Sigma(phi6, psi6, R_omega, T_CR, T_NR, T_CaR);
s7 = Sigma(phi7, psi7, R_omega, T_CR, T_NR, T_CaR);
s8 = Sigma(phi8, psi8, R_omega, T_CR, T_NR, T_CaR);
s9 = Sigma(phi9, psi9, R_omega, T_CR, T_NR, T_CaR);
s10 = Sigma(phi10, psi10, R_omega, T_CR, T_NR, T_CaR);
s11 = Sigma(phi11, psi11, R_omega, T_CR, T_NR, T_CaR);
s12 = Sigma(phi12, psi12, R_omega, T_CR, T_NR, T_CaR);
s13 = Sigma(phi13, psi13, R_omega, T_CR, T_NR, T_CaR);
s14 = Sigma(phi14, psi14, R_omega, T_CR, T_NR, T_CaR);
s15 = Sigma(phi15, psi15, R_omega, T_CR, T_NR, T_CaR);

e1 = [1;0;0];
e2 = [0;1;0];
Rotation = s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15;
Translation = s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*q + s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*q + s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*q + s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*q + s1*s2*s3*s4*s5*s6*s7*s8*q + s1*s2*s3*s4*s5*s6*s7*q + s1*s2*s3*s4*s5*s6*q + ...
    s1*s2*s3*s4*s5*q + s1*s2*s3*s4*q + s1*s2*s3*q + s1*s2*q + s1*q + q;
error = norm(Rotation*e1 - e1)^2 + norm(Rotation*e2 - e2)^2 + norm(Translation)^2;
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end

function sigma = Sigma(phi, psi, R_omega, T_CR, T_NR, T_CaR)
sigma = T_CR*R_omega*T_NR*R(phi)*T_CaR*R(psi);
end