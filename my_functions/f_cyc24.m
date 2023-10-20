% Function to compute the cyclic squared error
function error = f_cyc24(angles, R_omega, T_CR, T_NR, T_CaR, q)

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
phi21 = angles(41);
psi21 = angles(42);
phi22 = angles(43);
psi22 = angles(44);
phi23 = angles(45);
psi23 = angles(46);
phi24 = angles(47);
psi24 = angles(48);

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
s16 = Sigma(phi16, psi16, R_omega, T_CR, T_NR, T_CaR);
s17 = Sigma(phi17, psi17, R_omega, T_CR, T_NR, T_CaR);
s18 = Sigma(phi18, psi18, R_omega, T_CR, T_NR, T_CaR);
s19 = Sigma(phi19, psi19, R_omega, T_CR, T_NR, T_CaR);
s20 = Sigma(phi20, psi20, R_omega, T_CR, T_NR, T_CaR);
s21 = Sigma(phi21, psi21, R_omega, T_CR, T_NR, T_CaR);
s22 = Sigma(phi22, psi22, R_omega, T_CR, T_NR, T_CaR);
s23 = Sigma(phi23, psi23, R_omega, T_CR, T_NR, T_CaR);
s24 = Sigma(phi24, psi24, R_omega, T_CR, T_NR, T_CaR);

e1 = [1;0;0];
e2 = [0;1;0];
Rotation = s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*s19*s20*s21*s22*s23*s24;
Translation = s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*s19*s20*s21*s22*s23*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*s19*s20*s21*s22*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*s19*s20*s21*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*s19*s20*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*s19*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*s18*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*s17*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*s16*q + s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*s15*q + ...
    s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*s14*q + s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11*s12*s13*q + ...
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