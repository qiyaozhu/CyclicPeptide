% rama (0.25): Probability of backbone ϕ, ψ angles given the amino acid type & p_aa_pp (0.4): Probability of amino acid at Φ/Ψ
% @ angles = [phi_1; psi_1; phi_2; psi_2; ... phi_n; psi_n]
% @ A: coefficients for bicubic polynomials p_ij(x,y)=ΣΣa_mn*((x-xi)/delta)^m*((y-yj)/delta)^n
% @ X = [-180, -170, ..., 170, 180], Y = [-180, -170, ..., 170, 180], delta = 10 degree
% @ w = weight
function energy = E_rama(angles, A)

energy = 0;

X = linspace(-pi, pi, 37);
Y = linspace(-pi, pi, 37);
delta = deg2rad(10);

n = length(angles)/2;

for k = 1 : n
    phi = angles(k*2-1);
    psi = angles(k*2);
    i = floor((phi+pi)/delta)+1;
    j = floor((psi+pi)/delta)+1;
    if i == 37
        i = 36;
    end
    if j == 37
        j = 36;
    end
    xi = X(i);
    yj = Y(j);
    a = A(i,j,:);
    bicubic_interp = a(1) + a(2)*((phi-xi)/delta) + a(3)*((phi-xi)/delta).^2 + a(4)*((phi-xi)/delta).^3 + ...
    a(5)*((psi-yj)/delta) + a(6)*((phi-xi)/delta).*((psi-yj)/delta) + a(7)*((phi-xi)/delta).^2.*((psi-yj)/delta) + a(8)*((phi-xi)/delta).^3.*((psi-yj)/delta) + ...
    a(9)*((psi-yj)/delta).^2 + a(10)*((phi-xi)/delta).*((psi-yj)/delta).^2 + a(11)*((phi-xi)/delta).^2.*((psi-yj)/delta).^2 + a(12)*((phi-xi)/delta).^3.*((psi-yj)/delta).^2 + ...
    a(13)*((psi-yj)/delta).^3 + a(14)*((phi-xi)/delta).*((psi-yj)/delta).^3 + a(15)*((phi-xi)/delta).^2.*((psi-yj)/delta).^3 + a(16)*((phi-xi)/delta).^3.*((psi-yj)/delta).^3;
    energy = energy + bicubic_interp;
end
end