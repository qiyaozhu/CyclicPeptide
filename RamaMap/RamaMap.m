%% For L amino acids
clc;
clear all;

res = 'val';

fname = sprintf('%s.dat', res);
file = fopen(fname,'rt');
C = textscan(file, '%s%f%f%f%f');
fclose(file);

Phi = C{2};
Psi = C{3};
Energy = C{5};

F = reshape(Energy, 36, 36);
F = [F, F(:,1)];
F = [F; F(1,:)];
Energy = reshape(F, 37*37, 1);

P1 = reshape(Phi, 36, 36);
P1 = [P1, P1(:,1)];
P1 = [P1; -P1(1,:)];
Phi = reshape(P1, 37*37, 1);

P2 = reshape(Psi, 36, 36);
P2 = [P2, -P2(:,1)];
P2 = [P2; P2(1,:)];
Psi = reshape(P2, 37*37, 1);

% For each square in the grid, use the lower left corner as location
% indicator, and record if this square is within the good region
Matrix = (F<=10);
Map = zeros(37, 37);
for i = 1 : 36
    for j = 1 : 36
        if Matrix(i,j) == 1
            if Matrix(i+1,j) == 1 && Matrix(i,j+1) == 1 && Matrix(i+1,j+1) == 1
                Map(i,j) = 1;
            end
        end
    end
end

for i = 1 : 36
    if Matrix(i,37) == 1
        if Matrix(i+1,37) == 1
            Map(i,37) = 1;
        end
    end
end

for j = 1 : 36
    if Matrix(37,j) == 1
        if Matrix(37,j+1) == 1
            Map(37,j) = 1;
        end
    end
end

Map(37,37) = Matrix(37,37);
fname = sprintf('ramamap_%s.mat', res);
save(fname, 'Map', 'F');

Boundaries = rad2deg(load('C:\Users\Cathe\OneDrive\桌面\Cyclic Peptide\ramabin_glycine.mat').Boundaries.');

figure;
set(gcf,'color','w');
colormap autumn;
pointsize = 10;
hold on;
scatter(Phi(Energy<=10), Psi(Energy<=10), pointsize, Energy(Energy<=10), 'filled');
plot(Boundaries(1,:), Boundaries(2,:), 'k-');
hold off;
colorbar;
xlim([-180, 180]);
ylim([-180, 180]);
xticks([-180 -120 -60 0 60 120 180]);
yticks([-180 -120 -60 0 60 120 180]);
xlabel('$\phi$','interpreter','latex');
ylabel('$\psi$','interpreter','latex');
set(gca,'FontSize',16,'FontWeight','Bold');


%% For D amino acids
clc;
clear all;

res = 'val';

fname = sprintf('%s.dat', res);
file = fopen(fname,'rt');
C = textscan(file, '%s%f%f%f%f');
fclose(file);

Phi = C{2};
Psi = C{3};
Energy = C{5};

F = reshape(Energy, 36, 36);
F = [F, F(:,1)];
F = [F; F(1,:)];
E = zeros(37,37);
for i = 1 : 37
    for j = 1 : 37
        E(i,j) = F(38-i,38-j);
    end
end
Energy = reshape(E, 37*37, 1);

P1 = reshape(Phi, 36, 36);
P1 = [P1, P1(:,1)];
P1 = [P1; -P1(1,:)];
Phi = reshape(P1, 37*37, 1);

P2 = reshape(Psi, 36, 36);
P2 = [P2, -P2(:,1)];
P2 = [P2; P2(1,:)];
Psi = reshape(P2, 37*37, 1);

% For each square in the grid, use the lower left corner as location
% indicator, and record if this square is within the good region
Matrix = (E<=10);
Map = zeros(37, 37);
for i = 1 : 36
    for j = 1 : 36
        if Matrix(i,j) == 1
            if Matrix(i+1,j) == 1 && Matrix(i,j+1) == 1 && Matrix(i+1,j+1) == 1
                Map(i,j) = 1;
            end
        end
    end
end

for i = 1 : 36
    if Matrix(i,37) == 1
        if Matrix(i+1,37) == 1
            Map(i,37) = 1;
        end
    end
end

for j = 1 : 36
    if Matrix(37,j) == 1
        if Matrix(37,j+1) == 1
            Map(37,j) = 1;
        end
    end
end

Map(37,37) = Matrix(37,37);
fname = sprintf('ramamap_d%s.mat', res);
save(fname, 'Map', 'E');

Boundaries = rad2deg(load('C:\Users\Cathe\OneDrive\桌面\Cyclic Peptide\ramabin_glycine.mat').Boundaries.');

figure;
set(gcf,'color','w');
colormap autumn;
pointsize = 10;
hold on;
scatter(Phi(Energy<=10), Psi(Energy<=10), pointsize, Energy(Energy<=10), 'filled');
plot(Boundaries(1,:), Boundaries(2,:), 'k-');
hold off;
colorbar;
xlim([-180, 180]);
ylim([-180, 180]);
xticks([-180 -120 -60 0 60 120 180]);
yticks([-180 -120 -60 0 60 120 180]);
xlabel('$\phi$','interpreter','latex');
ylabel('$\psi$','interpreter','latex');
set(gca,'FontSize',16,'FontWeight','Bold');


%%
angles = [-2.0444
   -1.8766
   -1.9613
    2.1789
    1.0456
   -2.4087
   -0.8810
   -0.2430
   -0.8526
   -0.4267
    1.3858
   -1.3689
   -1.5966
   -0.0614];
angles = rad2deg(reshape(angles, [2,7]));

seq = ["lys", "tyr", "dpro", "glu", "glu", "dpro", "dlys"];

for i = 1 : 7
    res = char(seq(i));
    isD = false;
    if res(1) == 'd'
        res = res(2:end);
        isD = true;
    end
    fname = sprintf('%s.dat', res);
    file = fopen(fname,'rt');
    C = textscan(file, '%s%f%f%f%f');
    fclose(file);

    if isD
        Phi = -C{2};
        Psi = -C{3};
    else
        Phi = C{2};
        Psi = C{3};
    end

    Prob = C{4};
    Energy = C{5};

    F = reshape(Energy, 36, 36);
    F = [F, F(:,1)];
    F = [F; F(1,:)];
    Energy = reshape(F, 37*37, 1);

    P1 = reshape(Phi, 36, 36);
    P1 = [P1, P1(:,1)];
    P1 = [P1; -P1(1,:)];
    Phi = reshape(P1, 37*37, 1);

    P2 = reshape(Psi, 36, 36);
    P2 = [P2, -P2(:,1)];
    P2 = [P2; P2(1,:)];
    Psi = reshape(P2, 37*37, 1);

    figure;
    set(gcf,'color','w');
    colormap autumn;
    pointsize = 10;
    hold on;
    scatter(Phi(Energy<=10), Psi(Energy<=10), pointsize, Energy(Energy<=10), 'filled');
    plot(angles(1,i), angles(2,i), 'kx', 'LineWidth', 2, 'MarkerSize', 10);
    hold off;
    colorbar;
    xlim([-180, 180]);
    ylim([-180, 180]);
    xticks([-180 -120 -60 0 60 120 180]);
    yticks([-180 -120 -60 0 60 120 180]);
    xlabel('$\phi$','interpreter','latex');
    ylabel('$\psi$','interpreter','latex');
    set(gca,'FontSize',16,'FontWeight','Bold');
end