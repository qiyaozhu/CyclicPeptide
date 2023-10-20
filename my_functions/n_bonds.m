% Parameters for energy functions
% Calculate the number of bonds between any atom pair
% @ n is the number of residues, 7n atoms (N, H, CA, 1HA, 2HA, C, O)
% @ D is a 7nx7n matrix
function D = n_bonds(n)

A = zeros(7*n, 7*n); % adjacency matrix

for i = 1 : n
    A(7*i-6, 7*i-5) = 1; % N
    A(7*i-6, 7*i-4) = 1;
    if i == 1
        A(1, 7*n-1) = 1;
    else
        A(7*i-6, 7*i-8) = 1;
    end
    A(7*i-5, 7*i-6) = 1; % H
    A(7*i-4, 7*i-6) = 1; % Ca
    A(7*i-4, 7*i-3) = 1;
    A(7*i-4, 7*i-2) = 1;
    A(7*i-4, 7*i-1) = 1;
    A(7*i-3, 7*i-4) = 1; % 1HA
    A(7*i-2, 7*i-4) = 1; % 2HA
    A(7*i-1, 7*i-4) = 1; % C
    A(7*i-1, 7*i) = 1;
    if i == n
        A(7*n-1, 1) = 1;
    else
        A(7*i-1, 7*i+1) = 1;
    end
    A(7*i, 7*i-1) = 1; % O
end

G = graph(A);
D = distances(G);
end