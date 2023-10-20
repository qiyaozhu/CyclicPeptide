function Pnear_sampling_GA_20res(seq, name)

addpath("/scratch/qz886/my_functions");
n = 20;
seq = strsplit(seq, '-');

% Get corresponding Rama information for the given sequence
RamaMap = zeros(37,37,n);
for res = 1 : n
    RamaMap(:,:,res) = load("/scratch/qz886/RamaMap/ramamap_"+seq(res)+".mat").Map;
end

M = 50; % maximum number of iterations
stall_threshold = 3;
score_threshold = 0;

% filename = sprintf('After_SA_%s_Angles_20res.txt', name);
% Peptides = readmatrix(filename).';
% Scores = Peptides(3*n+1,:);
% RMSD = Peptides(end,:);
% Angles = Peptides(1:3*n,:);
% Angles = deg2rad(Angles);
% 
% N = size(Angles, 2);
% Coordinates = zeros(4*n, 3*N);
% for cand = 1 : N
%     angle = Angles(:,cand);
%     Coordinates(:,3*(cand-1)+1:3*cand) = get_coordinates(angle);
% end
% Coordinates = reshape(Coordinates, [4*n*3, N]);
% disp("A total population of "+N+" after SA.");
% 
% % Clustering to select lowest energy centers
% [Scores, bestIndices] = mink(Scores, 30000);
% Angles = Angles(:, bestIndices);
% RMSD = RMSD(bestIndices);
% Coordinates = Coordinates(:, bestIndices);
% 
% indices = clustering(Coordinates);
% Angles = Angles(:, indices);
% Scores = Scores(indices);
% RMSD = RMSD(indices);
% Coordinates = Coordinates(:, indices);
% disp("A total of "+size(Angles,2)+" cluster centers.");
% 
% candAngles = Angles(:,Scores<score_threshold);
% candScores = Scores(Scores<score_threshold);
% candRMSD = RMSD(Scores<score_threshold);
% candCoordinates = Coordinates(:,Scores<score_threshold);
% 
% filename = sprintf('After_SA_clustering_%s_20res.mat', name);
% save(filename, 'Angles', 'Coordinates', 'RMSD', 'Scores', 'candAngles', 'candScores', 'candRMSD', 'candCoordinates');

filename = sprintf('After_SA_clustering_%s_20res.mat', name);
data = load(filename);
Angles = data.Angles;
Coordinates = data.Coordinates;
RMSD = data.RMSD;
Scores = data.Scores;
candAngles = data.candAngles;
candScores = data.candScores;
candRMSD = data.candRMSD;
candCoordinates = data.candCoordinates;

population = 750;
[Scores, bestIndices] = mink(Scores, population);
Angles = Angles(:, bestIndices);
RMSD = RMSD(bestIndices);
Coordinates = Coordinates(:, bestIndices);
disp("Min RMSD is "+min(RMSD)+", max RMSD is "+max(RMSD)+", and mean RMSD is "+mean(RMSD)+".");
disp("Min Energy is "+min(Scores)+", max Energy is "+max(Scores)+", and mean Energy is "+mean(Scores)+".");

N = size(Coordinates, 2);
Coordinates = reshape(Coordinates, [4*n, 3*N]);

cyc = zeros(1,size(Coordinates,2)/3);
for ci = 1 : size(Coordinates,2)/3
    coordinates = Coordinates(:,3*ci-2:3*ci);
    cyc(ci) = norm(coordinates(end-1,:)-coordinates(1,:));
end
disp("Max cyclic error is "+max(cyc)+".");

minutes = 0;
seconds = 0;

i = 0;
stall = 0;
while i <= M && stall < stall_threshold
    fprintf('Generation %d\n', i);
    xover_region = [3,10];
    mut_region = [3,10];
    xover_size = 1.5*population;
    mut_size = 1.5*population;

    tic;
    [Coordinates_new, Labels] = GA(Angles, Coordinates, RamaMap, xover_region, mut_region, xover_size, mut_size);
    Coordinates2PDB(Coordinates_new, name, seq, i);
    N = size(Coordinates_new,2)/3;
    disp("Generation"+i+" has a population of "+N+". Of which "+length(Labels(Labels==1))+" are xovers, and "+length(Labels(Labels==2))+" are mutants.");
    seconds = seconds + toc;

    cyc = zeros(1,size(Coordinates_new,2)/3);
    for ci = 1 : size(Coordinates_new,2)/3
        coordinates = Coordinates_new(:,3*ci-2:3*ci);
        cyc(ci) = norm(coordinates(end-1,:)-coordinates(1,:));
    end
    disp("Max cyclic error is "+max(cyc)+".");

    FILE = fopen("relax_"+name+"_Generation"+i+"_20res.flags", 'w');
    fprintf(FILE, "-in:file:l "+name+"_Generation"+i+"_list_20res.txt\n");
    fprintf(FILE, '-in:file:fullatom\n');
    fprintf(FILE, "-in:file:native /scratch/qz886/Clustercenters_20res/Pnear_Candidates/output_"+name+".pdb\n");
    fprintf(FILE, '-parser:protocol relax_Cartesian.xml\n');
    fprintf(FILE, '-parser:script_vars Nres=20\n');
    fprintf(FILE, "-out:file:scorefile score_"+name+"_Generation"+i+"_20res.sc\n");
    fclose(FILE);

    FILE = fopen(name+"_Generation"+i+"_list_20res.txt", 'w');
    for cand = 1 : N
        fprintf(FILE, name+"_Generation"+i+"_cand"+cand+"_20res.pdb\n");
    end
    fclose(FILE);

    FILE = fopen("run_Relax_"+name+"_Generation"+i+"_20res.sbatch", 'w');
    fprintf(FILE, "#!/bin/bash\n");
    fprintf(FILE, "#SBATCH --nodes=1\n");
    fprintf(FILE, "#SBATCH --tasks-per-node=24\n");
    fprintf(FILE, "#SBATCH --mem=30GB\n");
    fprintf(FILE, "#SBATCH --time=24:00:00\n");
    fprintf(FILE, "#SBATCH --job-name=Generate"+i+"\n");
    fprintf(FILE, "#SBATCH --mail-user=qz886@nyu.edu\n");
    fprintf(FILE, "module purge\n");
    fprintf(FILE, "module rosetta/openmpi/intel/2020.46.61480\n");
    fprintf(FILE, "srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @relax_"+name+"_Generation"+i+"_20res.flags\n");
    fclose(FILE);

    [~, job_id] = system("sbatch run_Relax_"+name+"_Generation"+i+"_20res.sbatch");
    job_id = strsplit(job_id);
    job_id = job_id(4);
    while true
        [~, job_status] = system("squeue -j "+job_id);
        disp(job_status);
        job_status = strsplit(job_status);
        if length(job_status) > 10
            jtime = job_status(15);
            disp("current status "+jtime);
        else
            break
        end
        pause(20);
    end
    status = strsplit(jtime{1},":");
    minutes = minutes + str2double(status(1));
    seconds = seconds + str2double(status(2));

    tic;
    FILE = fopen("score_"+name+"_Generation"+i+"_20res.sc", 'r');
    data = textscan(FILE, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s', 'HeaderLines', 2);
    Scores = data{2}.';
    RMSD = data{23}.';
    Descriptions = data{25};
    fclose(FILE);
    [Coordinates, Angles, Numbers] = PDB2Coordinates(Descriptions);
    Labels = Labels(Numbers);
    disp("Generation"+i+" has a population of "+length(Scores)+" after relaxation. Of which "+length(Labels(Labels==1))+" are xovers, and "+length(Labels(Labels==2))+" are mutants.");

    cyc = zeros(1,size(Coordinates,2)/3);
    for ci = 1 : size(Coordinates,2)/3
        coordinates = Coordinates(:,3*ci-2:3*ci);
        cyc(ci) = norm(coordinates(end-1,:)-coordinates(1,:));
    end
    disp("Max cyclic error is "+max(cyc)+".");

    % Clustering to select lowest energy centers for next generation
    [Scores, bestIndices] = mink(Scores, length(Scores));
    Angles = Angles(:, bestIndices);
    RMSD = RMSD(bestIndices);
    Labels = Labels(bestIndices);
    N = size(Coordinates,2)/3;
    Coordinates = reshape(Coordinates, [4*n*3, N]);
    Coordinates = Coordinates(:, bestIndices);

    indices = clustering(Coordinates);
    Angles = Angles(:, indices);
    Scores = Scores(indices);
    RMSD = RMSD(indices);
    Labels = Labels(indices);
    Coordinates = Coordinates(:, indices);
    disp("Generation"+i+" has "+length(Scores)+" centers after clustering. Of which "+length(Labels(Labels==1))+" are xovers, and "+length(Labels(Labels==2))+" are mutants.");
    
    temp_len = length(candScores);
    candAngles = [candAngles, Angles(:,Scores<score_threshold)];
    candScores = [candScores, Scores(Scores<score_threshold)];
    candRMSD = [candRMSD, RMSD(Scores<score_threshold)];
    candCoordinates = [candCoordinates, Coordinates(:,Scores<score_threshold)];
    [candAngles, ia, ~] = unique(candAngles.', 'rows');
    candAngles = candAngles.';
    candScores = candScores(ia);
    candRMSD = candRMSD(ia);
    candCoordinates = candCoordinates(:,ia);
    pop_incre = length(candScores)-temp_len;
    disp("Generation"+i+" produces "+pop_incre+" new candidates.")

    % Of these centers, select the top ones with lowest Energy
    [Scores, bestIndices] = mink(Scores, population);
    Angles = Angles(:, bestIndices);
    RMSD = RMSD(bestIndices);
    Labels = Labels(bestIndices);
    Coordinates = Coordinates(:, bestIndices);
    N = size(Coordinates,2);
    Coordinates = reshape(Coordinates, [4*n, 3*N]);
    disp("Of top centers with lowest RMSD, "+length(Labels(Labels==1))+" are xovers, and "+length(Labels(Labels==2))+" are mutants.");
    disp("Min RMSD is "+min(RMSD)+", max RMSD is "+max(RMSD)+", and mean RMSD is "+mean(RMSD)+".");
    disp("Min Energy is "+min(Scores)+", max Energy is "+max(Scores)+", and mean Energy is "+mean(Scores)+".");
    seconds = seconds + toc;

    if pop_incre == 0
        stall = stall + 1;
    else
        stall = 0;
    end

    system("rm "+name+"_Generation"+i+"*_20res*.pdb");
    system("rm "+name+"_Generation"+i+"_list_20res.txt");
    system("rm relax_"+name+"_Generation"+i+"_20res.flags");
    system("rm run_Relax_"+name+"_Generation"+i+"_20res.sbatch");
    system("rm score_"+name+"_Generation"+i+"_20res.sc");
    system("rm slurm-"+job_id+".out");

    filename = sprintf('After_GA_%s_20res.mat', name);
    save(filename, 'Angles', 'Coordinates', 'RMSD', 'Scores', 'candAngles', 'candScores', 'candRMSD', 'candCoordinates');

    i = i + 1;
    rtime = minutes + ceil(seconds/60);
    disp("Running time now "+rtime+" minutes.");
    population = population - 5;
end
end


% Angles (size 3*n x N) and Coordinates (4*n x 3*N) form the initial population.
% Try all possible cross-overs using Kabsch algorithm, and rmsd decides if 
% succeeds. Random mutations compatible with the sequence.
function [Coordinates_new, Labels] = GA(Angles, Coordinates, RamaMap, xover_region, mut_region, xover_size, mut_size)

n = 20;
N = size(Angles,2); % population size
rmsd_threshold = 0.5;
e_threshold = 0.3;
Coordinates_new = [];
Labels = []; % 1 for xover, 2 for mutant

% Crossover
combinations = nchoosek(1:N, 2);
parents = randperm(size(combinations,1), size(combinations,1));
batchsize = ceil(xover_size/3);
batches = ceil(size(combinations,1)/batchsize);

for bat = 0 : batches-1
    bstart = bat*batchsize+1;
    if bat < batches-1
        bend = (bat+1)*batchsize;
    else
        bend = size(combinations,1);
    end
    parfor pair = bstart : bend
        bCoordinates = [];
        bLabels = [];
        
        p1 = combinations(parents(pair),1);
        p2 = combinations(parents(pair),2);
        xover_length = randi(xover_region);
        break1 = randi([0, n-1]);
        break2 = break1+xover_length;
        if break2 > n
            break2 = break2 - n;
        end

        coordinates1 = Coordinates(:, 3*(p1-1)+1:3*p1);
        coordinates2 = Coordinates(:, 3*(p2-1)+1:3*p2);
        [success, coordinates_x1, coordinates_x2] = xover(coordinates1, coordinates2, break1, break2, rmsd_threshold);

        if success
            bCoordinates = [coordinates_x1, coordinates_x2];
            bLabels = [1, 1];
        end

        Coordinates_new = [Coordinates_new, bCoordinates];
        Labels = [Labels, bLabels];
    end
    if length(Labels) >= xover_size
        break;
    end
end

% Mutation
mutants = randsample(N, 100*N, true);
batchsize = ceil(mut_size/3);
batches = ceil(100*N/batchsize);

for bat = 0 : batches-1
    bstart = bat*batchsize+1;
    if bat < batches-1
        bend = (bat+1)*batchsize;
    else
        bend = 100*N;
    end
    parfor choice = bstart : bend
        mut_length = randi(mut_region);
        mbreak = randi([0, n-1]);
        angles = Angles(:, mutants(choice));
        coordinates = Coordinates(:, 3*mutants(choice)-2:3*mutants(choice));

        NCa_align = norm(coordinates(2,:)-coordinates(1,:));
        CaC_align = norm(coordinates(3,:)-coordinates(2,:));
        NR_align = acos(dot((coordinates(1,:)-coordinates(end-1,:))/norm(coordinates(1,:)-coordinates(end-1,:)), (coordinates(1,:)-coordinates(2,:))/norm(coordinates(1,:)-coordinates(2,:))));
        CaR_align = acos(dot((coordinates(2,:)-coordinates(1,:))/norm(coordinates(2,:)-coordinates(1,:)), (coordinates(2,:)-coordinates(3,:))/norm(coordinates(2,:)-coordinates(3,:))));

        N_align = [0;0;0];
        rotation_align = T(NR_align)*R(angles(1));
        Ca_align = N_align + rotation_align*[NCa_align; 0; 0];
        rotation_align = rotation_align*T(CaR_align)*R(angles(2));
        C_align = Ca_align + rotation_align*[CaC_align; 0; 0];

        P = coordinates(1:3,:);
        Q = [N_align.'; Ca_align.'; C_align.'];
        [Rotation, translation, ~] = kabsch_algorithm(P, Q);
        coordinates = (Rotation * coordinates')' + translation';

        k0 = 10;
        break1 = mbreak;
        break2 = break1+mut_length;
        if break2 > n
            break2 = break2-n;
        end
        N_start = coordinates(break1*4+1,:).';
        if break2 == n
            N_end = coordinates(1,:).';
            Ca_end = coordinates(2,:).';
        else
            N_end = coordinates(break2*4+1,:).';
            Ca_end = coordinates(break2*4+2,:).';
        end

        if break1 >= break2
            angles_input = [angles; angles(1:3*break2)];
        else
            angles_input = angles(1:3*break2);
        end
        [success, mut_part, ~] = mutate(angles_input, N_start, N_end, Ca_end, break1, mut_length, e_threshold, RamaMap, k0);
        
        if break1 >= break2
            coordinates_mut = [mut_part(4*(n-break1)+1:end,:); coordinates(4*break2+1:4*break1,:); mut_part(1:4*(n-break1),:)];
        else
            coordinates_mut = [coordinates(1:4*break1,:); mut_part; coordinates(4*break2+1:end,:)];
        end

        if success
            Coordinates_new = [Coordinates_new, coordinates_mut];
            Labels = [Labels, 2];
        end
    end
    if length(Labels) >= xover_size+mut_size
        break;
    end
end
end


% Exchange parts of the peptides between break points, residues break1+1 to
% break2 of peptide1 and peptide2 are swaped. Using Kabsch algorithm, we
% align atoms Ca and C of break1, N and Ca of break1+1, Ca and C of break2, 
% N and Ca of break2+1, from the two peptides, to connect the parts. 
function [success, coordinates_x1, coordinates_x2] = xover(coordinates1, coordinates2, break1, break2, rmsd_threshold)

n = 20;

if break1 == 0
    atom_pos = [n*4-2,n*4-1,break1*4+1,break1*4+2,break2*4-2,break2*4-1,break2*4+1,break2*4+2];
elseif break2 == n
    atom_pos = [break1*4-2,break1*4-1,break1*4+1,break1*4+2,break2*4-2,break2*4-1,1,2];
else
    atom_pos = [break1*4-2,break1*4-1,break1*4+1,break1*4+2,break2*4-2,break2*4-1,break2*4+1,break2*4+2];
end

P = coordinates2(atom_pos,:);
Q = coordinates1(atom_pos,:);
[Rotation, translation, rmsd] = kabsch_algorithm(P, Q);
coordinates2_aligned = (Rotation * coordinates2')' + translation';
if break1 <= break2
    coordinates_x1 = [coordinates1(1:break1*4,:); coordinates2_aligned(break1*4+1:break2*4,:); coordinates1(break2*4+1:n*4,:)];
    coordinates_x2 = [coordinates2_aligned(1:break1*4,:); coordinates1(break1*4+1:break2*4,:); coordinates2_aligned(break2*4+1:n*4,:)];
else
    coordinates_x1 = [coordinates2_aligned(1:break2*4,:); coordinates1(break2*4+1:break1*4,:); coordinates2_aligned(break1*4+1:end,:)];
    coordinates_x2 = [coordinates1(1:break2*4,:); coordinates2_aligned(break2*4+1:break1*4,:); coordinates1(break1*4+1:end,:)];
end

if rmsd <= rmsd_threshold
    success = true;
else
    success = false;
end
end


% Mutate consecutive parts of the peptides between break points, torsion
% angles of residues break1+1 to break2. Use simulated annealing to make
% sure parts connect at break points.
function [success, mut_part, mut_part_angles] = mutate(angles, N_start, N_end, Ca_end, mbreak, mut_length, e_threshold, RamaMap, k0)
n = 20;
break1 = mbreak;
break2 = break1+mut_length;
angles_mut = angles;
perb = deg2rad(2*k0)*rand(2,1)-deg2rad(k0);
angles_mut(break1*3+1:break1*3+2) = angles_mut(break1*3+1:break1*3+2)+perb;
[error, mut_part] = N_error(N_start, N_end, Ca_end, angles_mut, mbreak, mut_length);

M = 1000;
t0_error = 1;
c_error = 20;
b = k0-1;

i = 1;
while error > e_threshold && i <= M
    k = k0/(1+b*i/M);
    perb = zeros(3,break2);
    perb(1:2,break1+2:break2) = deg2rad(2*k)*rand(2,mut_length-1)-deg2rad(k);
    perb = reshape(perb, [3*break2,1]);

    temp_angle = angles_mut+perb;
    temp_angle(temp_angle<-pi) = temp_angle(temp_angle<-pi)+2*pi;
    temp_angle(temp_angle>pi) = temp_angle(temp_angle>pi)-2*pi;
    locations = floor(rad2deg(temp_angle)/10) + 19;
    for res = break1+2 : break2
        if res > n
            resi = res-n;
        else
            resi = res;
        end
        if RamaMap(locations(3*res-2), locations(3*res-1), resi) ~= 1
            perb(3*res-2:3*res-1) = 0;
        end
    end    

    angles_mut_new = angles_mut+perb;
    [error_new, mut_part_new] = N_error(N_start, N_end, Ca_end, angles_mut_new, mbreak, mut_length);
    temp_error = t0_error/(1+c_error*i/M);
    if error_new <= error
        angles_mut = angles_mut_new;
        error = error_new;
        mut_part = mut_part_new;
    else
        prob = exp(1)^((error-error_new)/temp_error);
        if rand <= prob
            angles_mut = angles_mut_new;
            error = error_new;
            mut_part = mut_part_new;
        end
    end
    i = i+1;
end

if error <= e_threshold
    success = true;
else
    success = false;
end

mut_part_angles = angles_mut(break1*3+1:break2*3);
end


function [error, mut_part] = N_error(N_start, N_end, Ca_end, angles, mbreak, mut_length)
% bond angles
NR = deg2rad(121.7);
CaR = deg2rad(111.2);
CR = deg2rad(116.2);
COR = deg2rad(123.0);

% bond lengths
NCa = 1.458;
CaC = 1.524;
CN = 1.329;
CO = 1.231;

rotation = eye(3);
for i = 1 : mbreak
    rotation = rotation*T(NR)*R(angles(3*i-2));
    rotation = rotation*T(CaR)*R(angles(3*i-1));
    rotation = rotation*T(CR)*R(angles(3*i));
end

N = N_start;
rotation = rotation*T(NR)*R(angles(3*mbreak+1));
Ca = N + rotation*[NCa; 0; 0];
mut_part = zeros(4*mut_length,3);
mut_part(1,:) = N.';
mut_part(2,:) = Ca.';

for i = mbreak+1 : mbreak+mut_length
    rotation = rotation*T(CaR)*R(angles(3*i-1));
    C = Ca + rotation*[CaC; 0; 0];
    rotation = rotation*T(CR)*R(angles(3*i));
    O = C + rotation*[-CO*cos(pi-COR); CO*sin(pi-COR); 0];
    N = C + rotation*[CN; 0; 0];
    if i < mbreak+mut_length
        rotation = rotation*T(NR)*R(angles(3*i+1));
        Ca = N + rotation*[NCa; 0; 0];
    else
        rotation = rotation*T(NR);
        Ca = N + rotation*[NCa; 0; 0];
    end

    mut_part(4*(i-mbreak)-1,:) = C.';
    mut_part(4*(i-mbreak),:) = O.';
    if i < mbreak+mut_length
        mut_part(4*(i-mbreak)+1,:) = N.';
        mut_part(4*(i-mbreak)+2,:) = Ca.';
    end
end

error = norm(N-N_end)+norm(Ca-Ca_end);
end


% KABSCH_ALGORITHM calculates the optimal rigid body transformation
% that aligns two sets of 3D points (P and Q) using the Kabsch algorithm.
%
% Inputs:
%   P: Nx3 array of points to be aligned
%   Q: Nx3 array of target points
%
% Outputs:
%   R: 3x3 rotation matrix
%   t: 3x1 translation vector
%   rmsd: root-mean-square deviation between the aligned points
function [R, t, rmsd] = kabsch_algorithm(P, Q)

% Calculate the centroids of the two sets of points
centroid_P = mean(P, 1);
centroid_Q = mean(Q, 1);

% Center the points by subtracting their centroids
P_centered = P - centroid_P;
Q_centered = Q - centroid_Q;

% Calculate the covariance matrix of the centered points
covariance_matrix = P_centered' * Q_centered;

% Calculate the optimal rotation matrix using singular value decomposition (SVD)
[U, ~, V] = svd(covariance_matrix);
rotation_matrix = V * U';

% If the determinant of the rotation matrix is negative, we need to flip one axis
if det(rotation_matrix) < 0
    V(:, 3) = -V(:, 3);
    rotation_matrix = V * U';
end

% Calculate the translation vector
translation_vector = centroid_Q' - rotation_matrix * centroid_P';

% Apply the rotation and translation to the original set of points
P_aligned = (rotation_matrix * P')' + translation_vector';

% Calculate the root-mean-square deviation (RMSD) between the aligned points
rmsd = sqrt(sum(sum((Q - P_aligned).^2)) / size(P, 1));

% Output the rotation matrix, translation vector, and RMSD
R = rotation_matrix;
t = translation_vector;
end


function Coordinates2PDB(Coordinates, name, seq, generation)
n = 20;
seq = seq.upper;
d_aa_codes = containers.Map({'DALA', 'DARG', 'DASN', 'DASP', 'DGLU', 'DGLN', 'DHIS', 'DILE', 'DLEU', 'DLYS', 'DMET', 'DPHE', 'DPRO', 'DSER', 'DTHR', 'DTRP', 'DTYR', 'DVAL'}, ...
                            {'DAL', 'DAR', 'DAN', 'DAS', 'DGU', 'DGN', 'DHI', 'DIL', 'DLE', 'DLY', 'DME', 'DPH', 'DPR', 'DSE', 'DTH', 'DTR', 'DTY', 'DVA'});
for res = 1 : n
    amiaci = char(seq(res));
    if amiaci(1) == 'D'
        seq(res) = d_aa_codes(seq(res));
    end
end

atom_names = repmat(["N", "CA", "C", "O"], 1, n);
atom_number = 1:4*n;
res_names = repelem(seq, 4);
res_number = repelem([1:n], 4);
elements = repmat(["N", "C", "C", "O"], 1, n);

N = size(Coordinates,2)/3;
for cand = 1 : N
    peptide.X = Coordinates(:,3*cand-2);
    peptide.Y = Coordinates(:,3*cand-1);
    peptide.Z = Coordinates(:,3*cand);
    filename = sprintf(name+"_Generation%d_cand%d_20res.pdb", generation, cand);
    peptide.outfile = filename;
    peptide.atomName = atom_names;
    peptide.atomNum = atom_number;
    peptide.resName = res_names;
    peptide.resNum = res_number;
    peptide.element = elements;

    FILE = fopen(filename, 'w');
    fprintf(FILE, 'SEQRES   1 A  20\n');
    fclose(FILE);
    mat2pdb(peptide);
end
end


function [Coordinates, Angles, Numbers] = PDB2Coordinates(Descriptions)
n = 20;
Coordinates = zeros(4*n, 3*length(Descriptions));
Angles = zeros(3*n, length(Descriptions));
Numbers = zeros(1, length(Descriptions));

for cand = 1 : length(Descriptions)
    dname = strsplit(Descriptions{cand}, '_');
    dname = strsplit(dname{5}, 'cand');
    dname = str2double(dname{2});
    Numbers(cand) = dname;

    filename = Descriptions(cand)+".pdb";
    peptide = pdb2mat(filename);
    coordinate = [peptide.X.', peptide.Y.', peptide.Z.'];
    atomName = peptide.atomName;
    coorder = zeros(4*n,3);

    N_coor = coordinate(atomName=="N",:);
    Ca_coor = coordinate(atomName=="CA",:);
    C_coor = coordinate(atomName=="C",:);
    O_coor = coordinate(atomName=="O",:);

    for res = 1 : n
        coorder(4*res-3,:) = N_coor(res,:);
        coorder(4*res-2,:) = Ca_coor(res,:);
        coorder(4*res-1,:) = C_coor(res,:);
        coorder(4*res,:) = O_coor(res,:);
    end

    angle = Coordinates2Angles(coorder);
    Coordinates(:,3*cand-2:3*cand) = coorder;
    Angles(:,cand) = angle;
end
end


function angle = Coordinates2Angles(coor)
n = 20;
angle = zeros(3*n,1);
% atoms N, CA, C, O
for res = 2 : n-1
    phi = torsion(coor(4*res-5,:), coor(4*res-3,:), coor(4*res-2,:), coor(4*res-1,:));
    psi = torsion(coor(4*res-3,:), coor(4*res-2,:), coor(4*res-1,:), coor(4*res+1,:));
    omega = torsion(coor(4*res-2,:), coor(4*res-1,:), coor(4*res+1,:), coor(4*res+2,:));
    angle(3*res-2:3*res) = [phi; psi; omega];
end
phi = torsion(coor(4*n-1,:), coor(1,:), coor(2,:), coor(3,:));
psi = torsion(coor(1,:), coor(2,:), coor(3,:), coor(5,:));
omega = torsion(coor(2,:), coor(3,:), coor(5,:), coor(6,:));
angle(1:3) = [phi; psi; omega];
phi = torsion(coor(4*n-5,:), coor(4*n-3,:), coor(4*n-2,:), coor(4*n-1,:));
psi = torsion(coor(4*n-3,:), coor(4*n-2,:), coor(4*n-1,:), coor(1,:));
omega = torsion(coor(4*n-2,:), coor(4*n-1,:), coor(1,:), coor(2,:));
angle(3*n-2:3*n) = [phi; psi; omega];
angle = -angle;
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


function indices = clustering(Coordinates)
n = 20;
N = size(Coordinates,2);
lib = 1 : N;
indices = [];

while ~isempty(lib)
    members = [];
    center = lib(1);
    centercoor = Coordinates(:,center);
    parfor i = 1 : length(lib)
        coor = Coordinates(:,lib(i));
        [~, ~, rmsd] = kabsch_algorithm(reshape(coor, [4*n,3]), reshape(centercoor, [4*n,3]));
        if rmsd < 0.5
            members = [members, lib(i)];
        end
    end
    lib = setdiff(lib, members);
    indices = [indices, center];
end
end


function coordinates = get_coordinates(angles)
% bond angles
NR = deg2rad(121.7);
CaR = deg2rad(111.2);
CR = deg2rad(116.2);
COR = deg2rad(123.0);

% bond lengths
NCa = 1.458;
CaC = 1.524;
CN = 1.329;
CO = 1.231;

n = length(angles)/3;
coordinates = zeros(4*n, 3);
N = [0;0;0];
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n
    rotation = rotation*T(NR)*R(angles(3*i-2));
    Ca = N + rotation*[NCa; 0; 0];
    rotation = rotation*T(CaR)*R(angles(3*i-1));
    C = Ca + rotation*[CaC; 0; 0];
    rotation = rotation*T(CR)*R(angles(3*i));
    O = C + rotation*[-CO*cos(pi-COR); CO*sin(pi-COR); 0];
    N = C + rotation*[CN; 0; 0];

    coordinates(4*i-2,:) = Ca.';
    coordinates(4*i-1,:) = C.';
    coordinates(4*i,:) = O.';
    if i < n
        coordinates(4*i+1,:) = N.';
    end
end
end


function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end


function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end
