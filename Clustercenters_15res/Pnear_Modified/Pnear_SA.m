function Pnear_SA

FILE = fopen("Pnear_seq.txt", 'r');
data = textscan(FILE, '%s\t%s\n');
fclose(FILE);
CenterNames = data{2};

Names = {};
for i = 13 : 100
    name = CenterNames(i);
    name = name{1};
    name = split(name, "_");
    name = name{1};
    name = split(name, "clustercenter");
    name = name{2};
    Names{i-12} = name;
end

disp("Start jobs ...");
[~, job_id1] = system("sbatch run_Pnear_sampling_SA_"+Names(1)+".sbatch");
job_id1 = strsplit(job_id1);
job_id1 = job_id1(4);
disp("Job 1 for "+Names(1)+" submitted!");

[~, job_id2] = system("sbatch run_Pnear_sampling_SA_"+Names(2)+".sbatch");
job_id2 = strsplit(job_id2);
job_id2 = job_id2(4);
disp("Job 2 for "+Names(2)+" submitted!");

[~, job_id3] = system("sbatch run_Pnear_sampling_SA_"+Names(3)+".sbatch");
job_id3 = strsplit(job_id3);
job_id3 = job_id3(4);
disp("Job 3 for "+Names(3)+" submitted!");

cand = 4;
while cand <= length(Names)
    [~, job_status1] = system("squeue -j "+job_id1);
    job_status1 = strsplit(job_status1);
    if length(job_status1) > 10
        jtime = job_status1(15);
        disp("job 1 current status "+jtime);
    else
        [~, job_id1] = system("sbatch run_Pnear_sampling_SA_"+Names(cand)+".sbatch");
        job_id1 = strsplit(job_id1);
        job_id1 = job_id1(4);
        cand = cand + 1;
        disp("Job 1 for "+Names(cand)+" submitted!");
    end

    [~, job_status2] = system("squeue -j "+job_id2);
    job_status2 = strsplit(job_status2);
    if length(job_status2) > 10
        jtime = job_status2(15);
        disp("job 2 current status "+jtime);
    else
        [~, job_id2] = system("sbatch run_Pnear_sampling_SA_"+Names(cand)+".sbatch");
        job_id2 = strsplit(job_id2);
        job_id2 = job_id2(4);
        cand = cand + 1;
        disp("Job 2 for "+Names(cand)+" submitted!");
    end

    [~, job_status3] = system("squeue -j "+job_id3);
    job_status3 = strsplit(job_status3);
    if length(job_status3) > 10
        jtime = job_status3(15);
        disp("job 3 current status "+jtime);
    else
        [~, job_id3] = system("sbatch run_Pnear_sampling_SA_"+Names(cand)+".sbatch");
        job_id3 = strsplit(job_id3);
        job_id3 = job_id3(4);
        cand = cand + 1;
        disp("Job 3 for "+Names(cand)+" submitted!");
    end

    pause(30);
end
end