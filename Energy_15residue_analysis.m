clc;
clear all;

GoodAngles = [];
GoodScore = [];

for bat = 1 : 101
    filename = sprintf('GoodAngles15res_%d.mat', bat);
    file = load(filename);
    Angles = file.GoodAngles;
    Rama = file.GoodRama;
    Rep = file.GoodRep;
    Hbond = file.GoodHbond;

    Score = 0.45*Rama + Rep + Hbond;

    GoodAngles = [GoodAngles, Angles];
    GoodScore = [GoodScore, Score];
end

N_best = length(GoodScore);
[bestScores, bestIndices] = mink(GoodScore, N_best);
bestAngles = GoodAngles(:, bestIndices);

save('bestAngles_15res.mat', 'bestAngles', 'bestScores');
