clc;
clear all;

GoodAngles = [];
GoodScore = [];

for bat = 1 : 40
    filename = sprintf('GoodAngles7res_%d.mat', bat);
    file = load(filename);
    Angles = file.GoodAngles;
    Rama = file.GoodRama;
    Rep = file.GoodRep;
    Hbond = file.GoodHbond;
    Energy = file.GoodEnergy;

    Score = 0.45*Rama + Rep + Hbond + Energy;

    GoodAngles = [GoodAngles, Angles];
    GoodScore = [GoodScore, Score];
end

N_best = length(GoodScore);
[bestScores, bestIndices] = mink(GoodScore, N_best);
bestAngles = GoodAngles(:, bestIndices);

save('bestAngles_7res.mat', 'bestAngles', 'bestScores');