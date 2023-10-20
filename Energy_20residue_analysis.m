function Energy_20residue_analysis

GoodAngles = [];
GoodScore = [];
GoodCount = [];

for bat = 1 : 101
    filename = sprintf('GoodAngles20res_%d.mat', bat);
    file = load(filename);
    Angles = file.GoodAngles;
    Rama = file.GoodRama;
    Rep = file.GoodRep;
    Hbond = file.GoodHbond;
    Count = file.GoodCount;

    Score = 0.45*Rama + Rep + Hbond;

    GoodAngles = [GoodAngles, Angles];
    GoodScore = [GoodScore, Score];
    GoodCount = [GoodCount, Count];
end

N_best = length(GoodScore);
[bestScores, bestIndices] = mink(GoodScore, N_best);
bestAngles = GoodAngles(:, bestIndices);
bestCount = GoodCount(bestIndices);

save('bestAngles_20res.mat', 'bestAngles', 'bestScores', 'bestCount');
end