with open('/scratch/qz886/Clustercenters_7res/Pnear_list.txt') as f:
    lines = f.readlines()
candidates = []
scores = []
for l in lines:
    if l != "\n":
        candidates.append(l.split()[0])
        scores.append(l.split()[1])
        
with open('slurm-31708812.out', 'r') as f:
    plines = f.readlines()    
Pnear = []
for l in plines:
    if l[0:5] == 'PNear':
        Pnear.append(float(l.split()[1]))

Hbonds = []
for cand in candidates:
    with open('output_'+cand+'.pdb', 'r') as f:
        lines = f.readlines()
    for l in lines:
        if 'min_internal_hbonds' in l:
            Hbonds.append(float(l.split()[1]))

with open('Pnear_values_7res.txt', 'w') as f:
    for i in range(len(candidates)):
        f.write(candidates[i]+'\t'+scores[i]+f'\t{Pnear[i]:.3f}\t{Hbonds[i]:.0f}\n')
