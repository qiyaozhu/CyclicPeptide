# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

with open('/scratch/qz886/design_candidates.txt') as f:
    lines = f.readlines()
    
candidates = []

for l in lines:
    if l != '\n':
        candidates.append(l.split('.')[0])
        
for cand in candidates:
    with open(cand+'.pdb', 'r') as f:
        lines = f.readlines()    
    res = []
    
    for i in range(len(lines)):
        if lines[i].split()[0] == 'pose':
            res.append(lines[i+1].split(':')[0])
            for j in range(2,7):
                res.append(lines[i+j].split('_')[0])
            res.append(lines[i+7].split(':')[0])
            break
        
    with open(cand+'_seq.txt', 'w') as f:
        for r in res:
            f.write(r+' ')
