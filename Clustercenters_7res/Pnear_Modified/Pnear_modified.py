# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 12:56:45 2023

@author: Cathe
"""
import numpy as np

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

with open('Pnear_seq.txt', 'r') as f:
    lines = f.readlines()
candidates = []
for l in lines:
    if l != "\n":
        candidates.append(l.split()[1])

with open('Pnear_values.txt', 'w') as f1:
    f1.write('')

for cand in candidates:
    path = 'score_'+cand+'.sc'
    with open(path, 'r') as f:
        lines = f.readlines()
    Energy = []
    RMSD = []    
    for l in lines:
        if len(l.split()) > 1:
            if is_float(l.split()[1]):
                Energy.append(float(l.split()[1]))
                RMSD.append(float(l.split()[25]))
    
    lamda = 0.5
    kbT = 1
    p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))

    with open('Pnear_values.txt', 'a') as f1:
        f1.write(cand+f'\t{p:.3f}\n')
        
