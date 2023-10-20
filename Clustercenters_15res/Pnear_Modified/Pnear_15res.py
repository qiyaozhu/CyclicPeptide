# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:00:01 2023

@author: Cathe
"""

import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.io
import re
from matplotlib import rc
from scipy.stats import pearsonr

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
    

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 30}

aacode = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

with open('/scratch/qz886/Clustercenters_15res/Pnear_list.txt') as f:
    lines = f.readlines()
    
candidates = []
scores = []
seqs = []

for l in lines:
    if l != "\n":
        candidates.append(l.split()[0])
        scores.append(float(l.split()[1]))
        
        seq = l.split()[2]
        res1 = seq.split('[')[1].split(':')[0]
        resn = seq.split('[')[-1].split(':')[0]
        if res1 == 'HIS_D':
            res1 = 'DHIS'
        if resn == 'HIS_D':
            resn = 'DHIS'
        residues = [res1]

        lpos = [m.start() for m in re.finditer('\[', seq)]
        rpos = [m.start() for m in re.finditer('\]', seq)]
        seq = seq[rpos[0]+1:lpos[-1]-1]
        
        p = 0
        while p < len(seq):
            if p == len(seq)-1:
                residues.append(aacode[seq[p]])
                p += 1
            else:
                if seq[p+1] != '[':
                    residues.append(aacode[seq[p]])
                    p += 1
                else:
                    if seq[p+2:p+6] == 'HIS_':
                        residues.append('DHIS')
                        p += 8
                    else:
                        residues.append(seq[p+2:p+6])
                        p += 7
        residues.append(resn)
        
        s = residues[0]
        for r in range(1,len(residues)):
            s = s+'-'+residues[r]        
        seqs.append(s)



lamda = 1.5
kbT = 0.62
Pnear = []

with open('/scratch/qz886/Clustercenters_15res/Pnear_Candidates/Pnear_Cartesian_15res_1.out', 'r') as f:
    plines = f.readlines()
with open('/scratch/qz886/Clustercenters_15res/Pnear_Candidates/Pnear_Cartesian_15res_2.out', 'r') as f:
    plines2 = f.readlines()
with open('/scratch/qz886/Clustercenters_15res/Pnear_Candidates/Pnear_Cartesian_15res_3.out', 'r') as f:
    plines3 = f.readlines()
plines += plines2
plines += plines3
    
l = 0
while l < len(plines):
    if len(plines[l].split()) > 0:
        if plines[l].split()[0] == "MPI_worker_node":
            l += 1            
            Energy = []
            RMSD = []
            while plines[l].split()[0] != "End":
                Energy.append(float(plines[l].split()[3]))
                RMSD.append(float(plines[l].split()[2]))
                l += 1
            p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
            Pnear.append(p)
        else:
            l += 1
    else:
        l += 1



PnearM = []
for i in range(len(candidates)):
    path = 'After_GA_'+candidates[i]+'_15res.mat'
    if os.path.exists(path):
        mat = scipy.io.loadmat(path)
        Energy = mat["candScores"][0]
        RMSD = mat["candRMSD"][0]
        p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
        PnearM.append(p)

with open('Pnear_values_15res.txt', 'w') as f:
    f.write('Name\t\t\t\tEnergy\t\tRosetta\t\tGA\t\tSequence\n')
    for i in range(len(PnearM)):
        if Pnear[i] > 0.9:
            f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t\t{PnearM[i]:.3f}\t\t'+seqs[i]+'\n')
    f.write('\nName\t\t\t\tEnergy\t\tRosetta\t\tGA\t\tSequence\n')
    for i in range(len(PnearM)):
        if PnearM[i] > 0.9:
            f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t\t{PnearM[i]:.3f}\t\t'+seqs[i]+'\n')
    f.write('\nName\t\t\t\tEnergy\t\tRosetta\t\tGA\t\tSequence\n')
    for i in range(len(PnearM)):
        if abs(Pnear[i]-PnearM[i]) > 0.4:
            f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t\t{PnearM[i]:.3f}\t\t'+seqs[i]+'\n')
