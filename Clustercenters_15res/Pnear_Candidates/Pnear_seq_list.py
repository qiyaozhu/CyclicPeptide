# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import re

with open('/scratch/qz886/Clustercenters_15res/Pnear_list.txt') as f:
    lines = f.readlines()
    
candidates = []
seqs = []

for l in lines:
    if l != '\n':
        candidates.append(l.split()[0])
        seqs.append(l.split()[2])

aacode = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
      
for i in range(len(candidates)):
    seq = seqs[i]
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
    
    with open('Pnear_'+candidates[i]+'_seq.txt', 'w') as f:
        for r in residues:
            f.write(r+' ')
