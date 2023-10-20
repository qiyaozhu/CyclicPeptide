# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:23:02 2023

@author: Cathe
"""
import os
import re

with open('FastDesign_silent.txt', 'r') as f:
    lines = f.readlines()

os.system('mkdir Pnear')
    
Scores = []
Names = []
Seqs = []

i = 2
while i < len(lines):
    if lines[i].split()[0] == 'SCORE:':
        if float(lines[i].split()[1]) < -40:
            scoreline = lines[i]
            Scores.append(float(scoreline.split()[1]))
            name = scoreline.split()[-1]
            Names.append(name)
            i += 1

            while lines[i].split()[0] != 'SCORE:':
                if lines[i].split()[0] == 'REMARK':
                    i += 1
                else:
                    if lines[i].split()[0] == 'ANNOTATED_SEQUENCE:':
                        seq = lines[i].split()[1]
                        Seqs.append(seq)
                        seqcode = seq[0]
                        lpos = [m.start() for m in re.finditer('\[', seq)]
                        rpos = [m.start() for m in re.finditer('\]', seq)]
                        for pp in range(len(rpos)-1):
                            seqcode += seq[rpos[pp]+1 : lpos[pp+1]]
                        
                        if lines[i+1].split()[0] == 'NONCANONICAL_CONNECTION:':
                            with open('Pnear/Pnear_'+name+'.txt', 'w') as f:
                                f.write('SEQUENCE: '+seqcode+'\n')
                                f.write(lines[1])
                                f.write('REMARK BINARY SILENTFILE\n')
                                f.write(scoreline)
                                f.write(lines[i])
                                f.write(lines[i+1])
                            i += 2
                        else:
                            with open('Pnear/Pnear_'+name+'.txt', 'w') as f:
                                f.write('SEQUENCE: '+seqcode+'\n')
                                f.write(lines[1])
                                f.write(scoreline)
                                f.write('REMARK PROTEIN SILENTFILE\n')
                                f.write(lines[i])
                            i += 1
    
                    with open('Pnear/Pnear_'+name+'.txt', 'a') as f:
                        f.write(lines[i])
                    i += 1
                    if i >= len(lines):
                        break
        else:
            i +=1 
    else:
        i += 1
        
Survivors = {}
for i in range(len(Seqs)):
    if Seqs[i] in Survivors:
        if Scores[i] < Survivors[Seqs[i]][0]:
            Survivors[Seqs[i]] = [Scores[i], Names[i]]
    else:
        Survivors[Seqs[i]] = [Scores[i], Names[i]]

Survivors = {k:v for k, v in sorted(Survivors.items(), key=lambda item: item[1])}

with open('Pnear_list.txt', 'w') as f:
    for s in Survivors:
        f.write(Survivors[s][1]+f'\t{Survivors[s][0]:.3f}\t'+s+'\n')

os.system('mkdir Pnear_Candidates')

for s in Survivors:
    name = Survivors[s][1]
    os.system('mv Pnear/Pnear_'+name+'.txt Pnear_Candidates/Pnear_'+name+'.txt')
os.system('rm -r Pnear')
