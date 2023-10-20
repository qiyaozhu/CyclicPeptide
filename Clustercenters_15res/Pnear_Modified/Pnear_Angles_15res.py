# -*- coding: utf-8 -*-
"""
Created on Wed May 10 19:16:48 2023

@author: Cathe
"""

import os
import time
import subprocess
import numpy as np
import sys


Names = []
Seqs = []
with open("Pnear_seq.txt", 'r') as f:
    lines = f.readlines()
for i in range(12,100):
    Names.append(lines[i].split()[1].split('_')[0].split('clustercenter')[1])
    Seqs.append(lines[i].split()[0])

for l in range(len(Names)):
    name = Names[l]
    seq = Seqs[l]
    # Write the relaxed structures as torsion angle arrays into a text file
    Angles = []
    with open("After_SA_clustercenter"+name+"_0001_0001_15res.silent", 'r') as f:
        f.readline()
        f.readline()
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        if "SCORE:" in lines[i]:
            score = float(lines[i].split()[1])
            rmsd = float(lines[i].split()[23])
            i += 2
            temp_angle = []
            for j in range(15):
                i += 1
                temp_angle.append(float(lines[i].split()[2]))
                temp_angle.append(float(lines[i].split()[3]))
                temp_angle.append(float(lines[i].split()[4]))
            temp_angle.append(score)
            temp_angle.append(rmsd)
            Angles.append(temp_angle)
        i += 1
    
    np.savetxt('After_SA_clustercenter'+name+'_0001_0001_Angles_15res.txt', Angles, delimiter='\t', fmt='%f')