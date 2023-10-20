# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 23:52:15 2023

@author: Cathe
"""

with open('FastRelax_silent.txt', 'r') as f:
    lines = f.readlines()
    
Candidates = [lines[0],lines[1]]

n = 24

i = 2
while i < len(lines):
    if float(lines[i].split()[1]) < 0:
        for j in range(n+3):
            Candidates.append(lines[i+j])
    i = i+n+3

with open('FastRelax_filtered.txt', 'w') as f:
    for l in Candidates:
        f.write(l)
        
    
