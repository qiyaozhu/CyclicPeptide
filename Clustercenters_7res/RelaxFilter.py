# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 23:52:15 2023

@author: Cathe
"""

with open('FastRelax_silent.txt', 'r') as f:
    lines = f.readlines()
    
Candidates = [lines[0],lines[1]]

i = 2
while i < len(lines):
    if float(lines[i].split()[1]) < 8:
        for j in range(10):
            Candidates.append(lines[i+j])
    i += 10

with open('FastRelax_filtered.txt', 'w') as f:
    for l in Candidates:
        f.write(l)
        
    
