# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 21:38:43 2023

@author: Cathe
"""

with open('all.RamaProb', 'r') as f:
    lines = f.readlines()
    
i = 0
res = ''

while i < len(lines):
    if lines[i][0] != '#':
        name = lines[i].split()[0]
        if name != res:
            res = name
            with open(res+'.dat', 'w') as f:
                while i < len(lines):
                    if lines[i][0] != '#':
                        f.write(lines[i])
                    else:
                        break
                    i += 1
    else:
        i += 1
                