# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 13:24:56 2023

@author: Cathe
"""

with open('Pnear_seq.txt') as f:
    lines = f.readlines()
    
candidates = []

for l in lines:
    if l != '\n':
        candidates.append(l.split()[1])
        
for cand in candidates:
    with open(cand+'_rosetta.flags', 'w') as f:
        f.write('-in:file:silent Pnear_test_'+cand+'.txt\n')
        f.write('-in:file:fullatom\n')
        f.write('-in:file:native /scratch/qz886/Clustercenters_7res/Pnear_Candidates/output_'+cand+'.pdb\n')
        f.write('-parser:protocol relax_script.xml\n')
        f.write('-parser:script_vars Nres=7\n')
        f.write('-out:file:scorefile score_'+cand+'.sc\n')
        f.write('-out:file:silent '+cand+'_out.silent\n')
