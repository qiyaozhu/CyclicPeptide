# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 13:24:56 2023

@author: Cathe
"""

with open('/scratch/qz886/design_candidates.txt') as f:
    lines = f.readlines()
    
candidates = []

for l in lines:
    if l != '\n':
        candidates.append(l.split('.')[0])
        
for cand in candidates:
    with open(cand+'_rosetta.flags', 'w') as f:
        f.write('-nstruct 10000\n')
        f.write('-MPI_batchsize_by_level 20\n')
        f.write('-MPI_processes_by_level 1 47\n')
        f.write('-in:file:native '+cand+'.pdb\n')
        f.write('-out:file:silent '+cand+'_out.silent\n')
        f.write('-cyclic_peptide:sequence_file '+cand+'_seq.txt\n')
        f.write('-cyclic_peptide:MPI_output_fraction 1\n')
        f.write('-cyclic_peptide:MPI_choose_highest false\n')
        f.write('-cyclic_peptide:MPI_sort_by energy\n')
        f.write('-score:symmetric_gly_tables true\n')
        f.write('-cyclic_peptide:genkic_closure_attempts 250\n')
        f.write('-cyclic_peptide:genkic_min_solution_count 1\n')
        f.write('-cyclic_peptide:use_rama_filter true\n')
        f.write('-cyclic_peptide:rama_cutoff 3.0\n')
        f.write('-cyclic_peptide:min_genkic_hbonds 2\n')
        f.write('-cyclic_peptide:min_final_hbonds 2\n')
        f.write('-mute all\n')
        f.write('-unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary\n')
        
        
        