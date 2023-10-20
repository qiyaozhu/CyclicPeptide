# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 13:24:56 2023

@author: Cathe
"""

with open('/scratch/qz886/Clustercenters_20res/Pnear_list.txt') as f:
    lines = f.readlines()
    
candidates = []

for l in lines:
    if l != '\n':
        candidates.append(l.split()[0])
        
for cand in candidates:
    with open('Pnear_'+cand+'_Cartesian_rosetta.flags', 'w') as f:
        f.write('-nstruct 100000\n')
        f.write('-MPI_batchsize_by_level 20\n')
        f.write('-MPI_processes_by_level 1 95\n')
        f.write('-in:file:native output_'+cand+'.pdb\n')
        f.write('-out:file:silent Pnear_'+cand+'_Cartesian_out.silent\n')
        f.write('-cyclic_peptide:sequence_file Pnear_'+cand+'_seq.txt\n')
        f.write('-cyclic_peptide:MPI_output_fraction 1\n')
        f.write('-cyclic_peptide:MPI_choose_highest false\n')
        f.write('-cyclic_peptide:MPI_sort_by energy\n')
        f.write('-score:symmetric_gly_tables true\n')
        f.write('-cyclic_peptide:genkic_closure_attempts 100\n')
        f.write('-cyclic_peptide:genkic_min_solution_count 1\n')
        f.write('-cyclic_peptide:use_rama_filter true\n')
        f.write('-cyclic_peptide:rama_cutoff 3.0\n')
        f.write('-cyclic_peptide:min_genkic_hbonds 0\n')
        f.write('-cyclic_peptide:min_final_hbonds 0\n')
        f.write('-cyclic_peptide:cartesian_relax_rounds 3\n')
        f.write('-mute all\n')
        f.write('-unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary\n')
        
        
        
