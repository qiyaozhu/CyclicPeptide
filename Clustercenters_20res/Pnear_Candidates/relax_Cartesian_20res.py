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
with open("Pnear_seq.txt", 'r') as f:
    lines = f.readlines()
for l in lines:
    Names.append(l.split()[1])
    
for name in Names:
    with open(name+'_relax_Cartesian_rosetta_20res.flags', 'w') as f:
        f.write('-in:file:s output_'+name+'.pdb\n'+
                '-in:file:fullatom\n'+
                '-parser:protocol relax_Cartesian.xml\n'+
                '-parser:script_vars Nres=20\n'+
                '-out:file:scorefile output_'+name+'.sc\n')

with open('run_Relax_Cartesian_20res.sbatch', 'w') as f:
    f.write('#!/bin/bash\n'+
            '#SBATCH --nodes=1\n'+
            '#SBATCH --tasks-per-node=8\n'+
            '#SBATCH --mem=10GB\n'+
            '#SBATCH --time=72:00:00\n'+
            '#SBATCH --job-name=PnearRelax\n'+
            '#SBATCH --mail-type=ALL\n'+
            '#SBATCH --mail-user=qz886@nyu.edu\n\n'+
            'module purge\n'+
            'module rosetta/openmpi/intel/2020.46.61480\n\n')
    for name in Names:
        f.write('srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @'+name+'_relax_Cartesian_rosetta_20res.flags\n')

