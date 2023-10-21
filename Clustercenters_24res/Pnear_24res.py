# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:00:01 2023

@author: Cathe
"""

import matplotlib.pyplot as plt
import os
import numpy as np
import scipy.io
import re
from matplotlib import rc
from scipy.stats import pearsonr

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
    

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 30}

aacode = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

# with open('Pnear_list.txt') as f:
#     lines = f.readlines()
    
# candidates = []
# scores = []
# seqs = []

# for l in lines:
#     if l != "\n":
#         candidates.append(l.split()[0])
#         scores.append(float(l.split()[1]))
        
#         seq = l.split()[2]
#         res1 = seq.split('[')[1].split(':')[0]
#         resn = seq.split('[')[-1].split(':')[0]
#         if res1 == 'HIS_D':
#             res1 = 'DHIS'
#         if resn == 'HIS_D':
#             resn = 'DHIS'
#         residues = [res1]

#         lpos = [m.start() for m in re.finditer('\[', seq)]
#         rpos = [m.start() for m in re.finditer('\]', seq)]
#         seq = seq[rpos[0]+1:lpos[-1]-1]
        
#         p = 0
#         while p < len(seq):
#             if p == len(seq)-1:
#                 residues.append(aacode[seq[p]])
#                 p += 1
#             else:
#                 if seq[p+1] != '[':
#                     residues.append(aacode[seq[p]])
#                     p += 1
#                 else:
#                     if seq[p+2:p+6] == 'HIS_':
#                         residues.append('DHIS')
#                         p += 8
#                     else:
#                         residues.append(seq[p+2:p+6])
#                         p += 7
#         residues.append(resn)
        
#         s = residues[0]
#         for r in range(1,len(residues)):
#             s = s+'-'+residues[r]        
#         seqs.append(s)


# lamda = 2
# kbT = 0.62
# Pnear = []

# with open('Pnear_Cartesian_24res.out', 'r') as f:
#     plines = f.readlines()
    
# l = 0
# while l < len(plines):
#     if len(plines[l].split()) > 0:
#         if plines[l].split()[0] == "MPI_worker_node":
#             l += 1            
#             Energy = []
#             RMSD = []
#             while plines[l].split()[0] != "End":
#                 Energy.append(float(plines[l].split()[3]))
#                 RMSD.append(float(plines[l].split()[2]))
#                 l += 1
#             p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
#             Pnear.append(p)
#         else:
#             l += 1
#     else:
#         l += 1
        

# PnearM = []
# IndexM = []
# for i in range(len(candidates)):
#     path = 'After_GA_'+candidates[i]+'_24res.mat'
#     if os.path.exists(path):
#         mat = scipy.io.loadmat(path)
#         Energy = mat["candScores"][0]
#         RMSD = mat["candRMSD"][0]
#         p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
#         PnearM.append(p)
#         IndexM.append(i)

# mat = scipy.io.loadmat('Rg.mat')
# Rg = mat["Rg"][0]

# mat = scipy.io.loadmat('ReshapePnear.mat')
# Reshape = mat["Reshape"][0]

# Hbonds = []
# for i in range(len(candidates)):
#     with open('output_'+candidates[i]+'.pdb') as f:
#         pdblines = f.readlines()
#     for l in pdblines:
#         if len(l.split()) > 1:
#             if l.split()[0] == 'min_internal_hbonds':
#                 Hbonds.append(l.split()[1])

# with open('Pnear_values_24res.txt', 'w') as f:
#     f.write('Name\t\t\t\tEnergy\t\tRosetta\tGA\tReshape\tHbond\tRadius\tSequence\n')
#     for i in range(len(Pnear)):
#         f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t{PnearM[i]:.3f}\t{Reshape[i]:.3f}\t'+Hbonds[i]+f'\t{Rg[i]:.3f}\t'+seqs[i]+'\n')
#     for i in range(len(Pnear), len(PnearM)):
#         f.write(candidates[IndexM[i]].split('_')[0]+f'\t\t{scores[IndexM[i]]:.3f}\t\t\t{PnearM[i]:.3f}\t{Reshape[i]:.3f}\t'+Hbonds[IndexM[i]]+f'\t{Rg[IndexM[i]]:.3f}\t'+seqs[IndexM[i]]+'\n')


# plt.rc('font', **font)
# plt.figure(figsize=(10,10))
# plt.plot(Rg[:len(Pnear)], PnearM[:len(Pnear)], 'x', markersize=15, markeredgewidth=3)
# plt.axhline(y=0.9, color='red', linestyle='--', linewidth=3)
# plt.title('Structure Compactness')
# plt.xlabel(r'Backbone Radius ($\AA$)')
# plt.ylabel('Pnear GA')
# plt.savefig('Pnear_Rg_24res.pdf')


# # Pnear correlation plot
# plt.rc('font', **font)
# plt.figure(figsize=(10,10))
# plt.plot(Pnear, PnearM[:len(Pnear)], 'x', markersize=15, markeredgewidth=3)
# plt.axhline(y=0.9, color='red', linestyle='--', linewidth=3)
# plt.axvline(x=0.9, color='red', linestyle='-.', linewidth=3)
# plt.title('Pnear Correlation')
# plt.xlabel('Rosetta')
# plt.ylabel('Genetic Algorithm')
# plt.savefig('Pnear_correlation_24res.pdf')




########## Energy landscape comparison plots ##########

# lamda = 2
# kbT = 0.62

# count = 0
# l = 0

# plt.rc('font', **font)
# minE = -90
# xmax = 9
# maxE = 0

# while l < len(plines):
#     if len(plines[l].split()) > 0:
#         if plines[l].split()[0] == "MPI_worker_node":
#             l += 1
            
#             Energy = []
#             RMSD = []
#             while plines[l].split()[0] != "End":
#                 Energy.append(float(plines[l].split()[3]))
#                 RMSD.append(float(plines[l].split()[2]))
#                 l += 1
#             p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))               
            
#             fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(28,10))
#             ax1.plot(RMSD, Energy, '.', markersize=12)
#             ax1.set_title(f'Rosetta Pnear={p:.3f}')
#             ax1.set_xlabel(r'RMSD to design ($\AA$)')
#             ax1.set_ylabel('Energy (kcal/mol)')
#             ax1.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax1.transAxes, ha='right', va='bottom')
#             ax1.set_xlim([0, xmax])
#             ax1.set_ylim([minE, maxE])
            
#             Energy = []
#             RMSD = []
#             path = 'After_GA_'+candidates[count]+'_24res.mat'
#             if os.path.exists(path):
#                 mat = scipy.io.loadmat(path)
#                 Energy = mat["candScores"][0]
#                 RMSD = mat["candRMSD"][0]
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
                
#                 ax2.plot(RMSD, Energy, '.', markersize=12)
#                 ax2.set_title(f'GA Pnear={p:.3f}')
#                 ax2.set_xlabel(r'RMSD to design ($\AA$)')
#                 ax2.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax2.transAxes, ha='right', va='bottom')
#                 ax2.set_xlim([0, xmax])
#                 ax2.set_ylim([minE, maxE])
#                 plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
                
#                 Energy = []
#                 RMSD = []
#                 path = 'After_GA_'+candidates[count]+'_24res_reshape.mat'
#                 mat = scipy.io.loadmat(path)
#                 Energy = mat["Score"][0]
#                 RMSD = mat["RMSD"][0]
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
                
#                 ax3.plot(RMSD, Energy, '.', markersize=12)
#                 ax3.set_title(f'GA Reshape Pnear={p:.3f}')
#                 ax3.set_xlabel(r'RMSD to design ($\AA$)')
#                 ax3.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax2.transAxes, ha='right', va='bottom')
#                 ax3.set_xlim([0, xmax])
#                 ax3.set_ylim([minE, maxE])
#                 plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
#                 plt.savefig('Landscape_'+candidates[count].split('_')[0].split('clustercenter')[1]+'.png')


#             count += 1
#         else:
#             l += 1
#     else:
#         l += 1
        


########## Energy landscape reshape ##########

# lamda = 2
# kbT = 0.62

# count = 0
# l = 0

# plt.rc('font', **font)
# minE = -90
# xmax = 9
# maxE = 0

# for count in range(len(candidates)):
#         Energy = []
#         RMSD = []
#         path = 'After_GA_'+candidates[count]+'_24res.mat'
#         if os.path.exists(path):
#             mat = scipy.io.loadmat(path)
#             Energy = mat["candScores"][0]
#             RMSD = mat["candRMSD"][0]
#             p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
        
#             fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,10))
#             ax1.plot(RMSD, Energy, '.', markersize=12)
#             ax1.set_title(f'GA Pnear={p:.3f}')
#             ax1.set_xlabel(r'RMSD to design ($\AA$)')
#             ax1.set_ylabel('Energy (kcal/mol)')
#             ax1.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax1.transAxes, ha='right', va='bottom')
#             ax1.set_xlim([0, xmax])
#             ax1.set_ylim([minE, maxE])
        
#             Energy = []
#             RMSD = []
#             path = 'After_GA_'+candidates[count]+'_24res_reshape.mat'
#             mat = scipy.io.loadmat(path)
#             Energy = mat["Score"][0]
#             RMSD = mat["RMSD"][0]
#             p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
            
#             ax2.plot(RMSD, Energy, '.', markersize=12)
#             ax2.set_title(f'GA Reshape Pnear={p:.3f}')
#             ax2.set_xlabel(r'RMSD to design ($\AA$)')
#             ax2.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax2.transAxes, ha='right', va='bottom')
#             ax2.set_xlim([0, xmax])
#             ax2.set_ylim([minE, maxE])
#             plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
#             plt.savefig('Reshape_'+candidates[count].split('_')[0].split('clustercenter')[1]+'.png')



###################### S4 Symmetry Validation #############################
lamda = 2
kbT = 0.62

plt.rc('font', **font)
minE = -100
xmax = 9
maxE = 0

path = 'After_GA_clustercenterS4-Design_0001_0001_24res.mat'
mat = scipy.io.loadmat(path)
Energy = mat["candScores"][0]
RMSD = mat["candRMSD"][0]
p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))

plt.rc('font', **font)
plt.figure(figsize=(10,10))
plt.plot(RMSD, Energy, '.', markersize=12)
plt.title(f'GA Pnear={p:.3f}')
plt.xlabel(r'RMSD to design ($\AA$)')
plt.ylabel('Energy (kcal/mol)')
plt.text(0.98, 0.02, f"{len(RMSD):d} samples", ha='right', va='bottom', transform=plt.gca().transAxes)
plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.15)
plt.xlim([0, xmax])
plt.ylim([minE, maxE])
plt.savefig('Landscape_S4.png')