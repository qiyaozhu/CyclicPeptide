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

with open('Pnear_list.txt') as f:
    lines = f.readlines()
    
candidates = []
scores = []
seqs = []

for l in lines:
    if l != "\n":
        candidates.append(l.split()[0])
        scores.append(float(l.split()[1]))
        
        seq = l.split()[2]
        res1 = seq.split('[')[1].split(':')[0]
        resn = seq.split('[')[-1].split(':')[0]
        if res1 == 'HIS_D':
            res1 = 'DHIS'
        if resn == 'HIS_D':
            resn = 'DHIS'
        residues = [res1]

        lpos = [m.start() for m in re.finditer('\[', seq)]
        rpos = [m.start() for m in re.finditer('\]', seq)]
        seq = seq[rpos[0]+1:lpos[-1]-1]
        
        p = 0
        while p < len(seq):
            if p == len(seq)-1:
                residues.append(aacode[seq[p]])
                p += 1
            else:
                if seq[p+1] != '[':
                    residues.append(aacode[seq[p]])
                    p += 1
                else:
                    if seq[p+2:p+6] == 'HIS_':
                        residues.append('DHIS')
                        p += 8
                    else:
                        residues.append(seq[p+2:p+6])
                        p += 7
        residues.append(resn)
        
        s = residues[0]
        for r in range(1,len(residues)):
            s = s+'-'+residues[r]        
        seqs.append(s)


lamda = 1.5
kbT = 0.62
Pnear = []

with open('Pnear_Cartesian_15res_1.out', 'r') as f:
    plines = f.readlines()
with open('Pnear_Cartesian_15res_2.out', 'r') as f:
    plines2 = f.readlines()
with open('Pnear_Cartesian_15res_3.out', 'r') as f:
    plines3 = f.readlines()
plines += plines2
plines += plines3
    
l = 0
while l < len(plines):
    if len(plines[l].split()) > 0:
        if plines[l].split()[0] == "MPI_worker_node":
            l += 1            
            Energy = []
            RMSD = []
            while plines[l].split()[0] != "End":
                Energy.append(float(plines[l].split()[3]))
                RMSD.append(float(plines[l].split()[2]))
                l += 1
            p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
            Pnear.append(p)
        else:
            l += 1
    else:
        l += 1


PnearM = []
for i in range(len(candidates)):
    path = 'After_GA_'+candidates[i]+'_15res.mat'
    if os.path.exists(path):
        mat = scipy.io.loadmat(path)
        Energy = mat["candScores"][0]
        RMSD = mat["candRMSD"][0]
        p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
        PnearM.append(p)

# PnearHigh = []
# PnearMHigh = []
# PnearDiff = []
# with open('Pnear_values_15res.txt', 'w') as f:
#     f.write('Name\t\t\t\tEnergy\t\tRosetta\t\tGA\t\tSequence\n')
#     for i in range(len(PnearM)):
#         if Pnear[i] > 0.9:
#             PnearHigh.append(i)
#             f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t\t{PnearM[i]:.3f}\t\t'+seqs[i]+'\n')
#     f.write('\nName\t\t\t\tEnergy\t\tRosetta\t\tGA\t\tSequence\n')
#     for i in range(len(PnearM)):
#         if PnearM[i] > 0.9:
#             PnearMHigh.append(i)
#             f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t\t{PnearM[i]:.3f}\t\t'+seqs[i]+'\n')
#     f.write('\nName\t\t\t\tEnergy\t\tRosetta\t\tGA\t\tSequence\n')
#     for i in range(len(PnearM)):
#         if abs(Pnear[i]-PnearM[i]) > 0.4:
#             PnearDiff.append(i)
#             f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t\t{PnearM[i]:.3f}\t\t'+seqs[i]+'\n')

mat = scipy.io.loadmat('Rg.mat')
Rg = mat["Rg"][0]

Hbonds = []
for i in range(len(PnearM)):
    with open('output_'+candidates[i]+'.pdb') as f:
        pdblines = f.readlines()
    for l in pdblines:
        if len(l.split()) > 1:
            if l.split()[0] == 'min_internal_hbonds':
                Hbonds.append(l.split()[1])
                
with open('Pnear_values_15res.txt', 'w') as f:
    f.write('Name\t\t\t\tEnergy\t\tRosetta\tGA\tHbond\tRadius\tSequence\n')
    for i in range(len(PnearM)):
        f.write(candidates[i].split('_')[0]+f'\t\t{scores[i]:.3f}\t\t{Pnear[i]:.3f}\t{PnearM[i]:.3f}\t'+Hbonds[i]+f'\t{Rg[i]:.3f}\t'+seqs[i]+'\n')

plt.rc('font', **font)
plt.figure(figsize=(10,10))
plt.plot(Rg, PnearM, 'x', markersize=15, markeredgewidth=3)
plt.axhline(y=0.9, color='red', linestyle='--', linewidth=3)
plt.title('Structure Compactness')
plt.xlabel(r'Backbone Radius ($\AA$)')
plt.ylabel('Pnear GA')
plt.savefig('Pnear_Rg_15res.pdf')


# # Pnear correlation plot
# Pnear = Pnear[:len(PnearM)]
# scores = scores[:len(PnearM)]
# plt.rc('font', **font)
# plt.figure(figsize=(10,10))
# plt.plot(Pnear, PnearM, 'x', markersize=15, markeredgewidth=3)
# plt.axhline(y=0.9, color='red', linestyle='--', linewidth=3)
# plt.axvline(x=0.9, color='red', linestyle='-.', linewidth=3)
# plt.title('Pnear Correlation')
# plt.xlabel('Rosetta')
# plt.ylabel('Genetic Algorithm')
# plt.savefig('Pnear_correlation_15res.pdf')



# # Pnear distribution plot
# Pnear = np.array(Pnear)
# scores = np.array(scores)

# # Define score and Pnear bins
# score_bins = np.array([-43, -37, -36, -35, -34])
# Pnear_bins = np.array([0, 0.7, 0.9, 1])

# # Compute histogram frequencies
# freq = np.zeros((len(score_bins)-1, len(Pnear_bins)-1))
# prob = np.zeros((len(score_bins)-1, len(Pnear_bins)-1))
# for i in range(len(score_bins) - 1):
#     for j in range(len(Pnear_bins) - 1):
#         mask = np.logical_and(scores >= score_bins[i], scores < score_bins[i+1])
#         mask = np.logical_and(mask, Pnear >= Pnear_bins[j])
#         mask = np.logical_and(mask, Pnear < Pnear_bins[j+1])
#         freq[i,j] = np.sum(mask)
#     prob[i,:] = freq[i,:] / np.sum(freq[i,:])

# # Create stacked bar chart
# fig, ax = plt.subplots(figsize=(15,10))
# ax.bar(score_bins[:-1]+0.25, prob[:,0], label='0<Pnear<0.7', alpha=1, width=0.5)
# ax.bar(score_bins[:-1]+0.75, prob[:,1], label='0.7<Pnear<0.9', alpha=1, width=0.5)
# ax.bar(score_bins[:-1]+1.25, prob[:,2], label='0.9<Pnear', alpha=1, width=0.5)

# for i in range(len(score_bins) - 1):
#     for j in range(len(Pnear_bins) - 1):
#         ax.text(score_bins[i]+0.5*j+0.25, prob[i,j], str(int(freq[i,j])), ha='center', va='bottom')

# # Set plot title and axis labels
# ax.set_title("Rosetta Pnear distribution")
# ax.set_xlabel("Energy bins (kcal/mol)")
# ax.set_ylabel("Relative Frequency")
# ax.set_ylim(0,1)
# ax.set_xticks([-43, -37, -36, -35, -34])
# ax.legend()
# plt.savefig('Pnear_distribution_Rosetta.pdf')
# plt.show()


# # Pnear distribution plot
# PnearM = np.array(PnearM)

# # Define score and Pnear bins
# score_bins = np.array([-43, -37, -36, -35, -34])
# Pnear_bins = np.array([0, 0.7, 0.9, 1])

# # Compute histogram frequencies
# freq = np.zeros((len(score_bins)-1, len(Pnear_bins)-1))
# prob = np.zeros((len(score_bins)-1, len(Pnear_bins)-1))
# for i in range(len(score_bins) - 1):
#     for j in range(len(Pnear_bins) - 1):
#         mask = np.logical_and(scores >= score_bins[i], scores < score_bins[i+1])
#         mask = np.logical_and(mask, PnearM >= Pnear_bins[j])
#         mask = np.logical_and(mask, PnearM < Pnear_bins[j+1])
#         freq[i,j] = np.sum(mask)
#     prob[i,:] = freq[i,:] / np.sum(freq[i,:])

# # Create stacked bar chart
# fig, ax = plt.subplots(figsize=(15,10))
# ax.bar(score_bins[:-1]+0.25, prob[:,0], label='0<Pnear<0.7', alpha=1, width=0.5)
# ax.bar(score_bins[:-1]+0.75, prob[:,1], label='0.7<Pnear<0.9', alpha=1, width=0.5)
# ax.bar(score_bins[:-1]+1.25, prob[:,2], label='0.9<Pnear', alpha=1, width=0.5)

# for i in range(len(score_bins) - 1):
#     for j in range(len(Pnear_bins) - 1):
#         ax.text(score_bins[i]+0.5*j+0.25, prob[i,j], str(int(freq[i,j])), ha='center', va='bottom')

# # Set plot title and axis labels
# ax.set_title("GA Pnear distribution")
# ax.set_xlabel("Energy bins (kcal/mol)")
# ax.set_ylabel("Relative Frequency")
# ax.set_ylim(0,1)
# ax.set_xticks([-43, -37, -36, -35, -34])
# plt.savefig('Pnear_distribution_GA.pdf')
# plt.show()



# ########## Energy landscape comparison plots ##########
# # landscapes = {2:"169032", 8:"3114", 13:"116427"}

# lamda = 1.5
# kbT = 0.62

# count = 0
# l = 0

# plt.rc('font', **font)
# minE = -55
# xmax = 6
# maxE = 0

# while l < len(plines):
#     if len(plines[l].split()) > 0:
#         if plines[l].split()[0] == "MPI_worker_node":
#             if count in PnearHigh:
#                 l += 1
                
#                 Energy = []
#                 RMSD = []
#                 while plines[l].split()[0] != "End":
#                     Energy.append(float(plines[l].split()[3]))
#                     RMSD.append(float(plines[l].split()[2]))
#                     l += 1
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))               
                
#                 fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))
#                 ax1.plot(RMSD, Energy, '.', markersize=12)
#                 ax1.set_title(f'Rosetta Pnear={p:.3f}')
#                 ax1.set_xlabel(r'RMSD to design ($\AA$)')
#                 ax1.set_ylabel('Energy (kcal/mol)')
#                 ax1.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax1.transAxes, ha='right', va='bottom')
#                 ax1.set_xlim([0, xmax])
#                 ax1.set_ylim([minE, maxE])
                
#                 Energy = []
#                 RMSD = []
#                 path = 'After_GA_'+candidates[count]+'_15res.mat'
#                 mat = scipy.io.loadmat(path)
#                 Energy = mat["candScores"][0]
#                 RMSD = mat["candRMSD"][0]
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
                
#                 ax2.plot(RMSD, Energy, '.', markersize=12)
#                 ax2.set_title(f'GA Pnear={p:.3f}')
#                 ax2.set_xlabel(r'RMSD to design ($\AA$)')
#                 # ax2.set_ylabel('Energy (kcal/mol)')
#                 ax2.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax2.transAxes, ha='right', va='bottom')
#                 ax2.set_xlim([0, xmax])
#                 ax2.set_ylim([minE, maxE])
#                 plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
#                 plt.savefig('RosettaHigh_'+candidates[count].split('_')[0].split('clustercenter')[1]+'.png')
            
#             elif count in PnearMHigh:
#                 l += 1
                
#                 Energy = []
#                 RMSD = []
#                 while plines[l].split()[0] != "End":
#                     Energy.append(float(plines[l].split()[3]))
#                     RMSD.append(float(plines[l].split()[2]))
#                     l += 1
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))               
                
#                 fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))
#                 ax1.plot(RMSD, Energy, '.', markersize=12)
#                 ax1.set_title(f'Rosetta Pnear={p:.3f}')
#                 ax1.set_xlabel(r'RMSD to design ($\AA$)')
#                 ax1.set_ylabel('Energy (kcal/mol)')
#                 ax1.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax1.transAxes, ha='right', va='bottom')
#                 ax1.set_xlim([0, xmax])
#                 ax1.set_ylim([minE, maxE])
                
#                 Energy = []
#                 RMSD = []
#                 path = 'After_GA_'+candidates[count]+'_15res.mat'
#                 mat = scipy.io.loadmat(path)
#                 Energy = mat["candScores"][0]
#                 RMSD = mat["candRMSD"][0]
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
                
#                 ax2.plot(RMSD, Energy, '.', markersize=12)
#                 ax2.set_title(f'GA Pnear={p:.3f}')
#                 ax2.set_xlabel(r'RMSD to design ($\AA$)')
#                 # ax2.set_ylabel('Energy (kcal/mol)')
#                 ax2.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax2.transAxes, ha='right', va='bottom')
#                 ax2.set_xlim([0, xmax])
#                 ax2.set_ylim([minE, maxE])
#                 plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
#                 plt.savefig('GAHigh_'+candidates[count].split('_')[0].split('clustercenter')[1]+'.png')
            
#             elif count in PnearDiff:
#                 l += 1
                
#                 Energy = []
#                 RMSD = []
#                 while plines[l].split()[0] != "End":
#                     Energy.append(float(plines[l].split()[3]))
#                     RMSD.append(float(plines[l].split()[2]))
#                     l += 1
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))               
                
#                 fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))
#                 ax1.plot(RMSD, Energy, '.', markersize=12)
#                 ax1.set_title(f'Rosetta Pnear={p:.3f}')
#                 ax1.set_xlabel(r'RMSD to design ($\AA$)')
#                 ax1.set_ylabel('Energy (kcal/mol)')
#                 ax1.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax1.transAxes, ha='right', va='bottom')
#                 ax1.set_xlim([0, xmax])
#                 ax1.set_ylim([minE, maxE])
                
#                 Energy = []
#                 RMSD = []
#                 path = 'After_GA_'+candidates[count]+'_15res.mat'
#                 mat = scipy.io.loadmat(path)
#                 Energy = mat["candScores"][0]
#                 RMSD = mat["candRMSD"][0]
#                 p = sum(np.exp(-np.array(RMSD)**2 / lamda**2) * np.exp(-np.array(Energy) / kbT)) / sum(np.exp(-np.array(Energy) / kbT))
                
#                 ax2.plot(RMSD, Energy, '.', markersize=12)
#                 ax2.set_title(f'GA Pnear={p:.3f}')
#                 ax2.set_xlabel(r'RMSD to design ($\AA$)')
#                 # ax2.set_ylabel('Energy (kcal/mol)')
#                 ax2.text(0.98, 0.02, f"{len(RMSD):d} samples", transform=ax2.transAxes, ha='right', va='bottom')
#                 ax2.set_xlim([0, xmax])
#                 ax2.set_ylim([minE, maxE])
#                 plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
#                 plt.savefig('Difference_'+candidates[count].split('_')[0].split('clustercenter')[1]+'.png')
#             else:
#                 l += 1
#             count += 1
#         else:
#             l += 1
#     else:
#         l += 1