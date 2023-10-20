# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:28:46 2023

@author: Cathe
"""

import scipy.io
import sys

def writeSilent(name, seq, aacode, score, phi, psi, omega, x, y, z, output):
    olseq = ''
    for aa in seq:
        olseq += aacode[aa]
        
    annoseq = 'ANNOTATED_SEQUENCE: '+olseq[0]+'['+seq[0]+':protein_cutpoint_upper]'
    for i in range(1,len(seq)-1):
        if seq[i][0] == 'D':
            annoseq = annoseq+olseq[i]+'['+seq[i]+']'
        else:
            annoseq += olseq[i]
    annoseq = annoseq+olseq[-1]+'['+seq[-1]+':protein_cutpoint_lower] '+name+'\n'
    
    scoreline1 = 'SCORE:    score    description\n'
    scoreline2 = f'SCORE:    {score:.3f}    '+name+'\n'
    
    structure = ''
    for i in range(len(phi)):
        newres = f'{i+1:4d}  L{phi[i]:10.3f}{psi[i]:9.3f}{omega[i]:9.3f}{x[i]:9.3f}{y[i]:9.3f}{z[i]:9.3f}    0.000    0.000    0.000    0.000 '+name+'\n'
        structure += newres
    
    with open(output, 'a') as f:
        f.write('SEQUENCE: '+olseq+'\n')
        f.write(scoreline1)
        f.write(scoreline2)
        f.write('REMARK PROTEIN SILENTFILE\n')
        f.write(annoseq)
        f.write(structure)


aacode = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'DALA':'A', 'DARG':'R', 'DASN':'N', 'DASP':'D', 'DCYS':'C', 'DGLN':'Q', 'DGLU':'E', 'DGLY':'G', 'DHIS':'H', 'DILE':'I', 'DLEU':'L', 'DLYS':'K', 'DMET':'M', 'DPHE':'F', 'DPRO':'P', 'DSER':'S', 'DTHR':'T', 'DTRP':'W', 'DTYR':'Y', 'DVAL':'V'}

with open("Pnear_seq.txt", 'r') as f:
    lines = f.readlines()

for cand in range(12,100):
    name = lines[cand].split()[1]
    path = 'After_SA_'+name+'_15res.mat'
    mat = scipy.io.loadmat(path)
    candPhi = mat["candPhi"]
    candPsi = mat["candPsi"]
    candCoordinates_x = mat["candCoordinates_x"]
    candCoordinates_y = mat["candCoordinates_y"]
    candCoordinates_z = mat["candCoordinates_z"]
    candOmega = mat["candOmega"]
    candScores = mat["candScores"][0]

    seq = lines[cand].split()[0].split('-')
    seq = [s.upper() for s in seq]
    output = 'After_SA_'+name+'_15res.txt'

    for i in range(len(candScores)):
        name = "Pnear_sample"+str(i+1)
        score = candScores[i]
        phi = candPhi[:,i]
        psi = candPsi[:,i]
        omega = candOmega[:,i]
        x = candCoordinates_x[:,i]
        y = candCoordinates_y[:,i]
        z = candCoordinates_z[:,i]
        writeSilent(name, seq, aacode, score, phi, psi, omega, x, y, z, output)
    

