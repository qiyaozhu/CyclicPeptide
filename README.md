# CyclicPeptide
This is the design pipeline developed in my PhD thesis "Heuristic energy-based cyclic peptide design", with my advisor Dennis Shasha and our collaborator Vikram Mulligan. It is capable of designing cyclic peptides of mixed L- and D-amino acids, as well as non-canonical amino acids. The macrocyle sizes attemped are 7, 15, 20, and 24 amino acids, and stable designs that pass the 
\[$P_{Near}=\frac{\sum_{i=1}^N exp(-\frac{rmsd_i^2}{\lambda^2}) exp(-\frac{E_i}{k_B T})}{\sum_{i=1}^N exp(-\frac{E_i}{k_B T})}$\]
