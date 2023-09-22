# CyclicPeptide
This is the design pipeline developed in my PhD thesis _"Heuristic energy-based cyclic peptide design"_, with my advisor Dennis Shasha and our collaborator Vikram Mulligan. It is capable of designing cyclic peptides of mixed L- and D-amino acids, as well as non-canonical amino acids. The macrocyle sizes attemped are 7, 15, 20, and 24 amino acids, and stable designs that pass the Rosetta energy landscape $P_{Near}=\frac{\sum_i exp(-\frac{rmsd_i^2}{\lambda^2}) exp(-\frac{E_i}{k_B T})}{\sum_i exp(-\frac{E_i}{k_B T})}>0.9$ test are found.

To run the pipeline, the following steps are taken.
1. Determine simulated annealing backbone sampling parameters.
2. Run simulated annealing for different backbone initial configurations.
3. Cluster sampled good backbone candidates.
4. Run Rosetta's _FastRelax_ to relax the backbones, and select backbones with low energies for sequence design.
5. Run Rosetta's _FastDesign_ to add sidechains, and select designs with low energies for energy landscape stability test.
6. Generate energy landscapes for selected designs.
