# Cyclic Peptide Design Pipeline
This is the design pipeline developed in my PhD thesis _"Heuristic energy-based cyclic peptide design"_, with my advisor Dennis Shasha and our collaborator Vikram Mulligan. It is capable of designing cyclic peptides of mixed L- and D-amino acids, as well as non-canonical amino acids. The macrocyle sizes attemped are 7, 15, 20, and 24 amino acids, and stable designs that pass the Rosetta energy landscape $P_{Near}=\frac{\sum_i exp(-\frac{rmsd_i^2}{\lambda^2}) exp(-\frac{E_i}{k_B T})}{\sum_i exp(-\frac{E_i}{k_B T})}>0.9$ test are found.

To run the pipeline, the following steps are taken.
1. Determine simulated annealing backbone sampling parameters.
   - All other parameters can be determined using linear interpolation as described in the thesis.
   - For parameter ranges of 0.6-0.9 for $k_0$, 15-19 for $b$, 2-4 for $c_{rama}$, 12-18 for $c_{rep}$, 14-20 for $c_{cyc}$, 16-22 for $c_{hbond}$, and 4-10 for $c_{other}$ (for small macrocycles of 7-10 residues), divide the ranges evenly into 3-8 values.
   - Perform combinatorial design to select parameter settings (pivots can be used to select more settings, suggested pivots are $k_0$ and $b$). To install the combinatorial design program, please check this website https://cs.nyu.edu/~shasha/papers/comb.html.
3. Run simulated annealing for different backbone initial configurations.
4. Cluster sampled good backbone candidates.
5. Run Rosetta's _FastRelax_ to relax the backbones, and select backbones with low energies for sequence design.
6. Run Rosetta's _FastDesign_ to add sidechains, and select designs with low energies for energy landscape stability test.
7. Generate energy landscapes for selected designs.
