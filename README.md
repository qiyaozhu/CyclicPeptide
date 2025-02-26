# Cyclic Peptide Design Pipeline
This is the cyclic peptide design pipeline, _CyclicChamp_, developed in the paper _"Heuristic energy-based cyclic peptide design"_, by Qiyao Zhu, Vikram Mulligan, and Dennis Shasha. It is capable of designing cyclic peptides of mixed L- and D-amino acids, as well as non-canonical amino acids. The macrocyle sizes attemped in the paper were 7, 15, 20, and 24 amino acids, and stable designs that pass the Rosetta energy landscape $P_{Near}=\frac{\sum_i exp(-\frac{rmsd_i^2}{\lambda^2}) exp(-\frac{E_i}{k_B T})}{\sum_i exp(-\frac{E_i}{k_B T})}>0.9$ test were found and posted here.

This pipelines requires installation of Rosetta on your device. For installation instructions, see https://docs.rosettacommons.org/docs/latest/getting_started/Getting-Started

To run the pipeline for _n_ residue macrocycle design, the following steps are taken.
1. Determine simulated annealing backbone sampling parameters for Ramachandran energy ($rama$), repulsive energy ($rep$), cyclic backbone closure deviation ($cyc$), hydrogen bond energy ($hbond$), and number of hydrogen bonds ($count$).
   - Energy thresholds: $E_{thr,rama}=8n$, $E_{thr,rep}=10+\frac{(n-7)*10}{17}$, $E_{thr,cyc}=1$, $H_{thr,count}=\lceil n/3\rceil$
   - Good backbone candidate criteria: $E_{cri,rep}=5+\frac{(n-7)*10}{17}$, $E_{cri,cyc}=1$, $H_{cri,count}=\lceil n/3\rceil$
   - Initial temperatures for simulated annealing: $T_{0,rama}=10+\frac{(n-7)*20}{17}$, $T_{0,rep}=20+\frac{(n-7)*80}{17}$, $T_{0,cyc}=2+\frac{(n-7)*4}{17}$, $T_{0,hbond}=2+\frac{(n-7)*4}{17}$
   - Random move disk radius: $k_0 \in [0.5,1]$ and smaller value for larger $n$, $b \in [15,18]$ and larger value for larger $n$
   - Temperature dropping rates: $c_{rama} \sim 4$, $c_{rep} \sim 14$, $c_{cyc} \sim 18$, $c_{hbond} \sim 20$ [*]
2. Run simulated annealing Energy_*n*residue.mat to find good backbone candidates. Two 15-residue sample runs are displayed in SampleRun_15res_4.pse, SampleRun_15res_10.pse.
3. Run Clustering.m in Clustercenters_*n*res folders to cluster sampled good backbone candidates.
4. Run Rosetta's _FastRelax_ to relax the backbones, and select backbones with low energies for sequence design. Example scripts are provided as relax_script.xml.
5. Run Rosetta's _FastDesign_ to add sidechains, and select designs with low energies for energy landscape stability test. Example scripts are provided as design_script.xml.
6. Generate energy landscapes for selected designs.
   - For 7 residues, use _Ramachandran-stability filtering_ test in Clustercenters_7res/Pnear_Modified/Pnear_sampling.m
   - For 15-24 residues, use _ClusterGen_ test. Run Clustercenters_*n*res/Pnear_Modified/Pnear_sampling_SA_*n*res.m. Perform Rosetta's _FastRelax_ with designed sequences to form initial populations for genetic algorithm. Then run Clustercenters_*n*res/Pnear_Modified/Pnear_sampling_GA_*n*res.m.
  
The top designs ($P_{Near}>0.9$) generated by _CyclicChamp_ are provided in GoodDesigns.

Other useful resources:
1. Our energy functions are implemented in my_functions folder
2. The Ramachandran spaces used are stored in RamaMap folder
3. We conducted molecular dynamics simulations in the paper, and example MD scripts are provided in ExampleMD

---
[\*] For parameters $k_0$, $b$, $c_{rama}$, $c_{rep}$, $c_{cyc}$, $c_{hbond}$, and $c_{other}$, if you want higher success rate of finding good backbones, you could use combinatorial design to find the optimal parameter setting, as described in the paper. To install the combinatorial design program, please check this website https://cs.nyu.edu/~shasha/papers/comb.html. You input possible values for each parameter, and run the combinatorial design program to obtain well-spaced parameter settings. If too many parameter settings are generated, pivots can be used to reduce the number of settings. Suggested pivots are $k_0$ and $b$. Finally, you run Energy_*n*residue_parameter.m to find the optimal parameter setting.
