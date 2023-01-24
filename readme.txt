This code can entirely reprocude the results of https://www.biorxiv.org/content/10.1101/2023.01.23.525256v1, including figures and tables.

Results from the applied examples can be obtained by running:
fit_ungrouped_animal.R
fit_mdnd.R
fit_crushKO.R

Results from simulation studies can be reproduced by first simulating data via MDND_Sim_Data.R, then running sbatch jobArr.sh

Note that simulations were performed using SLURM on a computing cluster. To run simulations locally (not recommended), MDND_Sim_Fit.R can be altered and run directly.

Fits used for figure generation in paper are included in Out/Fits

Note that all files in Data/ and Fits/ are tracked via LFS and should be retrieved accordingly.
