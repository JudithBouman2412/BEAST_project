This github folder contains all code to reproduce the results in the manuscript ''. 

The file 'Simulate.sh' from the folder 'Code_run_cluster' can be ran on a cluster where the simulation tool from Eva Bons is installed: https://gitlab.ethz.ch/eva-bons/SeqSim where file multiple_compartments.py and simulations.py are replaced by the file with the same names from the main folder here. 
Running 'Simulate.sh' will create the simulated data used to produce the main result figures from the manuscript. 

The file 'Lemey_sequence_data.sh' from the folder 'Code_run_cluster' can be ran on the cluster where BEAST2 is installed: https://www.beast2.org/ and will run the deme-migration model described by Lemey et al (2009) for all simulated data sequence data from the previous step. 

The file 'ConstructRealTree.sh' from the folder 'Code_run_cluster' will create the real phylogenetic trees from the simulations when ran on the cluster on the simulated data. This file uses the R file 'Tree_backwards.R'.  
After creating these files, one can run 'Lemey_sequence_withtree.sh' to analyse the same simulations when using the real tree in BEAST2.

The results (saved in .log files) from the BEAST2 analyses of Lemey's model is then analysed with the file 'Creat_Figures.R' from folder 'Rcode' to create the main figures in the result section of the manuscript. This file uses the functions described in the file 'Functions.R' from the 'Rcode' folder. 

We did not save all .log files created with the analysis, as saving these files requires a large amount of memory and due to the possibility to recreate these data if neccasery, we prefer to save resources by not keeping these files on a public server. 

To create the figure that shows the distance between both compartments over the number of generations run 'Distance.sh' from the folder 'Code_run_cluster'. 

We also analysed the simulated data with the relaxed clock, for this use files: 'Lemey_sequence_data_relclock_200.sh' and 'Lemey_sequence_withtree_relclock_200.sh' from the folder 'Code_run_cluster'. 

For analysing the data while keeping the population fixed, use the original code with a narrow uniform prior on the population size. 

For the supplement, we used the structured coalescent as described by Vaughan et al. (2014). This is done with the 'Vaughan_sequence_data.sh' file and corresponding template. 

