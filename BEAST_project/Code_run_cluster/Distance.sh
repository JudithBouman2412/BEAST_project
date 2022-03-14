n=100
M=(0 0.05)
mu=(0.0000216)
S=100
P=(1000)
G=(1000)
s1=(100)
s2=(100)

name=distance

# loop over replicates
for ((j=0;j<n ;j+=1));
do
# loop over migration rates
for i in ${!M[*]};
do
# loop over sample size
for k in ${!S[*]};
do
# loop over population size
for g in ${!P[*]};
do

FastaNames_Neutral[i,$j]=/cluster/scratch/jbouman/distance/fasta_Neutral_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_Neutral[i,$j]=/cluster/scratch/jbouman/distance/IDs_Neutral_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_Neutral[i,$j]=/cluster/scratch/jbouman/distance/epis_Neutral_${name}_m${M[i]}n${j}P${P[g]}
treeFile_Neutral[i,$j]=/cluster/scratch/jbouman/distance/tree_Neutral_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_Neutral[i,$j]=/cluster/scratch/jbouman/distance/n_migration_Neutral_${name}_m${M[i]}n${j}mP${P[g]}
distanceFile_Neutral[i,$j]=/cluster/home/jbouman/distance/distance_Neutral_${name}_m${M[i]}n${j}mP${P[g]}

#Neutral
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng 500  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st 500  -sa $s1 $s2 -fastaFile ${FastaNames_Neutral[i,$j]} -names pop0 pop1 -o neutral_high -IDfile ${IDs_Neutral[i,$j]} -limit 1000 -s 1000 -epis_file ${Epi_files_Neutral[i,$j]} -treeFile ${treeFile_Neutral[i,$j]} -migrationFile ${migrationFile_Neutral[i,$j]} -mut ${mu} ${mu} -distanceFile ${distanceFile_Neutral[i,$j]}

FastaNames_MFED[i,$j]=/cluster/scratch/jbouman/distance/fasta_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_MFED[i,$j]=/cluster/scratch/jbouman/distance/IDs_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_MFED[i,$j]=/cluster/scratch/jbouman/distance/epis_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_MFED[i,$j]=/cluster/scratch/jbouman/distance/tree_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_MFED[i,$j]=/cluster/scratch/jbouman/distance/n_migration_MFED_${name}_m${M[i]}n${j}mP${P[g]}
distanceFile_MFED[I,$j]=/cluster/home/jbouman/distance/distance_MFED_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL[i,$j]=/cluster/scratch/jbouman/distance/fasta_CTL_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL[i,$j]=/cluster/scratch/jbouman/distance/IDs_CTL_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL[i,$j]=/cluster/scratch/jbouman/distance/epis_CTL_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL[i,$j]=/cluster/scratch/jbouman/distance/tree_CTL_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL[i,$j]=/cluster/scratch/jbouman/distance/n_migration_CTL_${name}_m${M[i]}n${j}mP${P[g]}
distanceFile_CTL[I,$j]=/cluster/home/jbouman/distance/distance_CTL_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL_MFED[i,$j]=/cluster/scratch/jbouman/distance/fasta_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL_MFED[i,$j]=/cluster/scratch/jbouman/distance/IDs_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL_MFED[i,$j]=/cluster/scratch/jbouman/distance/epis_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL_MFED[i,$j]=/cluster/scratch/jbouman/distance/tree_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL_MFED[i,$j]=/cluster/scratch/jbouman/distance/n_migration_CTL_MFED_${name}_m${M[i]}n${j}mP${P[g]}
distanceFile_CTL_MFED[I,$j]=/cluster/home/jbouman/distance/distance_CTL_MFED_${name}_m${M[i]}n${j}mP${P[g]}

#MFE
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G} -sa $s1 $s2 -fastaFile ${FastaNames_MFED[i,$j]} -names pop0 pop1 -o MFED_equal -IDfile ${IDs_MFED[i,$j]} -limit 1000 -s 1000 -epis_file ${Epi_files_MFED[i,$j]} -treeFile ${treeFile_MFED[i,$j]} -migrationFile ${migrationFile_MFED[i,$j]} -mut ${mu} ${mu} -distanceFile ${distanceFile_MFED[i,$j]}

#CTL
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G} -sa $s1 $s2 -fastaFile ${FastaNames_CTL[i,$j]} -names pop0 pop1 -o neutral_high -IDfile ${IDs_CTL[i,$j]} -limit 1000 -s 1 -epis_file ${Epi_files_CTL[i,$j]} -treeFile ${treeFile_CTL[i,$j]} -migrationFile ${migrationFile_CTL[i,$j]} -mut ${mu} ${mu} -distanceFile ${distanceFile_CTL[i,$j]}

#MFED & CTL
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G} -sa $s1 $s2 -fastaFile ${FastaNames_CTL_MFED[i,$j]} -names pop0 pop1 -o MFED_equal -IDfile ${IDs_CTL_MFED[i,$j]} -limit 1000 -s 1 -epis_file ${Epi_files_CTL_MFED[i,$j]} -treeFile ${treeFile_CTL_MFED[i,$j]} -migrationFile ${migrationFile_CTL_MFED[i,$j]} -mut ${mu} ${mu} -distanceFile ${distanceFile_CTL_MFED[i,$j]}


done
done
done
done
