# Simulate data

M=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
#M=(0.02)
mu=(0.0000216)
S=100
n=11
P=(1000)
G=(1000)
s1=(100)
s2=(100)

name=full_2

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

FastaNames_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_Neutral_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_Neutral_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_Neutral_${name}_m${M[i]}n${j}P${P[g]}
treeFile_Neutral[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_Neutral_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_Neutral_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_MFED_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_CTL_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_CTL_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_CTL_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_CTL_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_CTL_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_CTL_MFED_${name}_m${M[i]}n${j}mP${P[g]}

#Neutral
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G}  -sa $s1 $s2 -fastaFile ${FastaNames_Neutral[i,$j]} -names pop0 pop1 -o neutral_high -IDfile ${IDs_Neutral[i,$j]} -limit 1000 -s 1000 -epis_file ${Epi_files_Neutral[i,$j]} -treeFile ${treeFile_Neutral[i,$j]} -migrationFile ${migrationFile_Neutral[i,$j]} -mut ${mu} ${mu}

#MFE
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G} -sa $s1 $s2 -fastaFile ${FastaNames_MFED[i,$j]} -names pop0 pop1 -o MFED_equal -IDfile ${IDs_MFED[i,$j]} -limit 1000 -s 1000 -epis_file ${Epi_files_MFED[i,$j]} -treeFile ${treeFile_MFED[i,$j]} -migrationFile ${migrationFile_MFED[i,$j]} -mut ${mu} ${mu}

#CTL
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G} -sa $s1 $s2 -fastaFile ${FastaNames_CTL[i,$j]} -names pop0 pop1 -o neutral_high -IDfile ${IDs_CTL[i,$j]} -limit 1000 -s 1 -epis_file ${Epi_files_CTL[i,$j]} -treeFile ${treeFile_CTL[i,$j]} -migrationFile ${migrationFile_CTL[i,$j]} -mut ${mu} ${mu}

#MFED & CTL
bsub -W 20:00 -R "rusage[mem=10000]" python /cluster/home/jbouman/SimTool/sim_scripts/multiple_compartments.py -ng ${G}  -ninit ${P[g]} ${P[g]} -nc 2 -mig [[${M[i]},${M[i]}],[${M[i]},${M[i]}]] -R0 6 6 -maxpop ${P[g]} ${P[g]} -st ${G} -sa $s1 $s2 -fastaFile ${FastaNames_CTL_MFED[i,$j]} -names pop0 pop1 -o MFED_equal -IDfile ${IDs_CTL_MFED[i,$j]} -limit 1000 -s 1 -epis_file ${Epi_files_CTL_MFED[i,$j]} -treeFile ${treeFile_CTL_MFED[i,$j]} -migrationFile ${migrationFile_CTL_MFED[i,$j]} -mut ${mu} ${mu}


done
done
done 
done



