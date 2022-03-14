# Making the trees
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

declare XMLNames_Lemey
declare LogNames_Lemey
declare EpiFiles

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
newickFile_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/newick_Neutral_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_MFED_${name}_m${M[i]}n${j}mP${P[g]}
newickFile_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/newick_MFED_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_CTL_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_CTL_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_CTL_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_CTL_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_CTL_${name}_m${M[i]}n${j}mP${P[g]}
newickFile_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/newick_CTL_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/fasta_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/IDs_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/epis_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/tree_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/n_migration_CTL_MFED_${name}_m${M[i]}n${j}mP${P[g]}
newickFile_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/newick_CTL_MFED_${name}_m${M[i]}n${j}mP${P[g]}

# Neutral without CTL

bsub -W 120:00 -R "rusage[mem=10000]" "Rscript --vanilla /cluster/home/jbouman/Tree_backwards.R ${treeFile_Neutral[i,$j]} ${newickFile_Neutral[i,$j]} ${IDs_Neutral[i,$j]}"

#MFED substitution model

bsub -W 120:00 -R "rusage[mem=10000]" "Rscript --vanilla /cluster/home/jbouman/Tree_backwards.R ${treeFile_MFED[i,$j]} ${newickFile_MFED[i,$j]} ${IDs_MFED[i,$j]}"

# CTL with neutral background substitution model

bsub -W 120:00 -R "rusage[mem=10000]" "Rscript --vanilla /cluster/home/jbouman/Tree_backwards.R ${treeFile_CTL[i,$j]} ${newickFile_CTL[i,$j]} ${IDs_CTL[i,$j]}"

#MFED substitution model

bsub -W 120:00 -R "rusage[mem=10000]" "Rscript --vanilla /cluster/home/jbouman/Tree_backwards.R ${treeFile_CTL_MFED[i,$j]} ${newickFile_CTL_MFED[i,$j]} ${IDs_CTL_MFED[i,$j]}"


done
done
done
done
