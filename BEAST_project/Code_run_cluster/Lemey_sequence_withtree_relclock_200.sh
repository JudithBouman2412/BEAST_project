#M=(0.00001 0.0001 0.001 0.01 0.1)
#M=(0.002 0.004 0.006 0.008 0.02 0.04 0.06 0.08)
#M=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
#M=(0.02 0.03 0.04 0.06 0.07 0.08 0.09)
#M=(0.02)

M=(0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01)
mu=(0.000216)
S=100
n=20
P=(100)
G=(1000)
s1=(100)
s2=(100)

name=original_2

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

XMLNames_Lemey_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/Lemey_Neutral_with_relclock${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_Neutral[i,$j]=/cluster/scratch/jbouman/P100_dec2021/Lemey_Neutral_with_relclock_${name}_m${M[i]}n${j}P${P[g]}

XMLNames_Lemey_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/Lemey_MFED_with_relclock_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/Lemey_MFED_with_relclock_${name}_m${M[i]}n${j}P${P[g]}

XMLNames_Lemey_CTL[i,$j]=/cluster/home/jbouman/P100_dec2021/Lemey_CTL_with_relclock_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_CTL[i,$j]=/cluster/scratch/jbouman/P100_dec2021/Lemey_CTL_with_relclock_${name}_m${M[i]}n${j}P${P[g]}

XMLNames_Lemey_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/Lemey_CTL_MFED_with_relclock_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_CTL_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/Lemey_CTL_MFED_with_relclock_${name}_m${M[i]}n${j}P${P[g]}


# Neutral without CTL
/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Template_relaxed_clock_exp_Tim_Tree.xml --fastaFile ${FastaNames_Neutral[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_Neutral[i,$j]} --age ${IDs_Neutral[i,$j]} --MigrationModelType Lemey_tree --logFileBasename ${LogNames_Lemey_Neutral[i,$j]} --dateDirection forward --chainLength 10000000 --treeFile ${newickFile_Neutral[i,$j]}

bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_Neutral[i,$j]}"

# MFED
/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Template_relaxed_clock_exp_Tim_Tree.xml --fastaFile ${FastaNames_MFED[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_MFED[i,$j]} --age ${IDs_MFED[i,$j]} --MigrationModelType Lemey_tree --logFileBasename ${LogNames_Lemey_MFED[i,$j]} --dateDirection forward --chainLength 100000000 --treeFile ${newickFile_MFED[i,$j]}

bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_MFED[i,$j]}"

# CTL
/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Template_relaxed_clock_exp_Tim_Tree.xml --fastaFile ${FastaNames_CTL[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_CTL[i,$j]} --age ${IDs_CTL[i,$j]} --MigrationModelType Lemey_tree --logFileBasename ${LogNames_Lemey_CTL[i,$j]} --dateDirection forward --chainLength 100000000 --treeFile ${newickFile_CTL[i,$j]}

bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_CTL[i,$j]}"

# CTL MFED
/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Template_relaxed_clock_exp_Tim_Tree.xml --fastaFile ${FastaNames_CTL_MFED[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_CTL_MFED[i,$j]} --age ${IDs_CTL_MFED[i,$j]} --MigrationModelType Lemey_tree --logFileBasename ${LogNames_Lemey_CTL_MFED[i,$j]} --dateDirection forward --chainLength 100000000 --treeFile ${newickFile_CTL_MFED[i,$j]}

bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_CTL_MFED[i,$j]}"



done
done
done
done
