M=(0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01)
mu=(0.000216)
S=100
n=10
P=(100)
G=(1000)
s1=(100)
s2=(100)

name=rooted

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

FastaNames_CTL[i,$j]=/cluster/home/jbouman/P100_june2020/fasta_CTL_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL[i,$j]=/cluster/home/jbouman/P100_june2020/IDs_CTL_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL[i,$j]=/cluster/home/jbouman/P100_june2020/epis_CTL_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL[i,$j]=/cluster/scratch/jbouman/P100_june2020/tree_CTL_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL[i,$j]=/cluster/home/jbouman/P100_june2020/n_migration_CTL_${name}_m${M[i]}n${j}mP${P[g]}

FastaNames_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_june2020/fasta_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}.fasta
IDs_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_june2020/IDs_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
Epi_files_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_june2020/epis_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
treeFile_CTL_MFED[i,$j]=/cluster/scratch/jbouman/P100_june2020/tree_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}
migrationFile_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_june2020/n_migration_CTL_MFED_${name}_m${M[i]}n${j}mP${P[g]}


XMLNames_Lemey_Neutral[i,$j]=/cluster/home/jbouman/P100_dec2021/Vaughan_Neutral_new_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_Neutral[i,$j]=/cluster/scratch/jbouman/P100_dec2021/Vaughan_Neutral_new_${name}_m${M[i]}n${j}P${P[g]}

XMLNames_Lemey_MFED[i,$j]=/cluster/home/jbouman/P100_dec2021/Vaughan_MFED_new_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_MFED[i,$j]=/cluster/scratch/jbouman/P100_dec2021/Vaughan_MFED_new_${name}_m${M[i]}n${j}P${P[g]}

XMLNames_Lemey_CTL[i,$j]=/cluster/home/jbouman/P100_june2020/Vaughan_CTL_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_CTL[i,$j]=/cluster/scratch/jbouman/P100_june2020/Vaughan_CTL_${name}_m${M[i]}n${j}P${P[g]}

XMLNames_Lemey_CTL_MFED[i,$j]=/cluster/home/jbouman/P100_june2020/Vaughan_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}.xml
LogNames_Lemey_CTL_MFED[i,$j]=/cluster/scratch/jbouman/P100_june2020/Vaughan_CTL_MFED_${name}_m${M[i]}n${j}P${P[g]}


#Neutral
/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Vaughan_template.xml --fastaFile ${FastaNames_Neutral[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_Neutral[i,$j]} --age ${IDs_Neutral[i,$j]} --MigrationModelType Vaughan --logFileBasename ${LogNames_Lemey_Neutral[i,$j]} --dateDirection forward --chainLength 100000000 --mutationRate 0.000216

bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_Neutral[i,$j]}"

#MFED
#/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Vaughan_template.xml --fastaFile ${FastaNames_MFED[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_MFED[i,$j]} --age ${IDs_MFED[i,$j]} --MigrationModelType Vaughan --logFileBasename ${LogNames_Lemey_MFED[i,$j]} --dateDirection forward --chainLength 100000000 --mutationRate 0.0000216

#bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_MFED[i,$j]}"

#CTL
#/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Vaughan_template.xml --fastaFile ${FastaNames_CTL[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_CTL[i,$j]} --age ${IDs_CTL[i,$j]} --MigrationModelType Vaughan --logFileBasename ${LogNames_Lemey_CTL[i,$j]} --dateDirection forward --chainLength 100000000 --mutationRate 0.000216

#bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_CTL[i,$j]}"

#CTL MFED
#/cluster/home/jbouman/Beast/beast2-xml.py --templateFile /cluster/home/jbouman/templates/Vaughan_template.xml --fastaFile ${FastaNames_CTL_MFED[i,$j]} --ageIncluded --XML_FILENAME ${XMLNames_Lemey_CTL_MFED[i,$j]} --age ${IDs_CTL_MFED[i,$j]} --MigrationModelType Vaughan --logFileBasename ${LogNames_Lemey_CTL_MFED[i,$j]} --dateDirection forward --chainLength 100000000 --mutationRate 0.000216

#bsub -W 120:00 -R "rusage[mem=10000]" java -Xmx2G -Dbeast.load.jars=true -jar "/cluster/home/jbouman/Beast/beast_2.5.0/lib/beast.jar" -overwrite "${XMLNames_Lemey_CTL_MFED[i,$j]}"

done
done
done
done
