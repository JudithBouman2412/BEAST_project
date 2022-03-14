## Read libraries
library('ape')
library("phylobase")
#library(network)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
readFile = args[1]
treeFile = args[2]
SS_tips = args[3]

#readFile = '/Users/judit/Desktop/Solve_root/tree_Solve_root_m0.001n0P20' # file where data is saved
#treeFile = '/Users/judit/Desktop/Solve_root/solve_root_3.trees' # file where results are saved
#SS_tips = '/Users/judit/Desktop/Solve_root/IDs_Solve_root_m0.001n0P20'

# Read data of mothers from simulation
moms = as.matrix(read.table(readFile, header=FALSE, fill = TRUE))
moms[1,]='init'
moms_v = c(t(as.matrix(moms)))

# Calculate number of nodes = total number of sequences in the simulation
nseq = dim(moms)[2]/2
ngen = dim(moms)[1]
ncomp = 2
N = nseq*ngen*ncomp

# list of node names
node_names = c('init')
names_pop0 = rep('',nseq)
names_pop1 = rep('',nseq)

for (gen in seq(1,ngen,1)){ # can I make this more efficient? 
  for (seq in seq(1,nseq,1)){
    names_pop0[seq] = paste0(seq-1,'_pop0_',gen)
    names_pop1[seq] = paste0(seq-1,'_pop1_',gen)
  }
  node_names = c(node_names,names_pop0,names_pop1)
}

# make matrix with names per generation per compartment
node_names_matrix = matrix(node_names[2:(N+1)], ncol = nseq*ncomp, byrow=TRUE)

# get the tips that were sampled
tips_keep = read.table(SS_tips, header=FALSE, fill = TRUE, colClasses = "character")[,1]
n_tips = length(tips_keep)

# matrix with edges in tree. 
edge_matrix = matrix(0, nrow=ngen*n_tips, ncol=2)

for (i in seq(1,length(tips_keep),1)){
  # first tip 
  head = tips_keep[i]
  for (j in rev(seq(1,ngen,1))){
    if (strsplit(head,'_')[[1]][2]=='pop0'){
      tail = moms[j,as.numeric(strsplit(head,'_')[[1]][1])+1]
    } else if (strsplit(head,'_')[[1]][2]=='pop1'){
      tail = moms[j,as.numeric(strsplit(head,'_')[[1]][1])+nseq+1]
    }
    
    # create the edge of this connection
    edge_matrix[(i-1)*ngen+j,2] = head
    edge_matrix[(i-1)*ngen+j,1] = tail 
    
    head = tail 
    
  }
}

# remove edges that are mentioned multiple times
edge_unique <- unique(edge_matrix)

# get names of internal nodes 
internal_nodes <- unique(edge_matrix[,1])
# get number of internal nodes
n_int <- length(internal_nodes)

# all names, tips then internal nodes
nodes_tree_names = c(tips_keep,internal_nodes) # names of the nodes in the tree in correct order

# transfer node name to a number corresponding to the order of the 'nodes_tree_names'
nodes_tree_numbers = seq(1,length(nodes_tree_names),1)

edge_num = matrix(0, nrow=dim(edge_unique)[1], ncol=2)
for (j in seq(1,dim(edge_unique)[1],1)){
  edge_num[j,1] <- which(edge_unique[j,1]==nodes_tree_names)
  edge_num[j,2] <- which(edge_unique[j,2]==nodes_tree_names)
}

edge_length <- rep(1, dim(edge_num)[1])

edge_num <- rbind(edge_num, c(0,n_tips+1)) # DO we need this? 

# get location of internal nodes
locations <- rep('',length(internal_nodes))
locations[1] <- 'init'
for (i in seq(2,length(internal_nodes),1)){
  locations[i] <- strsplit(internal_nodes[i],'_')[[1]][2]
}

# Create tree
my_tree <- phylo4(edge_num, node.labels <- locations,  edge.length=c(rep(1,dim(edge_unique)[1]),0), tip.label=tips_keep) # , edge.length=rep(1,dim(edge_unique)[1])

# converting to phylo class
my_phylo <- as(my_tree,"phylo")

#collapse singleton nodes
my_new <- collapse.singles(my_phylo,root.edge = TRUE)

# Write tree to newick file
write.tree(my_new, file = treeFile, append = FALSE,
           digits = 3, tree.names = FALSE)


