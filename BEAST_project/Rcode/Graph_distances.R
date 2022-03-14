library(dichromat)

# calculated distances between the two compartments over time 

# create matrix to save all values per generation: 

n = 100 # number of replications
g = 1000 # number of generations

data = 'distance_paper_'

# neutral without migration 
raw_distances_neutral_0 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_Neutral_distance_paper_m0','n',i-1, 'mP100' ,sep="")
  
  # Read data from file -- lemey type
  raw_distances_neutral_0[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# neutral with migration 
raw_distances_neutral_01 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_Neutral_distance_paper_m0.001','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_neutral_01[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

raw_distances_neutral_1 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_Neutral_distance_paper_m0.01','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_neutral_1[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# MFED without migration
raw_distances_MFED_0 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_MFED_distance_paper_m0','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_MFED_0[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# MFED with migration 
raw_distances_MFED_1 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_MFED_distance_paper_m0.01','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_MFED_1[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

raw_distances_MFED_01 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_MFED_distance_paper_m0.001','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_MFED_01[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# CTL without migration
raw_distances_CTL_0 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_CTL_distance_paper_m0','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_CTL_0[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# CTL with migration 
raw_distances_CTL_1 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_CTL_distance_paper_m0.01','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_CTL_1[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

raw_distances_CTL_01 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_CTL_distance_paper_m0.001','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_CTL_01[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# CTL - MFED without migration
raw_distances_CTL_MFED_0 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_CTL_MFED_distance_paper_m0','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_CTL_MFED_0[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# CTL - MFED with migration 
raw_distances_CTL_MFED_1 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_CTL_MFED_distance_paper_m0.01','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_CTL_MFED_1[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

raw_distances_CTL_MFED_01 = matrix(ncol=g, nrow=n)
for (i in seq(1,n,1)){
  name = paste('/users/judith/polybox/BEAST_project/distance/distance_CTL_MFED_distance_paper_m0.001','n',i-1,'mP100',sep="")
  
  # Read data from file -- lemey type
  raw_distances_CTL_MFED_01[i,] = read.table(name, header=FALSE, fill = TRUE)[1:g,1]
}

# caculate mean distance per generation: 

distance_neutral_0 = apply(raw_distances_neutral_0, 2, mean)
distance_neutral_1 = apply(raw_distances_neutral_1, 2, mean)
distance_neutral_01 = apply(raw_distances_neutral_01, 2, mean)
distance_MFED_0 = apply(raw_distances_MFED_0, 2, mean)
distance_MFED_1 = apply(raw_distances_MFED_1, 2, mean)
distance_MFED_01 = apply(raw_distances_MFED_01, 2, mean)
distance_CTL_0 = apply(raw_distances_CTL_0, 2, mean)
distance_CTL_1 = apply(raw_distances_CTL_1, 2, mean)
distance_CTL_01 = apply(raw_distances_CTL_01, 2, mean)
distance_CTL_MFED_0 = apply(raw_distances_CTL_MFED_0, 2, mean)
distance_CTL_MFED_1 = apply(raw_distances_CTL_MFED_1, 2, mean)
distance_CTL_MFED_01 = apply(raw_distances_CTL_MFED_01, 2, mean)

x = seq(1,g,1)

# plot lines
pdf(file='/Users/judith/polybox/BEAST_project/Figures/Distances_main_p100.pdf', width = 12, height = 10)
colors = c('red','blue','green','orange')
par(mfrow=c(1,1))
plot(x, distance_neutral_0, 'l', col=colors[1], ylab="Distance between two compartments", xlab = 'generations',cex.lab=1.5, cex.axis=1.5, lwd=3, 
     ylim=c(0,1500))
lines(x, distance_neutral_1, col=colors[1], lty=3, lwd=3.5)
lines(x, distance_neutral_01, col=colors[1], lty=5, lwd=3.5)
lines(x[1:100], distance_MFED_0[1:100], col=colors[2], lwd=3.5)
lines(x, distance_MFED_1, col=colors[2],lty=3, lwd=3.5)
lines(x, distance_MFED_01, col=colors[2],lty=5, lwd=3.5)
lines(x, distance_CTL_0, col=colors[3], lwd=3.5)
lines(x, distance_CTL_1, col=colors[3],lty=3, lwd=3.5)
lines(x, distance_CTL_01, col=colors[3],lty=5, lwd=3.5)
lines(x[1:100], distance_CTL_MFED_0[1:100], col = colors[4], lwd=3.5)
lines(x, distance_CTL_MFED_1,col=colors[4], lty=3, lwd=3.5)
lines(x, distance_CTL_MFED_01,col=colors[4], lty=5, lwd=3.5)

legend('topleft', c('Neutral no mig', 'Neutral mig=0.01', 'Neutral mig=0.001', 'MFED no mig','MFED mig=0.01', 'MFED mig=0.001', 'CTL no mig', 'CTL mig=0.01', 
                    'CTL mig=0.001', 'CTL + MFED no mig', 'CTL + MFED mig=0.01', 'CTL + MFED mig=0.001'), 
       col = c(colors[1],colors[1],colors[1], colors[2], colors[2], colors[2], colors[3], colors[3], colors[3],colors[4],colors[4],colors[4]),
       lty=c(1,3,4,1,3,4,1,3,4) , lwd=2)
dev.off()
