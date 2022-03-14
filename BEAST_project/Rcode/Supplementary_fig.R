#### Figure on results with p=100 for Vaughan's structured coalescent model ####

library(easyGgplot2)
library(devtools)
library(ggplot2)
library(magrittr)
library(dplyr)
library(gridExtra)
library(reshape2)
library(grid)
#install.packages("patchwork")
library(patchwork)

####### Show results from n =0 m = 0.005 to exemplify why it does not work  ########

data_vaughan = read.table('/Users/judith/polybox/BEAST_project/P100_Vaughan/Vaughan_Neutral_new_original_2_m0.005n0P100.log', header=TRUE, fill = TRUE)

# Remove first 10 % of samples

size = dim(data_vaughan)[1]
data_vaughan <- data_vaughan[(size/10):size,]

# take as dataframe
data_vaughan <- as.data.frame(data_vaughan)


# Create ggplot for both migration rates 

plot_a <- ggplot(data_vaughan ) + geom_density(aes(x = migModel.t.fasta_Neutral_fullrun_2_m0.rateMatrix_pop0_pop1)) + 
  theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank()) +
  geom_vline(xintercept = 0.005, col = 'red') +
  xlab('Migration rate (0 -> 1)')  + labs(tag = "A")

plot_b <- ggplot(data_vaughan) + geom_density(aes(x = migModel.t.fasta_Neutral_fullrun_2_m0.rateMatrix_pop1_pop0)) + 
  theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank()) +
  geom_vline(xintercept = 0.005, col = 'red') +
  xlab('Migration rate (1 -> 0)')  + labs(tag = "B")


# Save plot as pdf

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/Suppl_vaughan.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 8)
grid.arrange(plot_a , plot_b, nrow=1, ncol=2)
dev.off()

n = 10 
M=c('0.001','0.002','0.003', '0.004', '0.005','0.006','0.007', '0.008','0.009','0.01') 
measure = 'migration'
method='absolute'
xlabel = "Simulated migration rate"
ylabel = "Estimated migration rate"
ylimit=c(0,0.03)

####### Supplement 2: P1000 ##############

logfilenames_with = matrix(0,
                           nrow = length(M),
                           ncol = n,
                           byrow = FALSE)
logfilenames_without = matrix(0,
                              nrow = length(M),
                              ncol = n,
                              byrow = FALSE)

############ neutral results #####
log_all_neutral_with <- data.frame(matrix(ncol = 5, nrow = 0))
names(log_all_neutral_with) = c("traitClockRate","real_mig","tree_status","x_axis",'treeStatus')       

log_all_neutral_without <- data.frame(matrix(ncol = 5 , nrow = 0))
names(log_all_neutral_without) = c("traitClockRate","real_mig","tree_status","x_axis",'treeStatus') 

means_neutral_with <- matrix(ncol=n, nrow = length(M))
var_neutral_with <- matrix(ncol=n, nrow = length(M))
means_neutral_without <- matrix(ncol=n, nrow = length(M))
var_neutral_without <- matrix(ncol=n, nrow = length(M))

for (i in 1:length(M)){
  for (j in 1:n){
    logfilenames_without[i,j] = paste('/users/judith/polybox/BEAST_project/P1000_forPaper/Lemey_Neutral_suppl_1000_m',M[i],'n',j-1,'P1000.log',sep="")
    
    # Read data from file -- lemey type
    log_tot_without = read.table(logfilenames_without[i,j], header=TRUE, fill = TRUE)
    
    # Remove first 10 %
    log_trim_without = log_tot_without[((dim(log_tot_without)[1]-1)/10+1):dim(log_tot_without)[1],]
    log_trim_frame_without = as.data.frame(log_trim_without)
    log_trim_frame_without["real_mig"] <- rep(as.numeric(M[i]),dim(log_trim_without)[1])
    
    log_trim_frame_without["treeStatus"]<- rep(0,dim(log_trim_frame_without)[1])
    log_trim_frame_without <- log_trim_frame_without %>% mutate(x_axis = as.factor(as.numeric(factor(real_mig))*100 + 10*treeStatus))
    
    levels(log_trim_frame_without$x_axis)<- paste(M[i],'b',sep='')
    
    # calculate relative difference between refered and real migration rate
    if (measure == 'migration' ){
      if (method == 'absolute'){
        log_trim_frame_without["rel"] <- (log_trim_frame_without$traitClockRate )
      } else {
        log_trim_frame_without["rel"] <- (log_trim_frame_without$traitClockRate - log_trim_frame_without$real_mig)/log_trim_frame_without$real_mig
      }
    } else if (measure == 'tree'){
      log_trim_frame_without["rel"] <- log_trim_frame_without$popSize
    } else if (measure == 'mut'){
      log_trim_frame_without["rel"] <- log_trim_frame_without$rateStat.Scenario_1.mean
    }
    
    
    means_neutral_without[i,j] <- mean(log_trim_frame_without$rel)
    var_neutral_without[i,j] <- var(log_trim_frame_without$rel)
    
    
    # combine two data frames
    log_all_neutral_without_final <- cbind(log_trim_frame_without$rel, log_trim_frame_without$real_mig, 
                                           log_trim_frame_without$treeStatus, log_trim_frame_without$x_axis) 
    log_all_neutral_without <- rbind(log_all_neutral_without,log_all_neutral_without_final)
  }
}

names(log_all_neutral_without) <- c('traitClockRate','real_mig','treeStatus','x_axis')
error_neutral_without <- data.frame(matrix(ncol = 4 , nrow = length(M)))
names(error_neutral_without) = c("mean","sd","upper","lower")
error_neutral_without["real_mig"]=as.factor(M)

for (i in 1:length(M) ){
  y = log_all_neutral_without$traitClockRate[log_all_neutral_without$real_mig==as.numeric(M[i])]
  error_neutral_without[i,1] = median(y)
  error_neutral_without[i,2] = sd(y)
  error_neutral_without[i,3] = sort(y)[round(length(y)/100*95)]
  error_neutral_without[i,4] = sort(y)[round(length(y)/100*5)]
}

#Melting data neutral
M2_neutral <- melt(log_all_neutral_without,id.vars= 2, measure.vars= 1 )
M2_neutral$treeStatus <- log_all_neutral_without[,3]
M2_neutral <- M2_neutral %>% mutate(x_axis = as.factor(as.numeric(factor(real_mig))*100 + 10*treeStatus))
levels(M2_neutral$x_axis)<- paste(M,'b',sep='')

M2_neutral$real_mig <- as.factor(M2_neutral$real_mig)
M2_neutral['mig'] <- as.numeric(M2_neutral$real_mig)/1000

error_neutral_without$x_axis <- levels(M2_neutral$x_axis)

for (i in 1:length(M)){
  for (j in 1:n){
    logfilenames_with[i,j] = paste('/users/judith/polybox/BEAST_project/P1000_forPaper/Lemey_Neutral_suppl_1000_m',M[i],'n',j-1,'P1000.log',sep="")
    
    # Read data from file -- lemey type
    log_tot_with = read.table(logfilenames_with[i,j], header=TRUE, fill = TRUE)
    
    # Remove first 10 %
    log_trim_with = log_tot_with[((dim(log_tot_with)[1]-1)/10+1):dim(log_tot_with)[1],]
    log_trim_frame_with = as.data.frame(log_trim_with)
    log_trim_frame_with["real_mig"] <- rep(as.numeric(M[i]),dim(log_trim_with)[1])
    
    log_trim_frame_with["treeStatus"] <- rep(1,dim(log_trim_frame_with)[1])
    log_trim_frame_with <- log_trim_frame_with %>% mutate(x_axis = as.factor(as.numeric(factor(real_mig))*100 + 10*treeStatus))
    
    levels(log_trim_frame_with$x_axis)<- paste(M[i],'a',sep='')
    
    # calculate relative difference between referred and real migration rate
    if (measure == 'migration' ){
      if (method == 'absolute'){
        log_trim_frame_with["rel"] <- log_trim_frame_with$traitClockRate
      } else {
        log_trim_frame_with["rel"] <- (log_trim_frame_with$traitClockRate - log_trim_frame_with$real_mig)/log_trim_frame_with$real_mig
      }
    } else if (measure == 'tree'){
      log_trim_frame_with["rel"] <- log_trim_frame_with$popSize
    } else if (measure == 'mut'){
      log_trim_frame_with["rel"] <- log_trim_frame_with$rateStat.Scenario_1.mean
    }
    
    
    means_neutral_with[i,j] <- mean(log_trim_frame_with$rel)
    var_neutral_with[i,j] <- var(log_trim_frame_with$rel)
    
    # combine two data frames
    log_all_neutral_with_final <- cbind(log_trim_frame_with$rel, log_trim_frame_with$real_mig, 
                                        log_trim_frame_with$treeStatus, log_trim_frame_with$x_axis)
    
    log_all_neutral_with <- rbind(log_all_neutral_with,log_all_neutral_with_final)
  }
}

names(log_all_neutral_with) <- c('traitClockRate','real_mig','treeStatus','x_axis')
error_neutral_with <- data.frame(matrix(ncol = 4 , nrow = length(M)))
names(error_neutral_with) = c("mean","sd","upper","lower")
error_neutral_with["real_mig"]=as.factor(M)

for (i in 1:length(M) ){
  x = log_all_neutral_with$traitClockRate[log_all_neutral_with$real_mig==as.numeric(M[i])]
  error_neutral_with[i,1] = median(x)
  error_neutral_with[i,2] = sd(x)
  error_neutral_with[i,3] = sort(x)[round(length(x)/100*95)]
  error_neutral_with[i,4] = sort(x)[round(length(x)/100*5)]
}

#Melting data neutral
M1_neutral <- melt(log_all_neutral_with[,1:2],id.vars= 'real_mig', variable.name = 'traitClockRate')
M1_neutral$treeStatus <- log_all_neutral_with[,3]
M1_neutral <- M1_neutral %>% mutate(x_axis = as.factor(as.numeric(factor(real_mig))*100 + 10*treeStatus))
levels(M1_neutral$x_axis)<- paste(M,'a',sep='')

M1_neutral$real_mig <- as.factor(M1_neutral$real_mig)
M1_neutral['mig'] <- as.numeric(M1_neutral$real_mig)/1000

error_neutral_with$x_axis <- levels(M1_neutral$x_axis)


  basis <- data.frame(matrix(ncol = 4 , nrow = 2*length(M)))
  names(basis) = c("traitClockRate","real_mig","tree_status","x_axis")    
  basis['traitClockRate']=as.numeric(c('0.001','0.001', '0.002', '0.002','0.003','0.003', '0.004', '0.004', '0.005', '0.005','0.006','0.006', '0.007', '0.007','0.008','0.008','0.009','0.009','0.01','0.01'))
  basis['real_mig']=as.numeric(c('0.001','0.001', '0.002', '0.002','0.003','0.003', '0.004', '0.004', '0.005', '0.005','0.006','0.006', '0.007', '0.007','0.008','0.008','0.009','0.009','0.01','0.01'))
  basis['x_axis']=c('0.001a','0.001b', '0.002a', '0.002b','0.003a','0.003b', '0.004a', '0.004b', '0.005a', '0.005b','0.006a','0.006b', '0.007a', '0.007b','0.008a','0.008b','0.009a','0.009b','0.01a','0.01b')
  
  basis.M <- melt(basis[,1:2],id.vars= 'real_mig', variable.name = 'traitClockRate')
  basis.M$real_mig <- as.factor(basis.M$real_mig)
  
  # Create Figure
  if (measure == 'migration') {
    een = as.numeric(c(rep(0.001,2), rep(0.002,2),rep(0.003,2),rep(0.004,2),rep(0.005,2),rep(0.006,2),rep(0.007,2),rep(0.008,2),rep(0.009,2),rep(0.01,2)))
  } else if (measure == 'tree'){
    een = rep(200, length(M)*2)
  } else if (measure == 'mut'){
    een = rep(2.16e-4, length(M)*2)
  }
  
  drie = c('0.001a','0.001b', '0.002a', '0.002b','0.003a','0.003b', '0.004a', '0.004b', '0.005a', '0.005b','0.006a','0.006b', '0.007a', '0.007b','0.008a','0.008b','0.009a','0.009b','0.01a','0.01b')
  true.rates = as.data.frame(cbind(een,drie))
  true.rates$een <- as.numeric(true.rates$een)
  
  p1 =  
    ggplot() +
    geom_violin(data=M1_neutral, aes(x = x_axis, y = value),colour="#1A2C56", fill="#1A2C56", alpha = .2, size=0.3) +
    geom_violin(data=M2_neutral, aes(x = x_axis, y = value),colour="#D1A683",fill="#D1A683", alpha = .2, size=0.3) +
    xlab(xlabel)+ ylab(ylabel) + ylim(ylimit) +
    ggtitle('Neutral') + 
    stat_summary(data=error_neutral_with,aes(x = x_axis, y = mean,colour="Using true genealogy"), geom="point", size=3) +
    stat_summary(data=error_neutral_without,aes(x = x_axis, y = mean,colour="Using sequence data"), geom="point", size=3) +
    geom_errorbar(data=error_neutral_with,aes(x=x_axis, ymax=upper, ymin=lower),colour="#1A2C56",size=2) +
    geom_errorbar(data=error_neutral_without,aes(x=x_axis, ymax=upper, ymin=lower),size=2,colour="#D1A683") +
    theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
          panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(), legend.background = element_blank()) +
    geom_point(data = true.rates, aes(x=drie, y = een), col = 'red') + 
    scale_x_discrete( labels= c('0.001','', '0.002', '','0.003','', '0.004', '', '0.005', '','0.006','', '0.007', '','0.008','', '0.009','','0.01','') )+
    scale_colour_manual("", values=c("Using True Genealogy" = "#1A2C56", "Using sequence data" ="#D1A683")) 

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/Suppl_P1000.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8)
p1
dev.off()

