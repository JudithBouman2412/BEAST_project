# Load libraries
library(easyGgplot2)
library(devtools)
library(ggplot2)
library(magrittr)
library(dplyr)
library(gridExtra)
library(reshape2)
library(grid)

# Load file with functions: 
source('/Users/judith/polybox/BEAST_project/Rcode/Functions.R')

# Rerun?
rerun = FALSE

# mutation rates + number of replicates
M=c('0.001','0.002','0.003', '0.004', '0.005','0.006','0.007', '0.008','0.009','0.01') 


##### Main result figure: migration rate, p = 100, relaxed clock #####

if (rerun){
  out_neutral_mig_rel <- create_subFig(selection_type='Neutral', M=M, n = 15, 
                                       data = 'relclock_original_2', folder = 'P100_forPaper_relclock')
  out_MFED_mig_rel <- create_subFig(selection_type='MFED', M=M, n = 15,
                                    data = 'relclock_original_2', folder = 'P100_forPaper_relclock')
  out_CTL_mig_rel <- create_subFig(selection_type='CTL', M=M, n = 15,
                                   data = 'relclock_original_2', folder = 'P100_forPaper_relclock')
  out_CTL_MFED_mig_rel <- create_subFig(selection_type='CTL_MFED', M=M, n = 15,
                                        data = 'relclock_original_2', folder = 'P100_forPaper_relclock')
} 

p1_rel <- out_neutral_mig_rel[[1]]
p2_rel <- out_MFED_mig_rel[[1]]
p3_rel <- out_CTL_mig_rel[[1]]
p4_rel <- out_CTL_MFED_mig_rel[[1]]

mylegend<-g_legend(p1_rel)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/Main_results_p100_relclock.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p1_rel+ theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_rel+ theme(legend.position="none", axis.title=element_blank()), 
                         p3_rel+ theme(legend.position="none"),
                         p4_rel+ theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

# Sig migration rates: 

means_neutral_with_rel <- out_neutral_mig_rel[[4]]
means_neutral_without_rel <- out_neutral_mig_rel[[5]]

means_MFED_with_rel <- out_MFED_mig_rel[[4]]
means_MFED_without_rel <- out_MFED_mig_rel[[5]]

means_CTL_with_rel <- out_CTL_mig_rel[[4]]
means_CTL_without_rel <- out_CTL_mig_rel[[5]]

means_CTL_MFED_with_rel <- out_CTL_MFED_mig_rel[[4]]
means_CTL_MFED_without_rel <- out_CTL_MFED_mig_rel[[5]]

m=seq(0.001,0.01,0.001)
n = 15

sig_results_mig_rel <- create_sig_fig(means_neutral_with=means_neutral_with_rel, means_neutral_without=means_neutral_without_rel,
                                  means_MFED_with=means_MFED_with_rel,means_MFED_without=means_MFED_without_rel,
                                  means_CTL_with=means_CTL_with_rel,means_CTL_without=means_CTL_without_rel,
                                  means_CTL_MFED_with=means_CTL_MFED_with_rel,means_CTL_MFED_without=means_CTL_MFED_without_rel, m=seq(0.001,0.01,0.001), n =15)

p2_A_rel <- sig_results_mig_rel[[1]]
p2_B_rel <- sig_results_mig_rel[[2]]
p2_C_rel <- sig_results_mig_rel[[3]]
p2_D_rel <- sig_results_mig_rel[[4]]

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/bias_detection_p100_relclock.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p2_A_rel + theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_B_rel + theme(legend.position="none", axis.title=element_blank()), 
                         p2_C_rel + theme(legend.position="none"),
                         p2_D_rel + theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

table_original_rel <- rbind( abs(out_neutral_mig_rel[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_MFED_mig_rel[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_CTL_mig_rel[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_CTL_MFED_mig_rel[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100)

rownames(table_original_rel) <- c('Neutral without','MFED without', 'CTL without', 'CTL MFED without')
colnames(table_original_rel) <- seq(0.001,0.01,0.001)


##### Effective population size, p = 100, relaxed clock #####

if (rerun){
  out_neutral_tree_rel <- create_subFig(selection_type='Neutral', M=M, measure='tree', n = 15,
                                        data = 'relclock_original_2', folder = 'P100_forPaper_relclock',
                                        xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
  out_MFED_tree_rel <- create_subFig(selection_type='MFED', M=M, measure='tree', n = 15,
                                     data = 'relclock_original_2', folder = 'P100_forPaper_relclock',
                                     xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
  out_CTL_tree_rel <- create_subFig(selection_type='CTL', M=M, measure='tree', n = 15,
                                    data = 'relclock_original_2', folder = 'P100_forPaper_relclock',
                                    xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
  out_CTL_MFED_tree_rel <- create_subFig(selection_type='CTL_MFED', M=M, measure='tree',
                                         data = 'relclock_original_2', n = 15, folder = 'P100_forPaper_relclock',
                                         xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
} 

p1_tree_rel <- out_neutral_tree_rel[[1]]
p2_tree_rel <- out_MFED_tree_rel[[1]]
p3_tree_rel <- out_CTL_tree_rel[[1]]
p4_tree_rel <- out_CTL_MFED_tree_rel[[1]]

mylegend<-g_legend(p1_tree_rel)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/PopSize_results_p100_relclock.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p1_tree_rel+ theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_tree_rel+ theme(legend.position="none", axis.title=element_blank()), 
                         p3_tree_rel+ theme(legend.position="none"),
                         p4_tree_rel+ theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

# Sig population sizes rates: 

means_neutral_with_tree_rel <- out_neutral_tree_rel[[4]]
means_neutral_without_tree_rel <- out_neutral_tree_rel[[5]]

means_MFED_with_tree_rel <- out_MFED_tree_rel[[4]]
means_MFED_without_tree_rel <- out_MFED_tree_rel[[5]]

means_CTL_with_tree_rel <- out_CTL_tree_rel[[4]]
means_CTL_without_tree_rel <- out_CTL_tree_rel[[5]]

means_CTL_MFED_with_tree_rel <- out_CTL_MFED_tree_rel[[4]]
means_CTL_MFED_without_tree_rel <- out_CTL_MFED_tree_rel[[5]]

sig_results_tree_rel <- create_sig_fig(means_neutral_with=means_neutral_with_tree_rel, means_neutral_without=means_neutral_without_tree_rel,
                                      means_MFED_with=means_MFED_with_tree_rel,means_MFED_without=means_MFED_without_tree_rel,
                                      means_CTL_with=means_CTL_with_tree_rel,means_CTL_without=means_CTL_without_tree_rel,
                                      means_CTL_MFED_with=means_CTL_MFED_with_tree_rel,means_CTL_MFED_without=means_CTL_MFED_without_tree_rel, m=seq(0.001,0.01,0.001), n =15, compare = rep(200, 10))

p2_A_tree_rel <- sig_results_tree_rel[[1]]
p2_B_tree_rel <- sig_results_tree_rel[[2]]
p2_C_tree_rel <- sig_results_tree_rel[[3]]
p2_D_tree_rel <- sig_results_tree_rel[[4]]

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/bias_detection_p100_popsize_relclock.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p2_A_tree_rel + theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_B_tree_rel + theme(legend.position="none", axis.title=element_blank()), 
                         p2_C_tree_rel + theme(legend.position="none"),
                         p2_D_tree_rel + theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

table_original_tree_rel <- rbind( abs(out_neutral_tree_rel[[3]][,1]-200)/200*100,
                             abs(out_MFED_tree_rel[[3]][,1]-200)/200*100,
                             abs(out_CTL_tree_rel[[3]][,1]-200)/200*100,
                             abs(out_CTL_MFED_tree_rel[[3]][,1]-200)/200*100)

rownames(table_original_tree_rel) <- c('Neutral without','MFED without', 'CTL without', 'CTL MFED without')
colnames(table_original_tree_rel) <- seq(0.001,0.01,0.001)


##### Mutation rate estimate, p = 100, relaxed clock #####

if (rerun){
  out_neutral_mut_rel <- create_subFig(selection_type='Neutral', M=M, measure='mut', n = 15,
                                       data = 'relclock_original_2', folder = 'P100_forPaper_relclock',
                                       xlabel = "Simulated migration rate", ylabel = "Mutation rate estimate", ylimit=c(0,5e-4), full =FALSE)
  out_MFED_mut_rel <- create_subFig(selection_type='MFED', M=M, measure='mut', n = 15,
                                    data = 'relclock_original_2', folder = 'P100_forPaper_relclock',
                                    xlabel = "Simulated migration rate", ylabel = "Mutation rate estimate", ylimit=c(0,5e-4), full =FALSE)
  out_CTL_mut_rel <- create_subFig(selection_type='CTL', M=M, measure='mut', n = 15,
                                   data = 'relclock_original_2', folder = 'P100_forPaper_relclock',
                                   xlabel = "Simulated migration rate", ylabel = "Mutation rate estimate", ylimit=c(0,5e-4), full =FALSE)
  out_CTL_MFED_mut_rel <- create_subFig(selection_type='CTL_MFED', M=M, measure='mut',
                                        data = 'relclock_original_2', n = 15, folder = 'P100_forPaper_relclock',
                                        xlabel = "Simulated migration rate", ylabel = "Mutation rate estimate", ylimit=c(0,5e-4), full =FALSE)
} 

p1_mut_rel <- out_neutral_mut_rel[[1]]
p2_mut_rel <- out_MFED_mut_rel[[1]]
p3_mut_rel <- out_CTL_mut_rel[[1]]
p4_mut_rel <- out_CTL_MFED_mut_rel[[1]]

mylegend<-g_legend(p1_mut_rel)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/mut_results_p100_relclock.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p1_mut_rel+ theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_mut_rel+ theme(legend.position="none", axis.title=element_blank()), 
                         p3_mut_rel+ theme(legend.position="none"),
                         p4_mut_rel+ theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

# Sig mutation rates: 

means_neutral_without_mut_rel <- out_neutral_mut_rel[[3]]

means_MFED_without_mut_rel <- out_MFED_mut_rel[[3]]

means_CTL_without_mut_rel <- out_CTL_mut_rel[[3]]

means_CTL_MFED_without_mut_rel <- out_CTL_MFED_mut_rel[[3]]

test_0_neutral_without <- vector("numeric", length(M))
test_0_CTL_MFED_without <- vector("numeric", length(M))
test_0_MFED_without <- vector("numeric", length(M))
test_0_CTL_without <- vector("numeric", length(M))

compare = rep(2.16e-4, 10)

for (i in 1:length(m)){

  test_0_neutral_without[i] = binom.test(sum(means_neutral_without_mut_rel[i,]-compare[i]>0),n,0.5)$p.value
  test_0_CTL_MFED_without[i] = binom.test(sum(means_CTL_MFED_without_mut_rel[i,]-compare[i]>0),n,0.5)$p.value
  test_0_MFED_without[i] = binom.test(sum(means_MFED_without_mut_rel[i,]-compare[i]>0),n,0.5)$p.value
  test_0_CTL_without[i] = binom.test(sum(means_CTL_without_mut_rel[i,]-compare[i]>0),n,0.5)$p.value
}

test_0_neutral_without.df <- as.data.frame(as.matrix(x=cbind(m,test_0_neutral_without)))
test_0_MFED_without.df <- as.data.frame(as.matrix(x=cbind(m,test_0_MFED_without)))
test_0_CTL_without.df <- as.data.frame(as.matrix(x=cbind(m,test_0_CTL_without)))
test_0_CTL_MFED_without.df <- as.data.frame(as.matrix(x=cbind(m,test_0_CTL_MFED_without)))

p2_A_mut_rel <- ggplot() + 
  geom_rect(mapping=aes(xmin=0, xmax=0.01, ymin=1e-6, ymax=0.005, fill=TRUE), color=NA, fill='red', alpha=0.2) +
  geom_point(data = test_0_neutral_without.df , aes(x=m, y = test_0_neutral_without), colour="#D1A683") +
  xlab("Simulated migration rate")+ ylab("P-value binomial test for bias") + ggtitle("Neutral") +
  theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank())+
  scale_y_continuous(trans = 'log10', limits = c(1e-6,1)) +
  scale_x_continuous( breaks = seq(0.001,0.01,0.001) ) + 
  scale_colour_manual("", values=c("Using sequence data" ="#D1A683"))

p2_B_mut_rel <- ggplot() + 
  geom_rect(mapping=aes(xmin=0, xmax=0.01, ymin=1e-6, ymax=0.005, fill=TRUE), color=NA, fill='red', alpha=0.2) +
  geom_point(data = test_0_MFED_without.df , aes(x=m, y = test_0_MFED_without), colour="#D1A683") +
  xlab("Simulated migration rate")+ ylab("P-value binomial test for bias") + ggtitle("MFED") +
  theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank())+
  scale_y_continuous(trans = 'log10', limits = c(1e-6,1)) +
  scale_x_continuous( breaks = seq(0.001,0.01,0.001) ) + 
  scale_colour_manual("", values=c("Using True Genealogy" = "#1A2C56","Using sequence data" ="#D1A683"))

p2_C_mut_rel <- ggplot() + 
  geom_rect(mapping=aes(xmin=0, xmax=0.01, ymin=1e-6, ymax=0.005, fill=TRUE), color=NA, fill='red', alpha=0.2) +
  geom_point(data = test_0_CTL_without.df , aes(x=m, y = test_0_CTL_without), colour="#D1A683") +
  xlab("Simulated migration rate")+ ylab("P-value binomial test for bias") + ggtitle("CTL") +
  theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank())+
  scale_y_continuous(trans = 'log10', limits = c(1e-6,1)) +
  scale_x_continuous( breaks = seq(0.001,0.01,0.001) ) + 
  scale_colour_manual("", values=c("Using True Genealogy" = "#1A2C56","Using sequence data" ="#D1A683"))

p2_D_mut_rel <- ggplot() + 
  geom_rect(mapping=aes(xmin=0, xmax=0.01, ymin=1e-6, ymax=0.005, fill=TRUE), color=NA, fill='red', alpha=0.2) +
  geom_point(data = test_0_CTL_MFED_without.df , aes(x=m, y = test_0_CTL_MFED_without), colour="#D1A683") +
  xlab("Simulated migration rate")+ ylab("P-value binomial test for bias") + ggtitle("CTL & MFED") +
  theme(axis.text = element_text(size=18),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom", 
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank())+
  scale_y_continuous(trans = 'log10', limits = c(1e-6,1)) +
  scale_x_continuous( breaks = seq(0.001,0.01,0.001) ) + 
  scale_colour_manual("", values=c("Using True Genealogy" = "#1A2C56","Using sequence data" ="#D1A683"))


mylegend<-g_legend(p2_A_mut_rel)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/bias_detection_p100_mutationrate_relclock.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p2_A_mut_rel + theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_B_mut_rel + theme(legend.position="none", axis.title=element_blank()), 
                         p2_C_mut_rel + theme(legend.position="none"),
                         p2_D_mut_rel + theme(legend.position="none",axis.title.y=element_blank()), nrow = 2), nrow=2, heights=c(10,1))
dev.off()

table_original_mut_rel <- rbind( abs(out_neutral_mut_rel[[3]][,1]-2.16e-4)/2.16e-4*100,
                                  abs(out_MFED_mut_rel[[3]][,1]-2.16e-4)/2.16e-4*100,
                                  abs(out_CTL_mut_rel[[3]][,1]-2.16e-4)/2.16e-4*100,
                                  abs(out_CTL_MFED_mut_rel[[3]][,1]-2.16e-4)/2.16e-4*100)

rownames(table_original_mut_rel) <- c('Neutral without','MFED without', 'CTL without', 'CTL MFED without')
colnames(table_original_mut_rel) <- seq(0.001,0.01,0.001)
