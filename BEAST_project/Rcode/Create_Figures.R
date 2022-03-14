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
rerun = TRUE

#if (!rerun){
#  load(file='')
#}

# mutation rates + number of replicates
M=c('0.001','0.002','0.003', '0.004', '0.005','0.006','0.007', '0.008','0.009','0.01') 

##### Main result figure: migration rate, p = 100, fixed clock #####
if (rerun){
  out_neutral_mig <- create_subFig(selection_type='Neutral', M=M, n = 15)
  out_MFED_mig <- create_subFig(selection_type='MFED', M=M, n = 15)
  out_CTL_mig <- create_subFig(selection_type='CTL', M=M, n = 15)
  out_CTL_MFED_mig <- create_subFig(selection_type='CTL_MFED', M=M, n = 15)
} else {
  load("~/polybox/BEAST_project/Rcode/OriginalModelData.RData")
}

p1 <- out_neutral_mig[[1]]
p2 <- out_MFED_mig[[1]]
p3 <- out_CTL_mig[[1]]
p4 <- out_CTL_MFED_mig[[1]]

mylegend<-g_legend(p1)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/Main_results_p100.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p1+ theme(legend.position="none",axis.title.x=element_blank()), 
                         p2+ theme(legend.position="none", axis.title=element_blank()), 
                         p3+ theme(legend.position="none"),
                         p4+ theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

##### Figure testing the significance of the bias, fixed clock, p =100 #####

means_neutral_with <- out_neutral_mig[[4]]
means_neutral_without <- out_neutral_mig[[5]]

means_MFED_with <- out_MFED_mig[[4]]
means_MFED_without <- out_MFED_mig[[5]]

means_CTL_with <- out_CTL_mig[[4]]
means_CTL_without <- out_CTL_mig[[5]]

means_CTL_MFED_with <- out_CTL_MFED_mig[[4]]
means_CTL_MFED_without <- out_CTL_MFED_mig[[5]]

m=seq(0.001,0.01,0.001)
n = 15

sig_results_mig <- create_sig_fig(means_neutral_with=means_neutral_with, means_neutral_without=means_neutral_without,
                           means_MFED_with=means_MFED_with,means_MFED_without=means_MFED_without,
                           means_CTL_with=means_CTL_with,means_CTL_without=means_CTL_without,
                           means_CTL_MFED_with=means_CTL_MFED_with,means_CTL_MFED_without=means_CTL_MFED_without, m=seq(0.001,0.01,0.001), n =15)

p2_A <- sig_results_mig[[1]]
p2_B <- sig_results_mig[[2]]
p2_C <- sig_results_mig[[3]]
p2_D <- sig_results_mig[[4]]

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/bias_detection_p100.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p2_A + theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_B + theme(legend.position="none", axis.title=element_blank()), 
                         p2_C + theme(legend.position="none"),
                         p2_D + theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

##### Effective population size, p = 100, fixed clock #####

if (rerun){
  out_neutral_tree <- create_subFig(selection_type='Neutral', M=M, measure='tree', n = 15,
                                    xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
  out_MFED_tree <- create_subFig(selection_type='MFED', M=M, measure='tree', n = 15,
                                 xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
  out_CTL_tree <- create_subFig(selection_type='CTL', M=M, measure='tree', n = 15,
                                xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
  out_CTL_MFED_tree <- create_subFig(selection_type='CTL_MFED', M=M, measure='tree', n = 15,
                                     xlabel = "Simulated migration rate", ylabel = "Effective population size", ylimit=c(50,300))
}

p1_tree <- out_neutral_tree[[1]]
p2_tree <- out_MFED_tree[[1]]
p3_tree <- out_CTL_tree[[1]]
p4_tree <- out_CTL_MFED_tree[[1]]

mylegend<-g_legend(p1_tree)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/PopSize_results_p100.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p1_tree+ theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_tree+ theme(legend.position="none", axis.title=element_blank()), 
                         p3_tree+ theme(legend.position="none"),
                         p4_tree+ theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

##### Figure testing the significance of the bias, fixed clock, p =100 #####

means_neutral_with_tree <- out_neutral_tree[[4]]
means_neutral_without_tree <- out_neutral_tree[[5]]

means_MFED_with_tree <- out_MFED_tree[[4]]
means_MFED_without_tree <- out_MFED_tree[[5]]

means_CTL_with_tree <- out_CTL_tree[[4]]
means_CTL_without_tree <- out_CTL_tree[[5]]

means_CTL_MFED_with_tree <- out_CTL_MFED_tree[[4]]
means_CTL_MFED_without_tree <- out_CTL_MFED_tree[[5]]

sig_results_mig_tree <- create_sig_fig(means_neutral_with=means_neutral_with_tree, means_neutral_without=means_neutral_without_tree,
                                  means_MFED_with=means_MFED_with_tree,means_MFED_without=means_MFED_without_tree,
                                  means_CTL_with=means_CTL_with_tree,means_CTL_without=means_CTL_without_tree,
                                  means_CTL_MFED_with=means_CTL_MFED_with_tree,means_CTL_MFED_without=means_CTL_MFED_without_tree, n =15, compare = rep(200, 10))

p2_A_tree <- sig_results_mig_tree[[1]]
p2_B_tree <- sig_results_mig_tree[[2]]
p2_C_tree <- sig_results_mig_tree[[3]]
p2_D_tree <- sig_results_mig_tree[[4]]

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/bias_detection_p100_popsize.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p2_A_tree + theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_B_tree + theme(legend.position="none", axis.title=element_blank()), 
                         p2_C_tree + theme(legend.position="none"),
                         p2_D_tree + theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()


###### Fill the overview result table #######

# Migration rates: 

table_original <- rbind( abs(out_neutral_mig[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_MFED_mig[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                             abs(out_CTL_mig[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                                 abs(out_CTL_MFED_mig[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100)

rownames(table_original) <- c('Neutral without','MFED without', 'CTL without', 'CTL MFED without')
colnames(table_original) <- seq(0.001,0.01,0.001)

table_original_tree <- rbind( abs(out_neutral_tree[[3]][,1]-200)/200*100,
                              abs(out_MFED_tree[[3]][,1]-200)/200*100,
                                  abs(out_CTL_tree[[3]][,1]-200)/200*100,
                                      abs(out_CTL_MFED_tree[[3]][,1]-200)/200*100)

rownames(table_original_tree) <- c('Neutral without','MFED without', 'CTL without', 'CTL MFED without')
colnames(table_original_tree) <- seq(0.001,0.01,0.001)


