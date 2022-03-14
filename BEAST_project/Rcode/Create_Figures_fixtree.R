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

# mutation rates + number of replicates
M=c('0.001','0.002','0.003', '0.004', '0.005','0.006','0.007', '0.008','0.009','0.01') 

##### Result figure: migration rate, p = 100, fixed clock, and fixed pop #####
if (rerun){
  out_neutral_mig_fixpop <- create_subFig(selection_type='Neutral', M=M, n = 15, folder = 'P100_forPaper_fixpop',
                                          data = 'fixpop_original_2')
  out_MFED_mig_fixpop <- create_subFig(selection_type='MFED', M=M, n = 15, folder = 'P100_forPaper_fixpop',
                                       data = 'fixpop_original_2')
  out_CTL_mig_fixpop <- create_subFig(selection_type='CTL', M=M, n = 15, folder = 'P100_forPaper_fixpop',
                                      data = 'fixpop_original_2')
  out_CTL_MFED_mig_fixpop <- create_subFig(selection_type='CTL_MFED', M=M, n = 15, folder = 'P100_forPaper_fixpop',
                                           data = 'fixpop_original_2')
} else {
  load("~/polybox/BEAST_project/Rcode/Data_fixedPop.RData")
}

p1_fixpop <- out_neutral_mig_fixpop[[1]]
p2_fixpop <- out_MFED_mig_fixpop[[1]]
p3_fixpop <- out_CTL_mig_fixpop[[1]]
p4_fixpop <- out_CTL_MFED_mig_fixpop[[1]]

mylegend<-g_legend(p1_fixpop)

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/Main_results_p100_fixpop.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p1_fixpop+ theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_fixpop+ theme(legend.position="none", axis.title=element_blank()), 
                         p3_fixpop+ theme(legend.position="none"),
                         p4_fixpop+ theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

#### Significance 

means_neutral_with_fixpop <- out_neutral_mig_fixpop[[4]]
means_neutral_without_fixpop <- out_neutral_mig_fixpop[[5]]

means_MFED_with_fixpop <- out_MFED_mig_fixpop[[4]]
means_MFED_without_fixpop <- out_MFED_mig_fixpop[[5]]

means_CTL_with_fixpop <- out_CTL_mig_fixpop[[4]]
means_CTL_without_fixpop <- out_CTL_mig_fixpop[[5]]

means_CTL_MFED_with_fixpop <- out_CTL_MFED_mig_fixpop[[4]]
means_CTL_MFED_without_fixpop <- out_CTL_MFED_mig_fixpop[[5]]

m=seq(0.001,0.01,0.001)
n = 15

sig_results_mig_fixpop <- create_sig_fig(means_neutral_with=means_neutral_with_fixpop, means_neutral_without=means_neutral_without_fixpop,
                                  means_MFED_with=means_MFED_with_fixpop,means_MFED_without=means_MFED_without_fixpop,
                                  means_CTL_with=means_CTL_with_fixpop,means_CTL_without=means_CTL_without_fixpop,
                                  means_CTL_MFED_with=means_CTL_MFED_with_fixpop,means_CTL_MFED_without=means_CTL_MFED_without_fixpop, m=seq(0.001,0.01,0.001), n =15)

p2_A_fixpop <- sig_results_mig_fixpop[[1]]
p2_B_fixpop <- sig_results_mig_fixpop[[2]]
p2_C_fixpop <- sig_results_mig_fixpop[[3]]
p2_D_fixpop <- sig_results_mig_fixpop[[4]]

pdf(file = "/Users/judith/polybox/BEAST_project/Figures/bias_detection_p100_fixpop.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 12)
grid.arrange(arrangeGrob(p2_A_fixpop + theme(legend.position="none",axis.title.x=element_blank()), 
                         p2_B_fixpop + theme(legend.position="none", axis.title=element_blank()), 
                         p2_C_fixpop + theme(legend.position="none"),
                         p2_D_fixpop + theme(legend.position="none",axis.title.y=element_blank()), nrow = 2),mylegend, nrow=2, heights=c(10,1))
dev.off()

table_fixpop <- rbind( abs(out_neutral_mig_fixpop[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_MFED_mig_fixpop[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_CTL_mig_fixpop[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100,
                         abs(out_CTL_MFED_mig_fixpop[[3]][,1]-seq(0.001,0.01,0.001))/seq(0.001,0.01,0.001)*100)

rownames(table_fixpop) <- c('Neutral without','MFED without', 'CTL without', 'CTL MFED without')
colnames(table_fixpop) <- seq(0.001,0.01,0.001)


