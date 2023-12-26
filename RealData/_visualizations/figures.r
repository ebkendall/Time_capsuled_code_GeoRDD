# ALL FIGURES ARE 1000 x 800

library(tidyverse, quietly = T)
library(gridExtra)

adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)

# Figure of Naive P-val  -------------------------------------------------------
load('../Output/origGridInfo/sim_orig_5.dat')
unadjPVal500 = data.frame("p" = na.omit(sim_orig$DATA$naive_pval))
realData_naive = ggplot(unadjPVal500, aes(x=p)) + 
                          geom_histogram(color="black", fill="white", bins = floor(sqrt(nrow(sim_orig$DATA)))) +
                          xlab("P-Values") + 
                          ylab("Frequency") + 
                          ggtitle("Histogram of p-values at buffer width 500 ft") + 
                          theme(text = element_text(size=8))
ggsave(filename = "Plots/realData_naive.png", plot = realData_naive, width = 1000, height = 800, units = "px")

load('../Output/origGridInfo/sim_orig_3.dat')
unadjPValTotal = data.frame(na.omit(sim_orig$DATA$naive_pval))
for (i in 1:length(adjust_val)) {
  load(paste0('../Output/origGridInfo/sim_orig_', i, '.dat'))
  unadjPValTotal = cbind(unadjPValTotal, data.frame(na.omit(sim_orig$DATA$naive_pval)))
}
colnames(unadjPValTotal) = as.character(adjust_val)

percRejection = data.frame("perc" = rep(1,length(adjust_val)), "buff" = adjust_val)
for (i in 1:length(adjust_val)) {
  percRejection[i,1] = sum(na.omit(unadjPValTotal[,i] < 0.05)) / sum(!is.na(unadjPValTotal[,i]))
}
percRejection$buff = as.factor(percRejection$buff)

realData_naive_total = ggplot(percRejection, aes(y=perc, x=buff)) + 
                          geom_bar(position="dodge", stat="identity") + 
                          ggtitle("Percentage of p-values less than 0.05") +
                          xlab("Buffer Width (100x in ft)") + 
                          ylab("Percent") +
                          ylim(0, 1) +
                          # scale_x_continuous(breaks = pretty(percRejection$buff, n = 8)) +
                          theme(text = element_text(size=8))
ggsave(filename = "Plots/realData_naive_total.png", plot = realData_naive_total, width = 1000, height = 800, units = "px")

# Figure of Individual results -----------------------------------------------------
perc_rejections_new = data.frame("orig" = c(1:length(adjust_val)), "adjusted" = c(1:length(adjust_val)))
load("../Output/p_vals_match_rel/p_val_df_new_stat_FINAL.dat")
for (i in 1:length(adjust_val)) {
    load(paste0('../Output/origGridInfo/sim_orig_', i, '.dat'))
    p_orig = mean(sim_orig$DATA$naive_pval < 0.05, na.rm = T)
    p_new = mean(p_val_df[[i]][1,] < 0.05, na.rm = T)
    perc_rejections_new[i, ] = c(p_orig, p_new)
}
# perc_rejections_new = perc_rejections_new[3:10,]

buff = rep(adjust_val, 2)
pValtype = c(rep("Naive", length(adjust_val)), rep("Corrected", length(adjust_val)))
yVal = c(perc_rejections_new[,1], perc_rejections_new[,2])
myData <- data.frame(buff,pValtype, yVal)

myData$buff = as.factor(myData$buff)

realData_pval_final = ggplot(myData, aes(fill=pValtype, y=yVal, x=buff)) + 
    geom_bar(position="dodge", stat="identity") + 
    ggtitle("Percent of p-values less than 0.05 (Arrest Data)") +
    xlab("Buffer Width (100x in ft)") + 
    ylab("Percent") +
    # scale_x_continuous(breaks = 3:10) +
    scale_y_continuous(breaks = pretty(myData$yVal, n = 10)) +
    geom_hline(yintercept=0.05, linetype="dashed", 
               color = "black", size = 0.5) +
    theme(text = element_text(size=8), legend.position="bottom",
          legend.title=element_blank(),
          legend.key.height= unit(1, 'mm'),
          legend.key.width= unit(4, 'mm'),
          legend.box.margin=margin(-10,-10,-10,-10)) +
    scale_fill_discrete(name = "P-Value")

ggsave(filename = "Plots/realData_pval_final.png", plot = realData_pval_final, width = 1000, height = 800, units = "px")

# Updated P-Value histograms ---------------------------------------------------
load("../Output/p_vals_match_rel/p_val_df_new_stat_FINAL.dat")
p = vector(mode = 'list', length = length(adjust_val))
for(i in 1:length(adjust_val)) {
    adjPVal500 = data.frame("p" = na.omit(p_val_df[[i]][1,]))
    p[[i]] = ggplot(adjPVal500, aes(x=p)) + 
        geom_histogram(color="black", fill="white", bins = sqrt(144)) +
        xlab("P-Values") + 
        ylab("Frequency") + 
        ggtitle(paste0("Corrected p-values at buffer ",i,"00 ft")) + 
        theme(text = element_text(size=8))
}
pdf("Plots/correctedHist.pdf", onefile = T)
grid.arrange(p[[1]], p[[2]], ncol = 1, nrow = 2)
grid.arrange(p[[3]], p[[4]], ncol = 1, nrow = 2)
grid.arrange(p[[5]], p[[6]], ncol = 1, nrow = 2)
grid.arrange(p[[7]], p[[8]], ncol = 1, nrow = 2)
grid.arrange(p[[9]], p[[10]], ncol = 1, nrow = 2)
grid.arrange(p[[11]], p[[12]], ncol = 1, nrow = 2)
dev.off()

df_new = data.frame("p" = na.omit(p_val_df[[5]][1, ]))
realData_pval_500 = ggplot(df_new, aes(x=p)) + 
    geom_histogram(color="black", fill="white", bins = sqrt(144)) +
    xlab("P-Values") + 
    ylab("Frequency") + 
    ggtitle(paste0("Corrected p-values at buffer 500 ft")) + 
    theme(text = element_text(size=8))
ggsave(filename = "Plots/realData_pval_500.png", plot = realData_pval_500, width = 1000, height = 800, units = "px")

# Figure of Global results -----------------------------------------------------
# load("../Data/indexList_MAIN.RData")
# 
# p_val_df = rep(NA, 10)
# 
# load('../Output/Global/global_t_stat_FINAL.dat')
# 
# for(k in 3:10) {
#     
#     load(paste0('../Output/origGridInfo/sim_orig_', k,".dat"))
#     which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
#     sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
#     sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
#     sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
#     sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1
# 
#     t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
#                       - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
#     
#     t_stat = max(t_stat_orig, na.rm = T)
#     
#     p_val_df[k] = mean(global_t_stat[[k]]$max_t_stat > t_stat)
#     
# }
# print("Global p-values for each buffer width")
# print(p_val_df[3:10])
# 
# globalPvals_new_2 = data.frame("buff" = 3:10,
#                                "p" = p_val_df[3:10])
# realData_global_total = ggplot(globalPvals_new_2, aes(y=p, x=buff)) + 
#                       geom_bar(position="dodge", stat="identity") + 
#                       ggtitle("P-Values for global test") +
#                       xlab("Buffer Width (100x in ft)") + 
#                       ylab("P-Values") +
#                       geom_hline(yintercept=0.05, linetype="dashed", 
#                                  color = "red", size = 0.5) +
#                       theme(text = element_text(size=8)) +
#                       scale_x_continuous(breaks= 2:10)
# 
# ggsave(filename = "Plots/realData_global_total.png", plot = realData_global_total, width = 1000, height = 800, units = "px")
# 
# ind = 1
# g_plots = vector(mode = "list", length = 3)
# for(k in c(3,6,9)) {
#     load(paste0('../Output/origGridInfo/sim_orig_', k,".dat"))
#     which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
#     sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
#     sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
#     sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
#     sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1
# 
#     t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
#                       - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
#     
#     t_stat = max(t_stat_orig, na.rm = T)
# 
#     globalEmpDist = data.frame("num" = global_t_stat[[k]]$max_t_stat)
# 
#     globalObsVal = data.frame("obs" = t_stat)
# 
#     g_plots[[ind]] = ggplot(globalEmpDist, aes(x=num)) + 
#                 geom_histogram(color="black", fill="white", bins = sqrt(nrow(globalEmpDist))) + 
#                 geom_vline(data=globalObsVal, aes(xintercept=obs, color="red"), size=1) + 
#                 ggtitle(paste0("Distribution of global test statistic (", k, "00 ft)")) +
#                 xlab("Test Statistic") + 
#                 ylab("Frequency") +
#                 theme(legend.position="none", text = element_text(size=6)) +
#                 theme(plot.margin = unit(c(.2,.2,.2,.2), "mm")) +
#                 theme(axis.title.x = element_text(vjust=2)) +
#                 theme(plot.title = element_text(vjust=-3)) 
#   
#     ind=ind+1
# }
# 
# realData_global_sep = grid.arrange(g_plots[[1]], g_plots[[2]], g_plots[[3]], nrow = 3)
# ggsave(filename = "Plots/realData_global_sep.png", plot = realData_global_sep, width = 1000, height = 800, units = "px")
# 
# # Newest plots to show importance of matching ----------------------------------
# library(splines)
# library(tidyverse)
# k = 5
# load(paste0('../Output/nullGridInfo/combinedMatchingSetup', k, ".dat"))
# load(paste0('../Output/origGridInfo/sim_orig_', k, '.dat'))
# 
# ## Now remove data points where these ratios are much different
# area_ratio = c(na.omit(sim_orig$DATA$area1 / sim_orig$DATA$area2))
# area_ratio[area_ratio < 1] = 1 / area_ratio[area_ratio < 1]
# wMax_a = max(area_ratio)
# wMin_a = min(area_ratio)
# 
# street_ratio = c(na.omit(sim_orig$DATA$streets1 / sim_orig$DATA$streets2))
# street_ratio[street_ratio < 1] = 1 / street_ratio[street_ratio < 1]
# wMax_s = max(street_ratio)
# wMin_s = min(street_ratio)
# 
# wRatioOk = which(combinedMatchingSetupFix$DATA$ratioArea > wMin_a &
#                      combinedMatchingSetupFix$DATA$ratioArea < wMax_a & 
#                      combinedMatchingSetupFix$DATA$ratioStreet > wMin_s &
#                      combinedMatchingSetupFix$DATA$ratioStreet < wMax_s)
# 
# combinedMatchingSetupFix2 = combinedMatchingSetupFix$DATA[wRatioOk,]
# 
# # Wherever there is a 0 for the offense count, everything gets scaled by 1
# which_zeros = which(combinedMatchingSetupFix2$n_off_1 == 0 | combinedMatchingSetupFix2$n_off_2 == 0)
# combinedMatchingSetupFix2$n_arr_1[which_zeros] = combinedMatchingSetupFix2$n_arr_1[which_zeros] + 1 
# combinedMatchingSetupFix2$n_arr_2[which_zeros] = combinedMatchingSetupFix2$n_arr_2[which_zeros] + 1 
# combinedMatchingSetupFix2$n_off_1[which_zeros] = combinedMatchingSetupFix2$n_off_1[which_zeros] + 1 
# combinedMatchingSetupFix2$n_off_2[which_zeros] = combinedMatchingSetupFix2$n_off_2[which_zeros] + 1
# 
# which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
# sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
# sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
# sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
# sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1
# 
# rat_off = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
# rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
# tot_off = combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2
# t_stat_new = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
#              - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)
# 
# plotting_info = data.frame("t_stat_new" = t_stat_new, "totLength" = tot_off, "ratioOff" = rat_off)
# plotting_info = plotting_info[!is.na(plotting_info$totLength), ]
# 
# ratio_breaks = seq(0, 60, 6)
# sum_breaks   = seq(1200, 0, -60)
# t_stat_plot = matrix(nrow = length(ratio_breaks)-1, ncol = length(sum_breaks)-1)
# var_stat_plot = matrix(nrow = length(ratio_breaks)-1, ncol = length(sum_breaks)-1)
# for(i in 2:length(ratio_breaks)) {
#     sub_rat = plotting_info[plotting_info$ratioOff <= ratio_breaks[i] &
#                                 plotting_info$ratioOff > ratio_breaks[i-1], , drop = F]
#     for(j in 2:length(sum_breaks)) {
#         sub_tot = sub_rat[sub_rat$totLength > sum_breaks[j] & sub_rat$totLength <= sum_breaks[j-1], ,drop=F]
#         if(nrow(sub_tot) > 0) {
#             t_stat_plot[i-1, j-1] = mean(sub_tot$t_stat_new, na.rm = T)
#             var_stat_plot[i-1, j-1] = var(sub_tot$t_stat_new, na.rm = T)
#         }
#     }
# }
# c_name = NULL
# r_name = NULL
# for(i in 2:length(ratio_breaks)) r_name = c(r_name, paste0(ratio_breaks[i-1], "-", ratio_breaks[i]))
# for(i in 2:length(sum_breaks)) c_name = c(c_name, paste0(sum_breaks[i], "-", sum_breaks[i-1]))
# 
# colnames(t_stat_plot) = c_name
# rownames(t_stat_plot) = r_name
# 
# colnames(var_stat_plot) = c_name
# rownames(var_stat_plot) = r_name
# 
# png(filename = "Plots/appendix_res.png", width = 1000, height = 800,
#     units = "px", pointsize = 12, bg = "white", res = NA)
# heatmap(t_stat_plot, Colv = NA, Rowv = NA, cexRow = 3, cexCol = 3, margins = c(12,2))
# dev.off()
# 
# png(filename = "Plots/appendix_res2.png", width = 1000, height = 800,
#     units = "px", pointsize = 12, bg = "white", res = NA)
# heatmap(var_stat_plot, Colv = NA, Rowv = NA, cexRow = 3, cexCol = 3, margins = c(12,2))
# dev.off()

