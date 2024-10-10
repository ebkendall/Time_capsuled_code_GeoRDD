set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)
load("../../NegControl/Output_tree_rev/match_count_list_old.dat")

perc_rejections_indiv = rep(NA, n_buff_width)
perc_rejections_indiv_surf = matrix(nrow = n_buff_width, ncol = length(adjust_val))
rownames(perc_rejections_indiv_surf) = paste0(seq(from = 300, to = 1000, by = 100), "ft")
colnames(perc_rejections_indiv_surf) = adjust_val

p_global = matrix(nrow = n_buff_width, ncol = 3)
colnames(p_global) = c("max", "mean", "med")
rownames(p_global) = paste0(seq(from = 300, to = 1000, by = 100), "ft")
p_global_surf = vector(mode = 'list', length = 3)
p_global_surf[[1]] = p_global_surf[[2]] = 
    p_global_surf[[3]] = matrix(nrow = n_buff_width, ncol = length(adjust_val))

global_null_plot = matrix(NA, ncol = 2, nrow = 2000)
global_obs_plot = rep(NA, 2)

for(k in 1:n_buff_width) {
    print(paste0("Buffer: ", k+2, " -------------------------"))
    load(paste0('../Output_rev/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output_rev/origGridInfo/origData_', k, '.dat'))
    
    ## Now remove data points where these ratios are much different
    area_ratio = c(na.omit(origData$str_info$area1 / origData$str_info$area2))
    area_ratio[area_ratio < 1] = 1 / area_ratio[area_ratio < 1]
    wMax_a = max(area_ratio)
    wMin_a = min(area_ratio)
    
    street_ratio = c(na.omit(origData$str_info$streets1 / origData$str_info$streets2))
    street_ratio[street_ratio < 1] = 1 / street_ratio[street_ratio < 1]
    wMax_s = max(street_ratio)
    wMin_s = min(street_ratio)
    
    # Which null streets meet these requirements
    wRatioOk = which(combinedMatchingSetupFix$DATA$ratioArea > wMin_a &
                         combinedMatchingSetupFix$DATA$ratioArea < wMax_a & 
                         combinedMatchingSetupFix$DATA$ratioStreet > wMin_s &
                         combinedMatchingSetupFix$DATA$ratioStreet < wMax_s)
    
    combinedMatchingSetupFix2 = combinedMatchingSetupFix$DATA[wRatioOk,]
    int_surface_info = combinedMatchingSetupFix$INT_SURFACE[wRatioOk,]
    off_surface_info = combinedMatchingSetupFix$OFF_SURFACE[wRatioOk,]
    
    # Wherever there is a 0 for the offense count, everything gets scaled by 1
    which_zeros = which(combinedMatchingSetupFix2$n_off_1 == 0 | combinedMatchingSetupFix2$n_off_2 == 0)
    combinedMatchingSetupFix2$n_arr_1[which_zeros] = combinedMatchingSetupFix2$n_arr_1[which_zeros] + 1 
    combinedMatchingSetupFix2$n_arr_2[which_zeros] = combinedMatchingSetupFix2$n_arr_2[which_zeros] + 1 
    combinedMatchingSetupFix2$n_off_1[which_zeros] = combinedMatchingSetupFix2$n_off_1[which_zeros] + 1 
    combinedMatchingSetupFix2$n_off_2[which_zeros] = combinedMatchingSetupFix2$n_off_2[which_zeros] + 1
    
    which_zeros_orig = which(origData$str_info$n_off_1_prec == 0 | origData$str_info$n_off_2_prec == 0)
    origData$str_info$n_arr_1_prec[which_zeros_orig] = origData$str_info$n_arr_1_prec[which_zeros_orig] + 1
    origData$str_info$n_arr_2_prec[which_zeros_orig] = origData$str_info$n_arr_2_prec[which_zeros_orig] + 1
    origData$str_info$n_off_1_prec[which_zeros_orig] = origData$str_info$n_off_1_prec[which_zeros_orig] + 1
    origData$str_info$n_off_2_prec[which_zeros_orig] = origData$str_info$n_off_2_prec[which_zeros_orig] + 1
    
    # Null locations
    null_sum = combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2
    null_ratio = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
    null_ratio[null_ratio < 1] = 1 / null_ratio[null_ratio < 1]
    
    sd1 = sd(null_sum, na.rm = T)
    sd2 = sd(null_ratio, na.rm = T)
    
    # Observed locations
    obs_sum = origData$str_info$n_off_1_prec[indexList_MAIN] + origData$str_info$n_off_2_prec[indexList_MAIN]
    obs_ratio = origData$str_info$n_off_1_prec[indexList_MAIN] / origData$str_info$n_off_2_prec[indexList_MAIN]
    obs_ratio[obs_ratio < 1] = 1 / obs_ratio[obs_ratio < 1]
    
    sd1_obs = sd(obs_sum, na.rm = T)
    sd2_obs = sd(obs_ratio, na.rm = T)
    
    # Test statistics: Theta
    t_stat = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
                 - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)
    t_stat_orig = abs(origData$str_info$n_arr_1_prec / origData$str_info$n_off_1_prec
                      - origData$str_info$n_arr_2_prec / origData$str_info$n_off_2_prec)
    
    # Test statistics: Tau
    t_stat_int_surface = int_surface_info
    t_stat_off_surface = off_surface_info
    t_stat_int_surface_orig = origData$str_surf$INT_SURFACE
    t_stat_off_surface_orig = origData$str_surf$OFF_SURFACE
    
    pval = rep(NA, nrow(origData$str_info))
    pval_int = matrix(NA, nrow = nrow(origData$str_info), ncol = length(adjust_val))
    
    match_vec = match_count_list[k,]
    global_null = matrix(nrow = 164, ncol = match_vec[1])
    global_null_int_surf= vector(mode = 'list', length = length(adjust_val))
    for(jj in 1:length(adjust_val)) {
        global_null_int_surf[[jj]] = matrix(nrow = 164, ncol = match_vec[2])
    }
    
    null_str_position = combinedMatchingSetupFix2[,c("precinct", "indigo", "juliet")]
    matched_null_str_loc = vector(mode = 'list', length = 164)
    
    for(ii in indexList_MAIN) {
        ## find matches
        off_temp = origData$str_info$n_off_1_prec[ii] + origData$str_info$n_off_2_prec[ii]
        ratio_temp = max(origData$str_info$n_off_1_prec[ii] / origData$str_info$n_off_2_prec[ii],
                         origData$str_info$n_off_2_prec[ii] / origData$str_info$n_off_1_prec[ii])
        
        # Remove relatively extreme values to improve the Mahalanobis distance
        remove_extreme1 = which((null_sum > (0.25)*off_temp) & (null_sum < (1/0.25)*off_temp))
        remove_extreme2 = which((null_ratio > (0.25)*ratio_temp) & (null_ratio < (1/0.25)*ratio_temp))
        remove_extreme = intersect(remove_extreme1, remove_extreme2)
        
        null_sum_ii = null_sum[remove_extreme]
        null_ratio_ii = null_ratio[remove_extreme]
        t_stat_ii = t_stat[remove_extreme]
        t_stat_int_surface_ii = t_stat_int_surface[remove_extreme, ]
        t_stat_off_surface_ii = t_stat_off_surface[remove_extreme, ]
        v1_ii = sd(null_sum_ii, na.rm = T)^2
        v2_ii = sd(null_ratio_ii, na.rm = T)^2
        null_str_position_ii = null_str_position[remove_extreme, ]
        
        dist_temp = sqrt(((off_temp - null_sum_ii)^2/v1_ii) + ((ratio_temp - null_ratio_ii)^2 / v2_ii))
        
        w50 = order(dist_temp)[1:match_vec[1]]
        
        null_dist = t_stat_ii[w50]
        global_null[ii,] = null_dist
        
        stat_temp = t_stat_orig[ii]
        pval[ii] = mean(null_dist > stat_temp, na.rm=TRUE)
        
        matched_null_str_loc[[ii]] = null_str_position_ii[w50, ]
        
        for(kk in 1:length(adjust_val)) {
            w50_k = order(dist_temp)[1:match_vec[2]]
            null_dist_int = abs(t_stat_int_surface_ii[w50_k, 3*kk-2] / t_stat_off_surface_ii[w50_k, 3*kk-2]
                                - t_stat_int_surface_ii[w50_k, 3*kk-1] / t_stat_off_surface_ii[w50_k, 3*kk-1])
            global_null_int_surf[[kk]][ii, ] = null_dist_int
            t_stat_int_surf_orig_kk = abs(t_stat_int_surface_orig[ii,3*kk-2] / t_stat_off_surface_orig[ii,3*kk-2] -
                                              t_stat_int_surface_orig[ii,3*kk-1] / t_stat_off_surface_orig[ii,3*kk-1])
            
            pval_int[ii, kk] = mean(null_dist_int > t_stat_int_surf_orig_kk, na.rm=TRUE)
        }
    }
    # Individual test ---------------------------------------------------------
    perc_rejections_indiv[k] = mean(pval < 0.05, na.rm=TRUE)
    
    tau_pval = apply(pval_int, 2, function(x){mean(x < 0.05, na.rm=TRUE)})
    perc_rejections_indiv_surf[k, ] = tau_pval
    
    
    # Global test -------------------------------------------------------------
    n_reps = 2000
    global_t_stat = data.frame("max_t_stat"  = rep(NA, n_reps),
                               "mean_t_stat" = rep(NA, n_reps),
                               "med_t_stat" = rep(NA, n_reps),
                               "max_loc" = rep(NA, n_reps))
    global_t_stat_int_surf = vector(mode='list',length=length(adjust_val))
    for(jj in 1:length(global_t_stat_int_surf)) {
        global_t_stat_int_surf[[jj]] = data.frame("max_t_stat"  = rep(NA, n_reps),
                                                  "mean_t_stat" = rep(NA, n_reps),
                                                  "med_t_stat" = rep(NA, n_reps),
                                                  "max_loc" = rep(NA, n_reps))
    }
    
    # Theta
    for (rep in 1:n_reps) {
        # This is the repetition to get the null distribution
        temp_loc = temp_max = rep(NA, 164)
        for(ii in indexList_MAIN) {
            rand_ind = sample(c(1:match_vec[1]), 1)
            temp_max[ii] = global_null[ii, rand_ind]
            temp_loc[ii] = rand_ind
        }
        global_t_stat[rep, 1] = max(temp_max, na.rm = T)
        global_t_stat[rep, 2] = mean(temp_max, na.rm = T)
        global_t_stat[rep, 3] = median(temp_max, na.rm = T)
        global_t_stat[rep, 4] = temp_loc[which.max(temp_max)]
    }
    
    t_stat_max = max(t_stat_orig, na.rm = T)
    t_stat_mean = mean(t_stat_orig, na.rm = T)
    t_stat_median = median(t_stat_orig, na.rm = T)
    
    p_val_global = rep(NA, 3)
    
    p_val_global[1] = mean(global_t_stat$max_t_stat > t_stat_max)
    p_val_global[2] = mean(global_t_stat$mean_t_stat > t_stat_mean)
    p_val_global[3] = mean(global_t_stat$med_t_stat > t_stat_median)
    
    p_global[k, ] = p_val_global
    
    # Tau
    global_match_loc = vector(mode = 'list', length = length(adjust_val))
    for(av in 1:length(adjust_val)) {
        global_match_loc[[av]] = matrix(0,nrow = n_reps, ncol = ncol(null_str_position))
        
        for (rep in 1:n_reps) {
            # This is the repetition to get the null distribution
            temp_loc = temp_max = rep(NA, 164)
            for(ii in indexList_MAIN) {
                rand_ind = sample(c(1:match_vec[2]), 1)
                temp_max[ii] = global_null_int_surf[[av]][ii, rand_ind]
                temp_loc[ii] = rand_ind
            }
            temp_loc_max = cbind(1:164, temp_loc, temp_max)
            temp_loc_max = temp_loc_max[indexList_MAIN, ]
            temp_loc_max = temp_loc_max[order(abs(temp_loc_max[,3] - median(temp_max, na.rm = T))), ]
            
            med_ii = as.numeric(temp_loc_max[1,1])
            match_ii = as.numeric(temp_loc_max[1,2])
            
            global_match_loc[[av]][rep, ] = as.numeric(matched_null_str_loc[[med_ii]][match_ii, ])
            
            global_t_stat_int_surf[[av]][rep, 1] = max(temp_max, na.rm = T)
            global_t_stat_int_surf[[av]][rep, 2] = mean(temp_max, na.rm = T)
            global_t_stat_int_surf[[av]][rep, 3] = median(temp_max, na.rm = T)
            global_t_stat_int_surf[[av]][rep, 4] = temp_loc[which.max(temp_max)]
        }
    }
    
    if(k == 4) {
        global_null_plot[,1] = global_t_stat_int_surf[[2]][,1]
        global_null_plot[,2] = global_t_stat_int_surf[[2]][,2]
    }
    
    p_val_global_surf = matrix(nrow = 3, ncol = length(adjust_val))
    global_match_loc_orig = rep(NA, length(adjust_val))
    
    for(av in 1:length(adjust_val)) {
        t_stat_scale = abs(t_stat_int_surface_orig[,3*av-2] / t_stat_off_surface_orig[,3*av-2] -
                            t_stat_int_surface_orig[,3*av-1] / t_stat_off_surface_orig[,3*av-1])
        if(av == 2) {
            global_obs_plot[1] = max(t_stat_scale, na.rm = T)
            global_obs_plot[2] = mean(t_stat_scale, na.rm = T)
        }
        for(m in 1:3) {
            if(m == 1) {
                t_stat_kk = max(t_stat_scale, na.rm = T)
                global_null_kk = global_t_stat_int_surf[[av]]$max_t_stat
            } else if(m == 2) {
                t_stat_kk = mean(t_stat_scale, na.rm = T)
                global_null_kk = global_t_stat_int_surf[[av]]$mean_t_stat
            } else {
                t_stat_kk = median(t_stat_scale, na.rm = T)
                temp_orig = cbind(1:164, t_stat_int_surface_orig[,av])
                temp_orig = temp_orig[indexList_MAIN, ]
                temp_orig = temp_orig[order(abs(temp_orig[,2] - t_stat_kk)), ]
                global_match_loc_orig[av] = as.numeric(temp_orig[1,1])
                
                global_null_kk = global_t_stat_int_surf[[av]]$med_t_stat
            }
            
            p_val_global_surf[m, av] = mean(global_null_kk > t_stat_kk, na.rm = T)
        }
    }
    
    colnames(p_val_global_surf) = adjust_val
    rownames(p_val_global_surf) = c("max", "mean", "med")
    
    p_global_surf[[1]][k, ] = p_val_global_surf[1,]
    p_global_surf[[2]][k, ] = p_val_global_surf[2,]
    p_global_surf[[3]][k, ] = p_val_global_surf[3,]
}

print("Individual test results ---------------------------------------------")
print("Buffer & theta & tau")
for(i in 1:8) {
    print(paste0(i+2, "00ft & ", round(perc_rejections_indiv[i], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,1], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,2], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,3], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,4], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,5], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,6], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,7], digits = 4),
                 " & ", round(perc_rejections_indiv_surf[i,8], digits = 4)))
}

print("Global test results -------------------------------------------------")
print("Buffer & theta & tau")
names(p_global_surf) = c("max", "mean", "med")
for(i in 1:8) {
    print(paste0(i+2, "00ft & (", round(p_global[i,1], digits = 3), ", ", round(p_global[i,2], digits = 3),
                 ") & (", round(p_global_surf$max[i,1], digits = 3),", ", round(p_global_surf$mean[i,1], digits = 3),
                 ") & (", round(p_global_surf$max[i,2], digits = 3),", ", round(p_global_surf$mean[i,2], digits = 3),
                 ") & (", round(p_global_surf$max[i,3], digits = 3),", ", round(p_global_surf$mean[i,3], digits = 3),
                 ") & (", round(p_global_surf$max[i,4], digits = 3),", ", round(p_global_surf$mean[i,4], digits = 3),
                 ") & (", round(p_global_surf$max[i,5], digits = 3),", ", round(p_global_surf$mean[i,5], digits = 3),
                 ") & (", round(p_global_surf$max[i,6], digits = 3),", ", round(p_global_surf$mean[i,6], digits = 3),
                 ") & (", round(p_global_surf$max[i,7], digits = 3),", ", round(p_global_surf$mean[i,7], digits = 3),
                 ") & (", round(p_global_surf$max[i,8], digits = 3),", ", round(p_global_surf$mean[i,8], digits = 3), ")"))
}

# Plotting results ------------------------------------------------------------
library(tidyverse, quietly = T)
library(gridExtra, quietly = T)
buff_val = 3:10
n_buff_width = length(buff_val)
# Figure of Naive P-val  -------------------------------------------------------
load(paste0('../Output_rev/origGridInfo/origData_4.dat'))
unadjPVal = data.frame("p" = na.omit(origData$str_info$naive_pval_prec))
realData_naive = ggplot(unadjPVal, aes(x=p)) + 
    geom_histogram(color="black", fill="white", bins = floor(sqrt(nrow(origData$str_info)))) +
    xlab("P-Values") + 
    ylab("Frequency") + 
    ggtitle(paste0("Histogram of p-values at buffer width 600 ft")) + 
    theme(text = element_text(size=8))
ggsave(filename = "../Plots_rev/realData_naive.png", plot = realData_naive, width = 1000, height = 800, units = "px")

load('../Output_rev/origGridInfo/origData_1.dat')
unadjPValTotal = data.frame(na.omit(origData$str_info$naive_pval_prec))
for (i in 2:n_buff_width) {
    load(paste0('../Output_rev/origGridInfo/origData_', i, '.dat'))
    unadjPValTotal = cbind(unadjPValTotal, data.frame(na.omit(origData$str_info$naive_pval_prec)))
}
colnames(unadjPValTotal) = as.character(buff_val)

percRejection = data.frame("perc" = rep(1,length(buff_val)), "buff" = buff_val)
for (i in 1:length(buff_val)) {
    percRejection[i,1] = sum(na.omit(unadjPValTotal[,i] < 0.05)) / sum(!is.na(unadjPValTotal[,i]))
}
percRejection$buff = as.factor(percRejection$buff)

print(percRejection)

realData_naive_total = ggplot(percRejection, aes(y=perc, x=buff)) + 
    geom_bar(position="dodge", stat="identity") + 
    ggtitle("Percentage of p-values less than 0.05") +
    xlab("Buffer width (100x in ft)") + 
    ylab("Percent") +
    ylim(0, 1) +
    # scale_x_continuous(breaks = pretty(percRejection$buff, n = 8)) +
    theme(text = element_text(size=8))
ggsave(filename = "../Plots_rev/realData_naive_total.png", plot = realData_naive_total, width = 1000, height = 800, units = "px")

# Individual test results ---------------------------------------------------
i_p = vector(mode = 'list', length = length(buff_val))
for(b in 1:length(buff_val)) {
    yVal = perc_rejections_indiv_surf[b,]
    myData <- data.frame(adjust_val, yVal)
    
    myData$adjust_val = as.factor(myData$adjust_val)
    
    i_p[[b]] = ggplot(myData, aes(y=yVal, x=adjust_val)) +
        geom_bar(position="dodge", stat="identity") +
        labs(title="Percent of p-values less than 0.05 (Arrest Data)",
             subtitle=paste0("Intensity surface correction, region size = ", b+2, "00 ft"))+
        xlab("Spatial smoothing multiplier") +
        ylab("Percent") +
        ylim(0,1)+
        # scale_x_continuous(breaks=1:length(adjust_val)) +
        # scale_y_continuous(breaks = seq(0.1,1,by=0.1)) +
        geom_hline(yintercept=0.05, linetype="dashed",
                   color = "red", linewidth = 0.5) +
        theme(text = element_text(size=8), legend.position="bottom",
              legend.title=element_blank(),
              legend.key.height= unit(1, 'mm'),
              legend.key.width= unit(4, 'mm'),
              legend.box.margin=margin(-10,-10,-10,-10)) +
        scale_fill_discrete(name = "P-Value")
}

pdf("../Plots_rev/int_surface_results.pdf", onefile = T)
grid.arrange(i_p[[1]], i_p[[2]], i_p[[3]], i_p[[4]], ncol = 2, nrow = 2)
grid.arrange(i_p[[5]], i_p[[6]], i_p[[7]], i_p[[8]], ncol = 2, nrow = 2)
dev.off()

ggsave(filename = "../Plots_rev/int_surface_single_real.png", plot = i_p[[4]], 
       width = 1000, height = 800, units = "px")

# Global test results ---------------------------------------------------
global_results_plot = data.frame("adjust" = as.factor(rep(adjust_val, 2)),
                                 "p" = c(p_global_surf$max[i,], p_global_surf$mean[i,]),
                                 "test_stat" = as.factor(c(rep("max", length(adjust_val)),
                                                           rep("mean", length(adjust_val)))))

realData_global_total = ggplot(global_results_plot, aes(y=p, x=adjust, fill = test_stat)) +
    geom_bar(position="dodge", stat="identity") +
    ggtitle("P-Values for global test") +
    xlab("Spatial smoothing multiplier") +
    ylab("P-Value") +
    geom_hline(yintercept=0.05, linetype="dashed",
               color = "red", linewidth = 0.5) +
    scale_fill_manual(name="Test statistic",
                        labels=c("max", "mean"),
                        values = c("#F8766D", "#00BFC4")) +
    theme(text = element_text(size=8),
          legend.title = element_text(size=5),
          legend.text = element_text(size=5),
          legend.key.size = unit(0.25, 'cm'))

ggsave(filename = "../Plots_rev/realData_global_total_surf.png", plot = realData_global_total, width = 1000, height = 800, units = "px")


# # Histograms for global test -------------------------------------------
# globalEmpDist = data.frame("num" = global_null_plot)
# g_plot_max = ggplot(globalEmpDist, aes(x=num)) +
#     geom_histogram(color="black", fill="white", bins = floor(sqrt(nrow(globalEmpDist)))) +
#     geom_vline(data=globalObsVal, aes(xintercept=obs, color="red"), size=1) +
#     ggtitle(paste0("Distr. of global test stat. (smoothing mult. = ", adjust_val[kk], ")")) +
#     xlab("Test statistic") +
#     ylab("Frequency") +
#     theme(legend.position="none", text = element_text(size=6)) +
#     theme(plot.margin = unit(c(.2,.2,.2,.2), "mm")) +
#     theme(axis.title.x = element_text(vjust=2)) +
#     theme(plot.title = element_text(vjust=-2))
# 
# realData_global_sep_surf = grid.arrange(g_plots_surf[[2]][[2]], g_plots_surf[[2]][[4]], g_plots_surf[[2]][[5]], nrow = 3)
# ggsave(filename = "Plots/realData_global_sep_surf.png", plot = realData_global_sep_surf, width = 1000, height = 800, units = "px")
