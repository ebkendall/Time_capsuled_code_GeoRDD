set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
load("../../NegControl/Output_tree_rev/match_count_list.dat")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4)
n_matches = 2000

# Theta ------------------------------------------------------------------------
indiv_results_theta = rep(NA, n_buff_width)
global_results_theta = matrix(NA, nrow = n_buff_width, ncol = 2)
colnames(global_results_theta) = c("max", "mean")

# Tau --------------------------------------------------------------------------
indiv_results = rep(NA, length(adjust_val))
global_results = matrix(NA, nrow = length(adjust_val), ncol = 2)
colnames(global_results) = c("max", "mean")

global_null_plot = matrix(nrow = n_matches, ncol = 2)
global_obs_plot = rep(NA, 2)

for(k in 1:n_buff_width) {
    print(paste0("Buffer ", k+2))
    
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
    
    # # Check no zeros for observed boundaries
    # which_zeros_orig = which(origData$str_info$n_off_1_prec == 0 | origData$str_info$n_off_2_prec == 0)
    # if(length(which_zeros_orig) > 0) print(paste0("Zeroes on offenses ", k))
    # 
    # # Remove zeros from the null streets
    # which_zeros = which(combinedMatchingSetupFix2$n_off_1 == 0 | combinedMatchingSetupFix2$n_off_2 == 0)
    # combinedMatchingSetupFix2 = combinedMatchingSetupFix2[-which_zeros, ]
    # int_surface_info = int_surface_info[-which_zeros, ]
    # off_surface_info = off_surface_info[-which_zeros, ]
    
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
    
    # Test statistics: Tau (NEW)
    t_stat_int_surface = int_surface_info
    t_stat_off_surface = off_surface_info
    t_stat_int_surface_orig = origData$str_surf$INT_SURFACE
    t_stat_off_surface_orig = origData$str_surf$OFF_SURFACE
    
    # Theta --------------------------------------------------------------------
    indPvalue_theta = rep(NA, length(indexList_MAIN))
    globalPvalue_theta = rep(NA, 2)
    Y_theta = rep(NA, length(indexList_MAIN))
    X_theta = matrix(NA, length(indexList_MAIN), 2)
    
    counter = 1
    for(ii in indexList_MAIN) {
        
        # Match on crime/offenses
        off_temp = origData$str_info$n_off_1_prec[ii] + origData$str_info$n_off_2_prec[ii]
        ratio_temp = max(origData$str_info$n_off_1_prec[ii] / origData$str_info$n_off_2_prec[ii],
                         origData$str_info$n_off_2_prec[ii] / origData$str_info$n_off_1_prec[ii])
        
        Y_theta[counter] = t_stat_orig[ii]
        X_theta[counter,] = c(ratio_temp, off_temp)
        
        counter = counter + 1
    }
    
    m_theta = match_count_list[["theta"]][k]
    store_theta = NULL
    repeat {
        total_match = rep(0, nrow(X_theta))
        store_theta = matrix(NA, nrow = length(indexList_MAIN), ncol = n_matches)
        for(ii in 1:nrow(X_theta)) {
            
            # Match on crime/offenses
            off_temp = X_theta[ii,2]
            ratio_temp = X_theta[ii,1]
            
            w1 = which(null_sum > m_theta*off_temp & null_sum < (1/m_theta)*off_temp)
            w2 = which(null_ratio >m_theta*ratio_temp & null_ratio < (1/m_theta)*ratio_temp)
            
            w3 = which(is.finite(t_stat))
            
            wAll = intersect(w1, w2)
            wAll = intersect(wAll, w3)
            
            total_match[ii] = length(wAll)
            
            if (length(wAll) > 10) {
                
                testStatsNULL = t_stat[wAll]
                
                testStatsNULL = testStatsNULL[which(testStatsNULL > 0)]
                
                store_theta[ii,] = sample(testStatsNULL, n_matches, replace=TRUE)
            }
        }
        
        # Ensure that at least 75% of the observed boundaries have more than 10 matches
        if(summary(total_match)[2] > 10) {
            break
        } else {
            print("Need more matches (theta)!")
            print(paste0("Current perc = ", m_theta))
            print(summary(total_match))
            m_theta = m_theta - 0.01
        }
    }
    
    for (jj in 1:nrow(X_theta)) {
        indPvalue_theta[jj] = mean(store_theta[jj,] > Y_theta[jj])
    }
    
    globalPvalue_theta[1] = mean(apply(store_theta, 2, max, na.rm=TRUE) > 
                                     max(Y_theta, na.rm=TRUE))
    globalPvalue_theta[2] = mean(apply(store_theta, 2, mean, na.rm=TRUE) > 
                                     mean(Y_theta, na.rm=TRUE))
    
    indiv_results_theta[k] = mean(indPvalue_theta < .05, na.rm=TRUE)
    global_results_theta[k,] = globalPvalue_theta
    
    # Tau ----------------------------------------------------------------------
    # use region of 600 ft for all computations for tau
    if(k == 4) {
        print("Tau")
        for(kk in 1:length(adjust_val)) {
            print(paste0("kk ", kk))
    
            indPvalue = rep(NA, length(indexList_MAIN))
            globalPvalue = rep(NA, 2)
                
            Ytest = rep(NA, length(indexList_MAIN))
            Xtest = matrix(NA, length(indexList_MAIN), 2)
            
            counter = 1
            for(ii in indexList_MAIN) {
                
                # Match on crime/offenses
                off_temp = origData$str_info$n_off_1_prec[ii] + origData$str_info$n_off_2_prec[ii]
                ratio_temp = max(origData$str_info$n_off_1_prec[ii] / origData$str_info$n_off_2_prec[ii],
                                origData$str_info$n_off_2_prec[ii] / origData$str_info$n_off_1_prec[ii])
                
                t_stat_int_surf_orig_kk = abs(t_stat_int_surface_orig[ii,3*kk-2] / t_stat_off_surface_orig[ii,3*kk-2] -
                                                t_stat_int_surface_orig[ii,3*kk-1] / t_stat_off_surface_orig[ii,3*kk-1])
                
                Ytest[counter] = t_stat_int_surf_orig_kk
                Xtest[counter,] = c(ratio_temp, off_temp)
                
                counter = counter + 1
            }
            
            m_tau = match_count_list[["tau"]][kk] # 0.96
            store = NULL
            repeat {
                total_match = rep(0, nrow(Xtest))
                store = matrix(NA, nrow = length(indexList_MAIN), ncol = n_matches)
                for(ii in 1:nrow(Xtest)) {
                    
                    # Match on crime/offenses
                    off_temp = Xtest[ii,2]
                    ratio_temp = Xtest[ii,1]
                    
                    w1 = which(null_sum > m_tau*off_temp & null_sum < (1/m_tau)*off_temp)
                    w2 = which(null_ratio > m_tau*ratio_temp & null_ratio < (1/m_tau)*ratio_temp)
                    
                    null_t_stat_surf = abs(t_stat_int_surface[, 3*kk-2] / t_stat_off_surface[, 3*kk-2]
                                           - t_stat_int_surface[, 3*kk-1] / t_stat_off_surface[, 3*kk-1])
                    
                    w3 = which(is.finite(null_t_stat_surf))
                    
                    wAll = intersect(w1, w2)
                    wAll = intersect(wAll, w3)
                    
                    total_match[ii] = length(wAll)
                    
                    if (length(wAll) > 10) {
                        
                        testStatsNULL = abs(t_stat_int_surface[wAll, 3*kk-2] / t_stat_off_surface[wAll, 3*kk-2]
                                            - t_stat_int_surface[wAll, 3*kk-1] / t_stat_off_surface[wAll, 3*kk-1])
                        
                        testStatsNULL = testStatsNULL[which(testStatsNULL > 0)]
                        
                        store[ii,] = sample(testStatsNULL, n_matches, replace=TRUE)
                    }
                    
                }
                
                # Ensure that at least 75% of the observed boundaries have more than 10 matches
                if(summary(total_match)[2] > 10) {
                    break
                } else {
                    print("Need more matches (tau)!")
                    print(paste0("Current perc = ", m_tau))
                    print(summary(total_match))
                    m_tau = m_tau - 0.01
                }
            }
            
            for (jj in 1:nrow(Xtest)) {
                indPvalue[jj] = mean(store[jj,] > Ytest[jj])
            }
            
            globalPvalue[1] = mean(apply(store, 2, max, na.rm=TRUE) > 
                                        max(Ytest, na.rm=TRUE))
            globalPvalue[2] = mean(apply(store, 2, mean, na.rm=TRUE) > 
                                        mean(Ytest, na.rm=TRUE))
            
            indiv_results[kk] = mean(indPvalue < .05, na.rm=TRUE)
            global_results[kk,] = globalPvalue
            
            if(kk == 2) {
                global_null_plot[,1] = apply(store, 2, max, na.rm=TRUE)
                global_null_plot[,2] = apply(store, 2, mean, na.rm=TRUE)
                global_obs_plot[1] = max(Ytest, na.rm=TRUE)
                global_obs_plot[2] = mean(Ytest, na.rm=TRUE)
            }
        }
    }
}

print("Theta")
print(indiv_results_theta)
print(global_results_theta)
print("Tau")
print(indiv_results)
print(global_results)

# Plotting results ------------------------------------------------------------
library(tidyverse, quietly = T)
library(gridExtra, quietly = T)
library(latex2exp, quietly = T)
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
ggsave(filename = "../Plots_rev/realData_pval.png", 
       plot = realData_naive, width = 1000, height = 800, units = "px")

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
    theme(text = element_text(size=8))
ggsave(filename = "../Plots_rev/realData_pval_total.png", 
       plot = realData_naive_total, width = 1000, height = 800, units = "px")

# Individual test results (theta) ----------------------------------------------
myData_theta <- data.frame(buff_val, indiv_results_theta)
myData_theta$buff_val = myData_theta$buff_val * 100
myData_theta$buff_val = as.factor(myData_theta$buff_val)

i_p_theta = ggplot(myData_theta, aes(y=indiv_results_theta, x=buff_val)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Percent of p-values less than 0.05 (Arrest Data)",
         subtitle=TeX(r'(Estimand: $\theta$)'))+
    xlab(TeX(r'($\delta$)')) +
    ylab("Percent") +
    ylim(0,1)+
    geom_hline(yintercept=0.05, linetype="dashed",
               color = "red", linewidth = 0.5) +
    theme(text = element_text(size=8), legend.position="bottom",
          legend.title=element_blank(),
          legend.key.height= unit(1, 'mm'),
          legend.key.width= unit(4, 'mm'),
          legend.box.margin=margin(-10,-10,-10,-10)) +
    scale_fill_discrete(name = "P-Value")

ggsave(filename = "../Plots_rev/realData_theta_single.png", plot = i_p_theta, 
       width = 1000, height = 800, units = "px")


# Individual test results (tau) ------------------------------------------------
myData_tau <- data.frame(adjust_val, indiv_results)
myData_tau$adjust_val = as.factor(myData_tau$adjust_val)

i_p_tau = ggplot(myData_tau, aes(y=indiv_results, x=adjust_val)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Percent of p-values less than 0.05 (Arrest Data)",
         subtitle=TeX(r'(Estimand: $\tau$)'))+
    xlab("Spatial smoothing multiplier") +
    ylab("Percent") +
    ylim(0,1)+
    geom_hline(yintercept=0.05, linetype="dashed",
               color = "red", linewidth = 0.5) +
    theme(text = element_text(size=8), legend.position="bottom",
          legend.title=element_blank(),
          legend.key.height= unit(1, 'mm'),
          legend.key.width= unit(4, 'mm'),
          legend.box.margin=margin(-10,-10,-10,-10)) +
    scale_fill_discrete(name = "P-Value")

ggsave(filename = "../Plots_rev/realData_tau_single.png", plot = i_p_tau, 
       width = 1000, height = 800, units = "px")

# Global test results (theta) --------------------------------------------------
global_results_plot_theta = data.frame("buff" = as.factor(rep(buff_val*100, 2)),
                                       "p" = c(global_results_theta),
                                       "test_stat" = as.factor(c(rep("max", n_buff_width),
                                                                 rep("mean", n_buff_width))))

global_plot_theta = ggplot(global_results_plot_theta, aes(y=p, x=buff, fill = test_stat)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="P-Values for global test (Arrest Data)")+
    xlab(TeX(r'($\delta$)')) +
    ylab("P-Value") +
    ylim(0,1)+
    geom_hline(yintercept=0.05, linetype="dashed",
               color = "red", linewidth = 0.5) +
    scale_fill_manual(name="Test statistic",
                      labels=c("max", "mean"),
                      values = c("#F8766D", "#00BFC4")) +
    theme(text = element_text(size=8),
          legend.title = element_text(size=5),
          legend.text = element_text(size=5),
          legend.key.size = unit(0.25, 'cm'))

ggsave(filename = "../Plots_rev/realData_theta_global.png", 
       plot = global_plot_theta, width = 1000, height = 800, units = "px")

# Global test results (tau) ----------------------------------------------------
global_results_plot = data.frame("adjust" = as.factor(rep(adjust_val, 2)),
                                 "p" = c(global_results),
                                 "test_stat" = as.factor(c(rep("max", length(adjust_val)),
                                                           rep("mean", length(adjust_val)))))

global_plot_tau = ggplot(global_results_plot, aes(y=p, x=adjust, fill = test_stat)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="P-Values for global test (Arrest Data)")+
    xlab("Spatial smoothing multiplier") +
    ylab("P-Value") +
    ylim(0,1)+
    geom_hline(yintercept=0.05, linetype="dashed",
               color = "red", linewidth = 0.5) +
    scale_fill_manual(name="Test statistic",
                      labels=c("max", "mean"),
                      values = c("#F8766D", "#00BFC4")) +
    theme(text = element_text(size=8),
          legend.title = element_text(size=5),
          legend.text = element_text(size=5),
          legend.key.size = unit(0.25, 'cm'))

ggsave(filename = "../Plots_rev/realData_tau_global.png", 
       plot = global_plot_tau, width = 1000, height = 800, units = "px")


# Histograms for global test -------------------------------------------
globalEmpDist_max = data.frame("num" = log(global_null_plot[,1]))
g_plot_max = ggplot(globalEmpDist_max, aes(x=num)) +
    geom_histogram(color="black", fill="white", bins = floor(sqrt(nrow(globalEmpDist_max)))) +
    geom_vline(xintercept=log(global_obs_plot[1]), color='red', linewidth=0.75)+
    labs(title=TeX(r'(Distribution of log($max_i$ $Z_i$) )'), 
         subtitle=TeX(r'(Spatial smoothing multiplier ($\sigma \times 1$) )')) +
    xlab(TeX(r'(log($max_i$ $Z_i$) )')) +
    ylab("Frequency") +
    theme(legend.position="none", text = element_text(size=6)) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "mm")) +
    theme(axis.title.x = element_text(vjust=2)) +
    theme(plot.title = element_text(vjust=-2))

globalEmpDist_mean = data.frame("num" = log(global_null_plot[,2]))
g_plot_mean = ggplot(globalEmpDist_mean, aes(x=num)) +
    geom_histogram(color="black", fill="white", bins = floor(sqrt(nrow(globalEmpDist_mean)))) +
    geom_vline(xintercept=log(global_obs_plot[2]), color='red', linewidth=0.75)+
    labs(title=TeX(r'(Distribution of log($\bar{Z}$) )'), 
         subtitle=TeX(r'(Spatial smoothing multiplier ($\sigma \times 1$) )')) +
    xlab(TeX(r'(log($\bar{Z}$) )')) +
    ylab("Frequency") +
    theme(legend.position="none", text = element_text(size=6)) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "mm")) +
    theme(axis.title.x = element_text(vjust=2)) +
    theme(plot.title = element_text(vjust=-2))

ggsave(filename = "../Plots_rev/realData_tau_global_hist.png", 
       plot = grid.arrange(g_plot_max, g_plot_mean, nrow = 2), width = 1000, height = 800, units = "px")

