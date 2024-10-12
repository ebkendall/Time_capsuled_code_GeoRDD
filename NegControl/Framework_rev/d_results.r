set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
load("../Output_tree_rev/match_count_list.dat")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4)

# Theta ------------------------------------------------------------------------
indiv_results_theta = rep(NA, n_buff_width)
indiv_pvalue_theta = matrix(nrow = n_buff_width, ncol = length(indexList_MAIN))
global_results_theta = matrix(NA, nrow = n_buff_width, ncol = 2)
colnames(global_results_theta) = c("max", "mean")

# Tau --------------------------------------------------------------------------
indiv_results = rep(NA, length(adjust_val))
indiv_pvalue = matrix(nrow = length(adjust_val), ncol = length(indexList_MAIN))
global_results = matrix(NA, nrow = length(adjust_val), ncol = 2)
colnames(global_results) = c("max", "mean")

for(k in 1:n_buff_width) {
    print(paste0("Buffer ", k+2))
    
    load(paste0('../Output_tree_rev/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output_tree_rev/origGridInfo/origData_', k, '.dat'))
    
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
    
    # Test statistics: Theta
    t_stat = abs(combinedMatchingSetupFix2$count1 / combinedMatchingSetupFix2$streets1
                 - combinedMatchingSetupFix2$count2 / combinedMatchingSetupFix2$streets2)
    t_stat_orig = abs(origData$str_info$count1 / origData$str_info$streets1
                      - origData$str_info$count2 / origData$str_info$streets2)
    
    # Test statistics: Tau (NEW)
    t_stat_int_surface = int_surface_info[,3*(1:8)]
    t_stat_int_surface_orig = origData$str_surf$INT_SURFACE[,3*(1:8)]

    rat_off = combinedMatchingSetupFix2$streets1 / combinedMatchingSetupFix2$streets2
    rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
   
    n_matches = 2000
    
    # Theta --------------------------------------------------------------------
    indPvalue_theta = rep(NA, length(indexList_MAIN))
    globalPvalue_theta = rep(NA, 2)
    Y_theta = rep(NA, length(indexList_MAIN))
    X_theta = matrix(NA, length(indexList_MAIN), 2)
    
    counter = 1
    for(ii in indexList_MAIN) {
        
        off_temp = origData$str_info$streets1[ii] + origData$str_info$streets2[ii]
        ratio_temp = max(origData$str_info$streets1[ii] / origData$str_info$streets2[ii],
                         origData$str_info$streets2[ii] / origData$str_info$streets1[ii])
        
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
            
            off_temp = X_theta[ii,2]
            ratio_temp = X_theta[ii,1]
            
            w1 = which((combinedMatchingSetupFix2$streets1 + 
                            combinedMatchingSetupFix2$streets2) > m_theta*off_temp &
                           (combinedMatchingSetupFix2$streets1 + 
                                combinedMatchingSetupFix2$streets2) < (1/m_theta)*off_temp)
            
            w2 = which(rat_off > m_theta*ratio_temp &
                           rat_off < (1/m_theta)*ratio_temp)
            
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
    indiv_pvalue_theta[k, ] = indPvalue_theta
    global_results_theta[k,] = globalPvalue_theta

    # Tau ----------------------------------------------------------------------
    # use region of 600 ft for all computations for tau
    if(k == 4) {
        print("Tau")
        for (kk in 1:length(adjust_val)) {
            print(paste0("kk ", kk))
            
            indPvalue = rep(NA, length(indexList_MAIN))
            globalPvalue = rep(NA, 2)
            
            Ytest = rep(NA, length(indexList_MAIN))
            Xtest = matrix(NA, length(indexList_MAIN), 2)
            
            counter = 1
            for(ii in indexList_MAIN) {
                
                off_temp = origData$str_info$streets1[ii] + origData$str_info$streets2[ii]
                ratio_temp = max(origData$str_info$streets1[ii] / origData$str_info$streets2[ii],
                                 origData$str_info$streets2[ii] / origData$str_info$streets1[ii])
                
                Ytest[counter] = t_stat_int_surface_orig[ii,kk]
                Xtest[counter,] = c(ratio_temp, off_temp)
                
                counter = counter + 1
            }
            
            m_tau = match_count_list[["tau"]][kk]
            store = NULL
            repeat {
                total_match = rep(0, nrow(Xtest))
                store = matrix(NA, nrow = length(indexList_MAIN), ncol = n_matches)
                for(ii in 1:nrow(Xtest)) {
                    
                    off_temp = Xtest[ii,2]
                    ratio_temp = Xtest[ii,1]
                    
                    w1 = which((combinedMatchingSetupFix2$streets1 + 
                                    combinedMatchingSetupFix2$streets2) > m_tau*off_temp &
                                   (combinedMatchingSetupFix2$streets1 + 
                                        combinedMatchingSetupFix2$streets2) < (1/m_tau)*off_temp)
                    
                    w2 = which(rat_off > m_tau*ratio_temp &
                                   rat_off < (1/m_tau)*ratio_temp)
                    w3 = which(is.finite(t_stat_int_surface[, kk]))
                    
                    wAll = intersect(w1, w2)
                    wAll = intersect(wAll, w3)
                    
                    total_match[ii] = length(wAll)
                    
                    if (length(wAll) > 10) {
                        
                        testStatsNULL = t_stat_int_surface[wAll, kk]
                        
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
            indiv_pvalue[kk, ] = indPvalue
            global_results[kk,] = globalPvalue
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
load(paste0('../Output_tree_rev/origGridInfo/origData_4.dat'))
unadjPVal = data.frame("p" = na.omit(origData$str_info$naive_pval))
realData_naive = ggplot(unadjPVal, aes(x=p)) + 
    geom_histogram(color="black", fill="white", bins = floor(sqrt(nrow(origData$str_info)))) +
    xlab("P-Values") + 
    ylab("Frequency") + 
    ggtitle(paste0("Histogram of p-values at buffer width 600 ft")) + 
    theme(text = element_text(size=8))
ggsave(filename = "../Plots_rev/negControl_pval.png", 
       plot = realData_naive, width = 1000, height = 800, units = "px")

load('../Output_tree_rev/origGridInfo/origData_1.dat')
unadjPValTotal = data.frame(na.omit(origData$str_info$naive_pval))
for (i in 2:n_buff_width) {
    load(paste0('../Output_tree_rev/origGridInfo/origData_', i, '.dat'))
    unadjPValTotal = cbind(unadjPValTotal, data.frame(na.omit(origData$str_info$naive_pval)))
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
ggsave(filename = "../Plots_rev/negControl_pval_total.png", 
       plot = realData_naive_total, width = 1000, height = 800, units = "px")

# Individual test results (theta) ----------------------------------------------
myData_theta <- data.frame(buff_val, indiv_results_theta)
myData_theta$buff_val = myData_theta$buff_val * 100
myData_theta$buff_val = as.factor(myData_theta$buff_val)

i_p_theta = ggplot(myData_theta, aes(y=indiv_results_theta, x=buff_val)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Percent of p-values less than 0.05 (Neg. Control)",
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

ggsave(filename = "../Plots_rev/negControl_theta_single.png", plot = i_p_theta, 
       width = 1000, height = 800, units = "px")

# Individual result p-value histograms (theta) ---------------------------------
p = list()
for(i in 1:length(buff_val)) {
    adjPVal500 = data.frame("p" = na.omit(indiv_pvalue_theta[i,]))
    p[[i]] = ggplot(adjPVal500, aes(x=p)) + 
        geom_histogram(color="black", fill="white", bins = sqrt(144)) +
        xlab("P-Values") + 
        ylab("Frequency") + 
        labs(title = "Corrected p-values (Neg. Control)",
             subtitle = substitute(paste("(", delta," = ",m, ")"),list(m=100*(i+2) ))) +
        theme(text = element_text(size=8))
}
pdf("../Plots_rev/negControl_theta_pHist.pdf", onefile = T)
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], nrow = 4, ncol = 2)
dev.off()

# Individual test results (tau) ------------------------------------------------
myData_tau <- data.frame(adjust_val, indiv_results)
myData_tau$adjust_val = as.factor(myData_tau$adjust_val)

i_p_tau = ggplot(myData_tau, aes(y=indiv_results, x=adjust_val)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="Percent of p-values less than 0.05 (Neg. Control)",
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

ggsave(filename = "../Plots_rev/negControl_tau_single.png", plot = i_p_tau, 
       width = 1000, height = 800, units = "px")

# Individual result p-value histograms (tau) ---------------------------------
p_tau = list()
for(i in 1:length(adjust_val)) {
    adjPVal500 = data.frame("p" = na.omit(indiv_pvalue[i,]))
    p_tau[[i]] = ggplot(adjPVal500, aes(x=p)) + 
        geom_histogram(color="black", fill="white", bins = sqrt(144)) +
        xlab("P-Values") + 
        ylab("Frequency") + 
        labs(title = "Corrected p-values (Neg. Control)",
             subtitle = substitute(paste("Spatial smoothing multiplier (", sigma," x ",m, ")"),list(m=adjust_val[i]))) +
        theme(text = element_text(size=8))
}
pdf("../Plots_rev/negControl_tau_pHist.pdf", onefile = T)
grid.arrange(p_tau[[1]], p_tau[[2]], p_tau[[3]], p_tau[[4]], 
             p_tau[[5]], p_tau[[6]], nrow = 3, ncol = 2)
dev.off()

# Global test results (theta) --------------------------------------------------
global_results_plot_theta = data.frame("buff" = as.factor(rep(buff_val*100, 2)),
                                 "p" = c(global_results_theta),
                                 "test_stat" = as.factor(c(rep("max", n_buff_width),
                                                           rep("mean", n_buff_width))))

global_plot_theta = ggplot(global_results_plot_theta, aes(y=p, x=buff, fill = test_stat)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="P-Values for global test (Neg. Control)")+
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

ggsave(filename = "../Plots_rev/negControl_theta_global.png", 
       plot = global_plot_theta, width = 1000, height = 800, units = "px")

# Global test results (tau) ----------------------------------------------------
global_results_plot = data.frame("adjust" = as.factor(rep(adjust_val, 2)),
                                 "p" = c(global_results),
                                 "test_stat" = as.factor(c(rep("max", length(adjust_val)),
                                                           rep("mean", length(adjust_val)))))

global_plot_tau = ggplot(global_results_plot, aes(y=p, x=adjust, fill = test_stat)) +
    geom_bar(position="dodge", stat="identity") +
    labs(title="P-Values for global test (Neg. Control)")+
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

ggsave(filename = "../Plots_rev/negControl_tau_global.png", 
       plot = global_plot_tau, width = 1000, height = 800, units = "px")


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

