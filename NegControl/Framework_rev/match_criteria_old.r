match_count <- seq(20, 1200, by = 20)
set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

indiv_results_theta = indiv_results_tau = matrix(NA, nrow = n_buff_width, 
                                                 ncol = length(match_count))

for(k in 1:n_buff_width) {
    print(paste0("buff ", k + 2))
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
    
    # Null locations
    null_sum = combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2
    null_ratio = combinedMatchingSetupFix2$streets1 / combinedMatchingSetupFix2$streets2
    null_ratio[null_ratio < 1] = 1 / null_ratio[null_ratio < 1]
    
    sd1 = sd(null_sum, na.rm = T)
    sd2 = sd(null_ratio, na.rm = T)
    
    # Observed locations
    obs_sum = origData$str_info$streets1[indexList_MAIN] + origData$str_info$streets2[indexList_MAIN]
    obs_ratio = origData$str_info$streets1[indexList_MAIN] / origData$str_info$streets2[indexList_MAIN]
    obs_ratio[obs_ratio < 1] = 1 / obs_ratio[obs_ratio < 1]
    
    sd1_obs = sd(obs_sum, na.rm = T)
    sd2_obs = sd(obs_ratio, na.rm = T)
    
    # Test statistics: Theta
    t_stat = abs(combinedMatchingSetupFix2$count1 / combinedMatchingSetupFix2$streets1
                 - combinedMatchingSetupFix2$count2 / combinedMatchingSetupFix2$streets2)
    t_stat_orig = abs(origData$str_info$count1 / origData$str_info$streets1
                      - origData$str_info$count2 / origData$str_info$streets2)
    
    # Test statistics: Tau (NEW)
    t_stat_int_surface = int_surface_info[,3*(1:8)]
    t_stat_int_surface_orig = origData$str_surf$INT_SURFACE[,3*(1:8)]
    
    pval = matrix(NA, nrow = nrow(origData$str_info), ncol = length(match_count))
    pval_int = matrix(NA, nrow = nrow(origData$str_info), ncol = length(match_count))
    kk = 2
    
    for(ii in indexList_MAIN) {
        print(ii)
        ## find matches
        off_temp = origData$str_info$streets1[ii] + origData$str_info$streets2[ii]
        ratio_temp = max(origData$str_info$streets1[ii] / origData$str_info$streets2[ii],
                         origData$str_info$streets2[ii] / origData$str_info$streets1[ii])
        
        # Remove relatively extreme values to improve the mahalanobis distance
        remove_extreme1 = which((null_sum > (0.25)*off_temp) & (null_sum < (1/0.25)*off_temp))
        remove_extreme2 = which((null_ratio > (0.25)*ratio_temp) & (null_ratio < (1/0.25)*ratio_temp))
        remove_extreme = intersect(remove_extreme1, remove_extreme2)
        
        null_sum_ii = null_sum[remove_extreme]
        null_ratio_ii = null_ratio[remove_extreme]
        t_stat_ii = t_stat[remove_extreme]
        t_stat_int_surface_ii = t_stat_int_surface[remove_extreme, kk]
        v1_ii = sd(null_sum_ii, na.rm = T)^2
        v2_ii = sd(null_ratio_ii, na.rm = T)^2
        
        dist_temp = sqrt(((off_temp - null_sum_ii)^2/v1_ii) + ((ratio_temp - null_ratio_ii)^2 / v2_ii))
        stat_temp = t_stat_orig[ii]
        stat_temp_int = t_stat_int_surface_orig[ii,kk]
        
        if(length(dist_temp) < 20) next
        
        for(m in 1:length(match_count)) {
            if(length(dist_temp) < match_count[m]) {
                w50 = order(dist_temp)
            } else {
                w50 = order(dist_temp)[1:match_count[m]]   
            }
            
            null_dist = t_stat_ii[w50]
            pval[ii, m] = mean(null_dist > stat_temp, na.rm=TRUE)
            
            null_dist_int = t_stat_int_surface_ii[w50]
            pval_int[ii, m] = mean(null_dist_int > stat_temp_int, na.rm=TRUE)
        }
    }
    
    indiv_results_theta[k, ] = apply(pval, 2, function(x){mean(x < 0.05, na.rm=TRUE)})
    indiv_results_tau[k, ] = apply(pval_int, 2, function(x){mean(x < 0.05, na.rm=TRUE)})
}

save(indiv_results_theta, file = '../Output_tree_rev/indiv_results_theta_old.rda')
save(indiv_results_tau, file = '../Output_tree_rev/indiv_results_tau_old.rda')


# Plotting the results
library(tidyverse, quietly = T)
library(gridExtra, quietly = T)
p_theta = p_tau = vector(mode = 'list', length = 8)
load('../Output_tree_rev/indiv_results_theta_old.rda')
load('../Output_tree_rev/indiv_results_tau_old.rda')
for(i in 1:8) {
    p_val_theta = data.frame(y = indiv_results_theta[i,],
                             x = match_count)
    p_theta[[i]] = ggplot(p_val_theta, aes( y=y, x=x)) +
        geom_point(color = "red", size = 2) +
        geom_smooth(method = "loess", formula = y ~ x, span=0.5) +
        ggtitle(paste0("Theta: Matching's Effect on Type I Error (",i+2,"00 ft)")) +
        xlab("Measure of match similarity") +
        ylab("Type I Error") +
        ylim(0,max(p_val_theta$y)) +
        scale_x_continuous(breaks = pretty(p_val_theta$x, n = 10)) +
        geom_hline(yintercept=0.05, linetype="dashed",
                   color = "black", size = 1.5) +
        theme(text = element_text(size=15))
    
    p_val_tau = data.frame(y = indiv_results_tau[i,],
                           x = match_count)
    p_tau[[i]] = ggplot(p_val_tau, aes( y=y, x=x)) +
        geom_point(color = "red", size = 2) +
        geom_smooth(method = "loess", formula = y ~ x, span=0.5) +
        ggtitle(paste0("Tau: Matching's Effect on Type I Error (",i+2,"00 ft)")) +
        xlab("Measure of match similarity") +
        ylab("Type I Error") +
        ylim(0,max(p_val_tau$y)) +
        scale_x_continuous(breaks = pretty(p_val_tau$x, n = 10)) +
        geom_hline(yintercept=0.05, linetype="dashed",
                   color = "black", size = 1.5) +
        theme(text = element_text(size=15))
}
pdf("../Plots_rev/match_theta_old.pdf", onefile = T)
grid.arrange(p_theta[[1]], p_theta[[2]], ncol = 1, nrow = 2)
grid.arrange(p_theta[[3]], p_theta[[4]], ncol = 1, nrow = 2)
grid.arrange(p_theta[[5]], p_theta[[6]], ncol = 1, nrow = 2)
grid.arrange(p_theta[[7]], p_theta[[8]], ncol = 1, nrow = 2)
dev.off()

pdf("../Plots_rev/match_tau_old.pdf", onefile = T)
grid.arrange(p_tau[[1]], p_tau[[2]], ncol = 1, nrow = 2)
grid.arrange(p_tau[[3]], p_tau[[4]], ncol = 1, nrow = 2)
grid.arrange(p_tau[[5]], p_tau[[6]], ncol = 1, nrow = 2)
grid.arrange(p_tau[[7]], p_tau[[8]], ncol = 1, nrow = 2)
dev.off()



match_count_list = matrix(nrow = n_buff_width, ncol = 2)
colnames(match_count_list) = c('theta', 'tau')
for(i in 1:8) {
    perc_rej_theta = indiv_results_theta[i,]
    perc_rej_tau = indiv_results_tau[i,]

    match_count_list[i, "theta"] = match_count[min(which(perc_rej_theta == min(perc_rej_theta)))]
    match_count_list[i, "tau"] = match_count[min(which(perc_rej_tau == min(perc_rej_tau)))]
}

save(match_count_list, file = '../Output_tree_rev/match_count_list_old.dat')


