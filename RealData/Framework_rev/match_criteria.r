match_percent = seq(0, 0.99, by = 0.01)
set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

indiv_results_theta = indiv_results_tau = matrix(NA, nrow = n_buff_width, 
                                                 ncol = length(match_percent))

for(k in 1:n_buff_width) {
    print(paste0("buff ", k + 2))
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
    
    # Test statistics: Tau (NEW)
    t_stat_int_surface = int_surface_info
    t_stat_off_surface = off_surface_info
    t_stat_int_surface_orig = origData$str_surf$INT_SURFACE
    t_stat_off_surface_orig = origData$str_surf$OFF_SURFACE
    
    n_null = nrow(combinedMatchingSetupFix2)
    n_matches = 2000
    
    Y_theta = Y_tau = rep(NA, length(indexList_MAIN))
    X = matrix(NA, length(indexList_MAIN), 2)
    
    counter = 1
    kk = 2
    for(ii in indexList_MAIN) {
        
        # Match on crime/offenses
        off_temp = origData$str_info$n_off_1_prec[ii] + origData$str_info$n_off_2_prec[ii]
        ratio_temp = max(origData$str_info$n_off_1_prec[ii] / origData$str_info$n_off_2_prec[ii],
                         origData$str_info$n_off_2_prec[ii] / origData$str_info$n_off_1_prec[ii])
        
        Y_theta[counter] = t_stat_orig[ii]
        Y_tau[counter] = abs(t_stat_int_surface_orig[ii,3*kk-2] / t_stat_off_surface_orig[ii,3*kk-2] -
                             t_stat_int_surface_orig[ii,3*kk-1] / t_stat_off_surface_orig[ii,3*kk-1])
        X[counter,] = c(ratio_temp, off_temp)
        
        counter = counter + 1
    }
    
    for(p in 1:length(match_percent)) {
        m_p = match_percent[p]; print(m_p)
        store_theta = store_tau = matrix(NA, nrow = length(indexList_MAIN), ncol = n_matches)
        ind_p_theta = ind_p_tau = rep(NA, length(indexList_MAIN))
        
        for(ii in 1:nrow(X)) {
            
            # Match on crime/offenses
            off_temp = X[ii,2]
            ratio_temp = X[ii,1]
            
            if(m_p == 0) {
                wAll = 1:nrow(combinedMatchingSetupFix2)
            } else {
                w1 = which(null_sum > m_p*off_temp & null_sum < (1/m_p)*off_temp)
                w2 = which(null_ratio > m_p*ratio_temp & null_ratio < (1/m_p)*ratio_temp)
                
                wAll = intersect(w1, w2)
            }
            
            if (length(wAll) > 10) {
                
                t_NULL_theta = t_stat[wAll]
                t_NULL_theta = t_NULL_theta[which(t_NULL_theta > 0)]
                store_theta[ii,] = sample(t_NULL_theta, n_matches, replace=TRUE)
                
                t_NULL_tau = abs(t_stat_int_surface[wAll, 3*kk-2] / t_stat_off_surface[wAll, 3*kk-2]
                                    - t_stat_int_surface[wAll, 3*kk-1] / t_stat_off_surface[wAll, 3*kk-1])
                t_NULL_tau = t_NULL_tau[which(t_NULL_tau > 0)]
                store_tau[ii,] = sample(t_NULL_tau, n_matches, replace=TRUE)
            }
            
        }
        
        for (jj in 1:nrow(X)) {
            ind_p_theta[jj] = mean(store_theta[jj,] > Y_theta[jj])
            ind_p_tau[jj] = mean(store_tau[jj,] > Y_tau[jj])
        }
        
        indiv_results_theta[k,p] = mean(ind_p_theta < .05, na.rm=TRUE)
        indiv_results_tau[k,p] = mean(ind_p_tau < .05, na.rm=TRUE)
    }
}
save(indiv_results_theta, file = '../Output_rev/indiv_results_theta.rda')
save(indiv_results_tau, file = '../Output_rev/indiv_results_tau.rda')


# Plotting the results
library(tidyverse, quietly = T)
library(gridExtra, quietly = T)
p_theta = p_tau = vector(mode = 'list', length = 8)
load('../Output_rev/indiv_results_theta.rda')
load('../Output_rev/indiv_results_tau.rda')
for(i in 1:8) {
    p_val_theta = data.frame(y = indiv_results_theta[i,],
                             x = match_percent)
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
                             x = match_percent)
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
pdf("../_visualizations/Plots_rev/match_theta.pdf", onefile = T)
grid.arrange(p_theta[[1]], p_theta[[2]], ncol = 1, nrow = 2)
grid.arrange(p_theta[[3]], p_theta[[4]], ncol = 1, nrow = 2)
grid.arrange(p_theta[[5]], p_theta[[6]], ncol = 1, nrow = 2)
grid.arrange(p_theta[[7]], p_theta[[8]], ncol = 1, nrow = 2)
dev.off()

pdf("../_visualizations/Plots_rev/match_tau.pdf", onefile = T)
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
    
    diff_theta = abs(perc_rej_theta - 0.05)
    diff_tau = abs(perc_rej_tau - 0.05)
    
    match_count_list[i, "theta"] = match_percent[min(which(diff_theta == min(diff_theta)))]
    match_count_list[i, "tau"] = match_percent[min(which(diff_tau == min(diff_tau)))]
}

save(match_count_list, file = '../Output_rev/match_count_list.dat')


