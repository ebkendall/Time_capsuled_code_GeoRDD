set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

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

for(k in 1:n_buff_width) {
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
    # t_stat_int_surface = int_surface_info[,3*(1:8)]
    # t_stat_off_surface = off_surface_info[,3*(1:8)]
    # t_stat_int_surface_orig = origData$str_surf$INT_SURFACE[,3*(1:8)]
    # t_stat_off_surface_orig = origData$str_surf$OFF_SURFACE[,3*(1:8)]
    
    pval = rep(NA, nrow(origData$str_info))
    pval_int = matrix(NA, nrow = nrow(origData$str_info), ncol = length(adjust_val))
    
    null_str_position = combinedMatchingSetupFix2[,c("precinct", "indigo", "juliet")]
    matched_null_str_loc = vector(mode = 'list', length = 164)
    
    n_matches = 2000
    store = array(NA, dim=c(164, n_matches, length(adjust_val)))
    
    for(kk in 1:length(adjust_val)) {
        y_obs = rep(NA, 164)
        y_obs_scale = rep(NA, 164)
        
        x_obs = matrix(NA, nrow = 164, ncol = 2)
        
        # Extract information at the observed boundaries -----------------------
        for(ii in indexList_MAIN) {
            y_obs[ii] = origData$str_surf$INT_SURFACE[ii, 3*kk]
            
            
            off_temp = origData$str_info$n_off_1_prec[ii] + origData$str_info$n_off_2_prec[ii]
            ratio_temp = max(origData$str_info$n_off_1_prec[ii] / origData$str_info$n_off_2_prec[ii],
                             origData$str_info$n_off_2_prec[ii] / origData$str_info$n_off_1_prec[ii])
            x_obs[ii,] = c(ratio_temp, off_temp)
        }
        
        
    }
    
    for(ii in indexList_MAIN) {
        ## find matches
        off_temp = origData$str_info$n_off_1_prec[ii] + origData$str_info$n_off_2_prec[ii]
        ratio_temp = max(origData$str_info$n_off_1_prec[ii] / origData$str_info$n_off_2_prec[ii],
                         origData$str_info$n_off_2_prec[ii] / origData$str_info$n_off_1_prec[ii])
        
        # Remove relatively extreme values to improve the mahalanobis distance
        remove_extreme = which((null_ratio < 2*ratio_temp) & (null_ratio > (1/2)*ratio_temp))
        null_sum_ii = null_sum[remove_extreme]
        null_ratio_ii = null_ratio[remove_extreme]
        t_stat_ii = t_stat[remove_extreme]
        t_stat_int_surface_ii = t_stat_int_surface[remove_extreme, ]
        v1_ii = sd(null_sum_ii, na.rm = T)^2
        v2_ii = sd(null_ratio_ii, na.rm = T)^2
        null_str_position_ii = null_str_position[remove_extreme, ]
        
        
        dist_temp = sqrt(((off_temp - null_sum_ii)^2/v1_ii) + ((ratio_temp - null_ratio_ii)^2 / v2_ii))
        
        w50 = order(dist_temp)[1:n_matches]
        
        null_dist = t_stat_ii[w50]
        global_null[ii,] = null_dist
        
        ind_order = order(dist_temp)[1:min(5000, length(dist_temp))] # sample from the best 5,000
        t_stat_ind = t_stat_ii[ind_order]
        prob_w50 = exp(-dist_temp[ind_order])
        prob_w50 = prob_w50 / sum(prob_w50)
        prob_dist = cbind(ind_order, t_stat_ind, prob_w50)
        colnames(prob_dist) = c("w50_ind", "w50_tstat", "w50_prob")
        global_samp[[ii]] = prob_dist
        
        stat_temp = t_stat_orig[ii]
        pval[ii] = mean(null_dist > stat_temp, na.rm=TRUE)
        
        matched_null_str_loc[[ii]] = null_str_position_ii[w50, ]
        
        for(kk in 1:length(adjust_val)) {
            null_dist_int = t_stat_int_surface_ii[w50, ]
            global_null_int_surf[[kk]][ii, ] = null_dist_int[,kk]
            
            t_stat_ind_kk = t_stat_int_surface_ii[ind_order, kk]
            prob_dist_kk = cbind(ind_order, t_stat_ind_kk, prob_w50)
            colnames(prob_dist_kk) = c("w50_ind", "w50_tstat", "w50_prob")
            global_samp_int_surf[[kk]][[ii]] = prob_dist_kk
            
            pval_int[ii, kk] = mean(null_dist_int[,kk] > t_stat_int_surface_orig[ii,kk], na.rm=TRUE)
        }
    }
    # Individual test ---------------------------------------------------------
    print(paste0("Buffer: ", k+2, "00 ft Individual -------------------------"))
    print(paste0("Theta: ", round(mean(pval < 0.05, na.rm=TRUE), digits = 4)))
    print("Tau")
    tau_pval = apply(pval_int, 2, function(x){mean(x < 0.05, na.rm=TRUE)})
    print(round(tau_pval, digits = 4))    
    
    perc_rejections_indiv[k] = mean(pval < 0.05, na.rm=TRUE)
    perc_rejections_indiv_surf[k, ] = tau_pval
    
    
    # Global test -------------------------------------------------------------
    print(paste0("Buffer: ", k+2, "00 ft Global -----------------------------"))
    n_reps = 300
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
            rand_ind = sample(c(1:n_matches), 1)
            temp_max[ii] = global_null[ii, rand_ind]
            temp_loc[ii] = rand_ind
            # temp_max[ii] = sample(x = global_samp[[ii]][,"w50_tstat"], size = 1,
            #                          prob = global_samp[[ii]][,"w50_prob"])
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
    print("Theta")
    print(paste0("Max = ", round(p_val_global[1], digits = 4), 
                 ", Mean = ", round(p_val_global[2], digits = 4),
                 ", Med. = ", round(p_val_global[3], digits = 4)))
    
    p_global[k, ] = p_val_global
    
    # Tau
    global_match_loc = vector(mode = 'list', length = length(adjust_val))
    for(av in 1:length(adjust_val)) {
        global_match_loc[[av]] = matrix(0,nrow = n_reps, ncol = ncol(null_str_position))
        
        for (rep in 1:n_reps) {
            # This is the repetition to get the null distribution
            temp_loc = temp_max = rep(NA, 164)
            for(ii in indexList_MAIN) {
                rand_ind = sample(c(1:n_matches), 1)
                temp_max[ii] = global_null_int_surf[[av]][ii, rand_ind]
                temp_loc[ii] = rand_ind
                # temp_max[ii] = sample(x = global_samp_int_surf[[av]][[ii]][,"w50_tstat"], size = 1,
                #                          prob = global_samp_int_surf[[av]][[ii]][,"w50_prob"])
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
    
    p_val_global_surf = matrix(nrow = 3, ncol = length(adjust_val))
    global_match_loc_orig = rep(NA, length(adjust_val))
    
    for(av in 1:length(adjust_val)) {
        for(m in 1:3) {
            if(m == 1) {
                t_stat_kk = max(t_stat_int_surface_orig[,av], na.rm = T)
                global_null_kk = global_t_stat_int_surf[[av]]$max_t_stat
            } else if(m == 2) {
                t_stat_kk = mean(t_stat_int_surface_orig[,av], na.rm = T)
                global_null_kk = global_t_stat_int_surf[[av]]$mean_t_stat
            } else {
                t_stat_kk = median(t_stat_int_surface_orig[,av], na.rm = T) 
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
    print("Tau")
    print(round(p_val_global_surf, digits = 4))
    p_global_surf[[1]][k, ] = p_val_global_surf[1,]
    p_global_surf[[2]][k, ] = p_val_global_surf[2,]
    p_global_surf[[3]][k, ] = p_val_global_surf[3,]
}

print("Individual test results ---------------------------------------------")
print("Theta")
print(round(perc_rejections_indiv, digits = 4))
print("Tau")
print(round(perc_rejections_indiv_surf, digits = 4))

print("Global test results -------------------------------------------------")
print("Theta")
print(round(p_global, digits = 4))
print("Tau")
names(p_global_surf) = c("max", "mean", "med")
print(p_global_surf)

curr_results = list(perc_rejections_indiv = perc_rejections_indiv,
                    perc_rejections_indiv_surf = perc_rejections_indiv_surf,
                    p_global = p_global,
                    p_global_surf = p_global_surf)