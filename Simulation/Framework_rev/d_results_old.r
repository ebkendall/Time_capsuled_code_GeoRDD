trialNum = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(trialNum)

file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))

tau = 0.5 

options(warn=1)
load("../Data/indexList_MAIN.RData")
n_buff_width = 8

perc_rejections_indiv = matrix(NA, nrow = n_buff_width, ncol = 4)

p_global = vector(mode = 'list', length = 4)
p_global_mat = matrix(nrow = n_buff_width, ncol = 3)
colnames(p_global_mat) = c("max", "mean", "med")
rownames(p_global_mat) = paste0(seq(from = 300, to = 1000, by = 100), "ft")
p_global[[1]] = p_global[[2]] = p_global[[3]] = p_global[[4]] = p_global_mat

load("../../NegControl/Output_tree_rev/match_count_list_old.dat")

test_stats <- function(gridPointValues, combinedMatchingSetupFix, w50) {
    
    t_stat_df = data.frame("t_stat" = c(-1),
                           "naive_pval" = c(-1))
    
    rowInd = 1
    
    for (jj in w50) {
        
        s1 = combinedMatchingSetupFix$DATA$streets1[jj]
        s2 = combinedMatchingSetupFix$DATA$streets2[jj]
        
        area1 = combinedMatchingSetupFix$DATA$area1[jj]
        area2 = combinedMatchingSetupFix$DATA$area2[jj]
        
        gridVals_ind_1 = combinedMatchingSetupFix$GRID_IND_1[[jj]]
        gridVals_ind_2 = combinedMatchingSetupFix$GRID_IND_2[[jj]]
        
        gridValues1 = gridPointValues[gridVals_ind_1]
        gridValues2 = gridPointValues[gridVals_ind_2]
        
        arr1 <- sum(gridValues1)
        arr2 <- sum(gridValues2)
        
        count1 = rpois(1, arr1)
        count2 = rpois(1, arr2)
        
        n = count1 + count2
        p = 0.5
        pval = NA
        
        if (count1 <= n/2) {
            pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
        } else {
            pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
        }
        
        t_stat_df[rowInd, ] = c(abs(count1 / area1  - count2 / area2), pval)
        rowInd = rowInd + 1
    }
    
    return(c(t_stat_df[,1]))
}

test_stats_orig <- function(gridPointValues, sim_orig, ii) {
    
    t_stat_df = matrix(nrow = 1, ncol = 2)
    colnames(t_stat_df) = c("t_stat", "naive_pval")
    
    s1 = sim_orig$DATA$streets1[ii]
    s2 = sim_orig$DATA$streets2[ii]
    
    area1 = sim_orig$DATA$area1[ii]
    area2 = sim_orig$DATA$area2[ii]
    
    gridValues1 = gridPointValues[sim_orig$GRID_IND_1[[ii]]]
    gridValues2 = gridPointValues[sim_orig$GRID_IND_2[[ii]]]
    
    arr1 <- sum(gridValues1)
    arr2 <- sum(gridValues2)
    
    count1 = rpois(1, arr1)
    count2 = rpois(1, arr2)
    
    n = count1 + count2
    p = 0.5
    
    if (count1 <= n/2) {
        pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
    } else {
        pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
    }
    
    t_stat_df[1, ] = c(abs(count1 / area1 - count2 / area2), pval)
    
    return(c(t_stat_df[,1]))
}

for(k in 1:n_buff_width) {
    load(paste0('../Output_noWater_rev/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output_noWater_rev/origGridInfo/sim_orig_', k, '.dat'))

    ## Now remove data points where these ratios are much different
    area_ratio = c(na.omit(sim_orig$DATA$area1 / sim_orig$DATA$area2))
    area_ratio[area_ratio < 1] = 1 / area_ratio[area_ratio < 1]
    wMax_a = max(area_ratio)
    wMin_a = min(area_ratio)
    
    street_ratio = c(na.omit(sim_orig$DATA$streets1 / sim_orig$DATA$streets2))
    street_ratio[street_ratio < 1] = 1 / street_ratio[street_ratio < 1]
    wMax_s = max(street_ratio)
    wMin_s = min(street_ratio)
    
    # Which null streets meet these requirements
    wRatioOk = which(combinedMatchingSetupFix$DATA$ratioArea > wMin_a &
                         combinedMatchingSetupFix$DATA$ratioArea < wMax_a & 
                         combinedMatchingSetupFix$DATA$ratioStreet > wMin_s &
                         combinedMatchingSetupFix$DATA$ratioStreet < wMax_s)
    
    combinedMatchingSetupFix2 = combinedMatchingSetupFix
    combinedMatchingSetupFix2$DATA = combinedMatchingSetupFix$DATA[wRatioOk,]
    combinedMatchingSetupFix2$GRID_IND_1 = combinedMatchingSetupFix$GRID_IND_1[wRatioOk]
    combinedMatchingSetupFix2$GRID_IND_2 = combinedMatchingSetupFix$GRID_IND_2[wRatioOk]
    
    # Null locations
    null_sum = combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2
    null_ratio = combinedMatchingSetupFix2$DATA$area1 / combinedMatchingSetupFix2$DATA$area2
    null_ratio[null_ratio < 1] = 1 / null_ratio[null_ratio < 1]
    
    sd1 = sd(null_sum, na.rm = T)
    sd2 = sd(null_ratio, na.rm = T)
    
    # Observed locations
    obs_sum = sim_orig$DATA$area1[indexList_MAIN] + sim_orig$DATA$area2[indexList_MAIN]
    obs_ratio = sim_orig$DATA$area1[indexList_MAIN] / sim_orig$DATA$area2[indexList_MAIN]
    obs_ratio[obs_ratio < 1] = 1 / obs_ratio[obs_ratio < 1]
    
    sd1_obs = sd(obs_sum, na.rm = T)
    sd2_obs = sd(obs_ratio, na.rm = T)
    
    # # Test statistics: Theta
    # t_stat = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
    #              - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)
    # t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
    #                   - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
    
    
    pval = matrix(NA, nrow = nrow(sim_orig$DATA), ncol = 4)
    stat_orig = matrix(NA, nrow = nrow(sim_orig$DATA), ncol = 4)
    
    match_vec = match_count_list[k,]
    global_null = vector(mode = 'list', length = 4)
    global_null[[1]] = global_null[[2]] = global_null[[3]] = global_null[[4]] = 
        matrix(NA, nrow = 164, ncol = match_vec[1])
    
    null_str_position = combinedMatchingSetupFix2$DATA[,c("precinct", "indigo", "juliet")]
    matched_null_str_loc = vector(mode = 'list', length = 164)
    
    for(ii in indexList_MAIN) {
        ## find matches
        off_temp = sim_orig$DATA$area1[ii] + sim_orig$DATA$area2[ii]
        ratio_temp = max(sim_orig$DATA$area1[ii] / sim_orig$DATA$area2[ii],
                         sim_orig$DATA$area2[ii] / sim_orig$DATA$area1[ii])
        
        # Remove relatively extreme values to improve the mahalanobis distance
        remove_extreme1 = which((null_sum > (0.25)*off_temp) & (null_sum < (1/0.25)*off_temp))
        remove_extreme2 = which((null_ratio > (0.25)*ratio_temp) & (null_ratio < (1/0.25)*ratio_temp))
        remove_extreme = intersect(remove_extreme1, remove_extreme2)
        
        null_sum_ii = null_sum[remove_extreme]
        null_ratio_ii = null_ratio[remove_extreme]
        v1_ii = sd(null_sum_ii, na.rm = T)^2
        v2_ii = sd(null_ratio_ii, na.rm = T)^2
        null_str_position_ii = null_str_position[remove_extreme, ]
        
        cms_ii = combinedMatchingSetupFix2
        cms_ii$DATA = combinedMatchingSetupFix2$DATA[remove_extreme,]
        cms_ii$GRID_IND_1 = combinedMatchingSetupFix2$GRID_IND_1[remove_extreme]
        cms_ii$GRID_IND_2 = combinedMatchingSetupFix2$GRID_IND_2[remove_extreme]
        
        dist_temp = sqrt(((off_temp - null_sum_ii)^2/v1_ii) + ((ratio_temp - null_ratio_ii)^2 / v2_ii))
        
        w50 = order(dist_temp)[1:match_vec[1]]
        matched_null_str_loc[[ii]] = null_str_position_ii[w50, ]
        
        stat_temp = p_val_temp = rep(NA, 4)
        for(s_name in 1:4) {
            load(file_names[s_name])
            gridPointValues = NULL
            
            if (s_name == 1) {
                gridPointValues = gridPointValues_hotspot * tau
            } else if (s_name == 2) {
                gridPointValues = gridPointValues_uniform
            } else if (s_name == 3) {
                gridPointValues = gridPointValues_cov_r
            } else if (s_name == 4) {
                gridPointValues = gridPointValues_cov_c_big
            } else {
                print("Incorrect input to start")
            }
            
            global_null[[s_name]][ii, ] = test_stats(gridPointValues, cms_ii, w50)
            stat_temp[s_name] = test_stats_orig(gridPointValues, sim_orig, ii)
            p_val_temp[s_name] = mean(global_null[[s_name]][ii, ] > stat_temp[s_name], na.rm=TRUE)
        }
        
        pval[ii, ] = p_val_temp
        stat_orig[ii, ] = stat_temp
        
        # t_stat_ii = t_stat[remove_extreme]
        # null_dist = t_stat_ii[w50]
        # global_null[ii,] = null_dist
        # stat_temp = t_stat_orig[ii]
        # pval[ii] = mean(null_dist > stat_temp, na.rm=TRUE)
    }
    # Individual test ---------------------------------------------------------
    print(paste0("Buffer: ", k+2, "00 ft Individual -------------------------"))
    perc_rej_k = apply(pval, 2, function(x){mean(x < 0.05, na.rm=TRUE)})
    perc_rejections_indiv[k,] = perc_rej_k
    print(perc_rej_k)
    
    
    # Global test -------------------------------------------------------------
    print(paste0("Buffer: ", k+2, "00 ft Global -----------------------------"))
    n_reps = 2000
    global_t_stat = vector(mode = 'list', length = 4)
    global_t_stat[[1]] = global_t_stat[[2]] = 
        global_t_stat[[3]] = global_t_stat[[4]] = data.frame("max_t_stat"  = rep(NA, n_reps),
                                                             "mean_t_stat" = rep(NA, n_reps),
                                                             "med_t_stat" = rep(NA, n_reps),
                                                             "max_loc" = rep(NA, n_reps))
    
    # Theta
    for (rep in 1:n_reps) {
        # This is the repetition to get the null distribution
        temp_loc = rep(NA, 164)
        temp_max = matrix(NA, nrow = 164, ncol = 4)
        for(ii in indexList_MAIN) {
            rand_ind = sample(c(1:match_vec[1]), 1)
            temp_max[ii,1] = global_null[[1]][ii, rand_ind]
            temp_max[ii,2] = global_null[[2]][ii, rand_ind]
            temp_max[ii,3] = global_null[[3]][ii, rand_ind]
            temp_max[ii,4] = global_null[[4]][ii, rand_ind]
            temp_loc[ii] = rand_ind
        }
        for(s in 1:4) {
            global_t_stat[[s]][rep, 1] = max(temp_max[,s], na.rm = T)
            global_t_stat[[s]][rep, 2] = mean(temp_max[,s], na.rm = T)
            global_t_stat[[s]][rep, 3] = median(temp_max[,s], na.rm = T)
            global_t_stat[[s]][rep, 4] = temp_loc[which.max(temp_max[,s])]
        }
    }
    
    t_stat_max = apply(stat_orig, 2, max, na.rm = T)
    t_stat_mean = colMeans(stat_orig, na.rm = T)
    t_stat_median = apply(stat_orig, 2, median, na.rm = T)
    
    p_val_global = matrix(NA, nrow = 4, ncol = 3)
    colnames(p_val_global) = c("max", "mean", "median")
    for(s in 1:4) {
        p_val_global[s,1] = mean(global_t_stat[[s]]$max_t_stat > t_stat_max[s])
        p_val_global[s,2] = mean(global_t_stat[[s]]$mean_t_stat > t_stat_mean[s])
        p_val_global[s,3] = mean(global_t_stat[[s]]$med_t_stat > t_stat_median[s])
    }
    
    print(p_val_global)
    
    p_global[[1]][k, ] = p_val_global[1,]
    p_global[[2]][k, ] = p_val_global[2,]
    p_global[[3]][k, ] = p_val_global[3,]
    p_global[[4]][k, ] = p_val_global[4,]
}

print("Individual test results ---------------------------------------------")
print(round(perc_rejections_indiv, digits = 4))

print("Global test results -------------------------------------------------")
print(p_global)

indiv_results_theta = perc_rejections_indiv
global_results_theta = p_global

save(indiv_results_theta, file = paste0('../Output_noWater_rev/test_results/indiv_test_', trialNum, '_old.rda'))
save(global_results_theta, file = paste0('../Output_noWater_rev/test_results/global_test_', trialNum, '_old.rda'))


all_indiv = array(NA, dim=c(1000, n_buff_width, 4))

all_global = vector(mode = 'list', length = 4)
all_global[[1]] = all_global[[2]] =
    all_global[[3]] = all_global[[4]] = vector(mode = 'list', length = 8)

for(t in 1:1000) {

    load(paste0('../Output_noWater_rev/test_results/indiv_test_', t, '_old.rda'))
    all_indiv[t,,] = indiv_results_theta

    load(paste0('../Output_noWater_rev/test_results/global_test_', t, '_old.rda'))

    for(i in 1:4) {
        for(j in 1:8) {
            all_global[[i]][[j]] = rbind(all_global[[i]][[j]],
                                         global_results_theta[[i]][j,,drop = F])
        }
    }

}

print("Individual Results")
avg_perc_rej = apply(all_indiv, 2:3, mean)
for(j in 1:8) {
    print(paste0((j+2)*100, "ft & ", round(avg_perc_rej[j,2], 3),
                 " & ", round(avg_perc_rej[j,3], 3),
                 " & ", round(avg_perc_rej[j,4], 3),
                 " & ", round(avg_perc_rej[j,1], 3)))
}

print("Global Results")
for(j in 1:8) {
    perc_rej_j = apply(all_global[[2]][[j]], 2, function(x){mean(x < 0.05)})
    perc_rej_j = c(perc_rej_j, apply(all_global[[3]][[j]], 2, function(x){mean(x < 0.05)}))
    perc_rej_j = c(perc_rej_j, apply(all_global[[4]][[j]], 2, function(x){mean(x < 0.05)}))
    perc_rej_j = c(perc_rej_j, apply(all_global[[1]][[j]], 2, function(x){mean(x < 0.05)}))

    print(paste0((j+2)*100, "ft & (", round(perc_rej_j[1], 3), ", ",
                 round(perc_rej_j[2], 3), ") & (",
                 round(perc_rej_j[4], 3), ", ",
                 round(perc_rej_j[5], 3), ") & (",
                 round(perc_rej_j[7], 3), ", ",
                 round(perc_rej_j[8], 3), ") & (",
                 round(perc_rej_j[10], 3), ", ",
                 round(perc_rej_j[11], 3), ")"))
}
