load("../Data/indexList_MAIN.RData")
load("../../NegControl/Output_tree/combination/match_count_list.dat")
set.seed(2023)

`%notin%` <- Negate(`%in%`)

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)
buff_val = 3:10
# ------------------------------------------------------------------------------
# Step 1 -----------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Run this for the different surface types
global_null = vector(mode = "list", length = length(buff_val))
global_null_int_surf = vector(mode = "list", length = length(buff_val))

for (k in 1:length(buff_val)) {
    
    match_vec = match_count_list[[k]]

    print(k)
    n_matches = match_vec[[1]]

    global_null[[k]] = matrix(nrow = 164, ncol = n_matches)
    global_null_int_surf[[k]] = vector(mode = 'list', length = length(adjust_val))
    for(jj in 1:length(global_null_int_surf[[k]])) {
        global_null_int_surf[[k]][[jj]] = matrix(nrow = 164, ncol = match_vec[[2]][jj])
    }
  
    load(paste0('../Output/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output/origGridInfo/sim_orig_', k, '.dat'))

    ## Now remove data points where these ratios are much different
    area_ratio = c(na.omit(sim_orig$DATA$area1 / sim_orig$DATA$area2))
    area_ratio[area_ratio < 1] = 1 / area_ratio[area_ratio < 1]
    wMax_a = max(area_ratio)
    wMin_a = min(area_ratio)
    
    street_ratio = c(na.omit(sim_orig$DATA$streets1 / sim_orig$DATA$streets2))
    street_ratio[street_ratio < 1] = 1 / street_ratio[street_ratio < 1]
    wMax_s = max(street_ratio)
    wMin_s = min(street_ratio)
    
    wRatioOk = which(combinedMatchingSetupFix$DATA$ratioArea > wMin_a &
                         combinedMatchingSetupFix$DATA$ratioArea < wMax_a & 
                         combinedMatchingSetupFix$DATA$ratioStreet > wMin_s &
                         combinedMatchingSetupFix$DATA$ratioStreet < wMax_s)
    
    combinedMatchingSetupFix2 = combinedMatchingSetupFix$DATA[wRatioOk,]
    int_surface_info = combinedMatchingSetupFix$INT_SURFACE[wRatioOk,]
    
    # Wherever there is a 0 for the offense count, everything gets scaled by 1
    which_zeros = which(combinedMatchingSetupFix2$n_off_1 == 0 | combinedMatchingSetupFix2$n_off_2 == 0)
    combinedMatchingSetupFix2$n_arr_1[which_zeros] = combinedMatchingSetupFix2$n_arr_1[which_zeros] + 1 
    combinedMatchingSetupFix2$n_arr_2[which_zeros] = combinedMatchingSetupFix2$n_arr_2[which_zeros] + 1 
    combinedMatchingSetupFix2$n_off_1[which_zeros] = combinedMatchingSetupFix2$n_off_1[which_zeros] + 1 
    combinedMatchingSetupFix2$n_off_2[which_zeros] = combinedMatchingSetupFix2$n_off_2[which_zeros] + 1
    
    which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
    sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1
    
    v1 = sd(combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2, na.rm=TRUE)^2
    # v1 = sd(combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2, na.rm=TRUE)^2
    
    rat_off = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
    rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
    v2 = sd(rat_off, na.rm=TRUE)^2
    # v2 = sd(combinedMatchingSetupFix2$ratioStreet, na.rm=TRUE)^2
    
    # =====================================================================

    t_stat = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
                 - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)

    t_stat_int_surface = int_surface_info[,1:8]
    
    for(ii in indexList_MAIN) {
        ## find matches
        off_temp = sim_orig$DATA$n_off_1_prec[ii] + sim_orig$DATA$n_off_2_prec[ii]
        ratio_temp = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
                         sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
        
        dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2))^2/v1) +
                             ((ratio_temp - rat_off)^2 / v2))

        w50 = order(dist_temp)[1:n_matches]

        null_dist = t_stat[w50]

        global_null[[k]][ii,] = null_dist
        
        for(kk in 1:length(adjust_val)) {
            w50_k = order(dist_temp)[1:match_vec[[2]][kk]]
            null_dist_int = c(t_stat_int_surface[w50_k, kk])
            global_null_int_surf[[k]][[kk]][ii, ] = null_dist_int
        }
    }

}

# save(global_null, file = paste0("../Output/Global/global_null_FINAL.dat"))
# save(global_null_int_surf, file = paste0("../Output/Global/global_null_int_surf_FINAL.dat"))


# ------------------------------------------------------------------------------
# Step 2 -----------------------------------------------------------------------
# ------------------------------------------------------------------------------

global_t_stat <- vector(mode = "list", length = length(buff_val))
global_t_stat_int_surf <- vector(mode = "list", length = length(buff_val))

for(k in 1:length(buff_val)) {

    match_vec = match_count_list[[k]]
    
    print(k)
    n_matches = match_vec[[1]]

    global_t_stat[[k]] = data.frame("max_t_stat"  = rep(NA, n_matches),
                                    "mean_t_stat" = rep(NA, n_matches),
                                    "med_t_stat" = rep(NA, n_matches),
                                    "max_loc" = rep(NA, n_matches))
    global_t_stat_int_surf[[k]] = vector(mode='list',length=length(adjust_val))
    for(jj in 1:length(global_t_stat_int_surf[[k]])) {
        match_jj = match_vec[[2]][jj]
        global_t_stat_int_surf[[k]][[jj]] = data.frame("max_t_stat"  = rep(NA, match_jj),
                                                       "mean_t_stat" = rep(NA, match_jj),
                                                       "med_t_stat" = rep(NA, match_jj),
                                                       "max_loc" = rep(NA, match_jj))
    }
    
    for (rep in 1:n_matches) {
        # This is the repetition to get the null distribution
        temp_loc = temp_max = c()
        myInd = 1
        for(ii in indexList_MAIN) {
            rand_ind = sample(c(1:n_matches), 1)

            temp_loc[myInd] = ii
            temp_max[myInd] = global_null[[k]][ii, rand_ind]
            myInd = myInd + 1
        }
        global_t_stat[[k]][rep, 1] = max(temp_max, na.rm = T)
        global_t_stat[[k]][rep, 2] = mean(temp_max, na.rm = T)
        global_t_stat[[k]][rep, 3] = median(temp_max, na.rm = T)
        global_t_stat[[k]][rep, 4] = temp_loc[which.max(temp_max)]
    }
    
    for(av in 1:length(adjust_val)) {
        print(paste0("av: ", av))
        match_av = match_vec[[2]][av]
        for (rep in 1:match_av) {
            # This is the repetition to get the null distribution
            temp_loc = temp_max = c()
            myInd = 1
            for(ii in indexList_MAIN) {
                rand_ind = sample(c(1:match_av), 1)
                
                temp_loc[myInd] = ii
                temp_max[myInd] = global_null_int_surf[[k]][[av]][ii, rand_ind]
                myInd = myInd + 1
            }
            global_t_stat_int_surf[[k]][[av]][rep, 1] = max(temp_max, na.rm = T)
            global_t_stat_int_surf[[k]][[av]][rep, 2] = mean(temp_max, na.rm = T)
            global_t_stat_int_surf[[k]][[av]][rep, 3] = median(temp_max, na.rm = T)
            global_t_stat_int_surf[[k]][[av]][rep, 4] = temp_loc[which.max(temp_max)]
        }
    }

}

# save(global_t_stat, file = paste0("../Output/Global/global_t_stat_FINAL.dat"))
# save(global_t_stat_int_surf, file = paste0("../Output/Global/global_t_stat_int_surf_FINAL.dat"))

# Computing the p-value
p_val_df = matrix(nrow =length(buff_val), ncol = 3)
colnames(p_val_df) = c("max", "mean", "median")
p_val_df_int_surf = vector(mode = 'list', length = 3)
names(p_val_df_int_surf) = c("max", "mean", "median")
p_val_df_int_surf[[1]] = p_val_df_int_surf[[2]] = p_val_df_int_surf[[3]] = 
    matrix(nrow=length(buff_val), ncol = length(adjust_val))

pdf('test_global_tStat_tau.pdf')
par(mfrow = c(3,1))
for(k in 1:length(buff_val)) {
    
    load(paste0('../Output/origGridInfo/sim_orig_', k,".dat"))
    which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
    sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1
    
    t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
                      - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
    
    t_stat_max = max(t_stat_orig, na.rm = T)
    t_stat_mean = mean(t_stat_orig, na.rm = T)
    t_stat_median = median(t_stat_orig, na.rm = T)
    
    p_val_df[k, 'max'] = mean(global_t_stat[[k]]$max_t_stat > t_stat_max)
    p_val_df[k, 'mean'] = mean(global_t_stat[[k]]$mean_t_stat > t_stat_mean)
    p_val_df[k, 'median'] = mean(global_t_stat[[k]]$med_t_stat > t_stat_median)
    
    
    t_stat_int_surface_orig = sim_orig$INT_SURFACE[,1:8]
    for(kk in 1:length(adjust_val)) {
        for(m in 1:3) {
            if(m == 1) {
                t_stat_kk = max(t_stat_int_surface_orig[,kk], na.rm = T)
                global_null_kk = global_t_stat_int_surf[[k]][[kk]]$max_t_stat
            } else if(m == 2) {
                t_stat_kk = mean(t_stat_int_surface_orig[,kk], na.rm = T)
                global_null_kk = global_t_stat_int_surf[[k]][[kk]]$mean_t_stat
            } else {
                t_stat_kk = median(t_stat_int_surface_orig[,kk], na.rm = T)   
                global_null_kk = global_t_stat_int_surf[[k]][[kk]]$med_t_stat
            }
            
            p_val_df_int_surf[[m]][k, kk] = mean(global_null_kk > t_stat_kk, na.rm = T)
            hist(global_null_kk, main = paste0(colnames(p_val_df)[m], ': buff(', 
                                               buff_val[k], '), adjust(', adjust_val[kk], ')'), 
                 breaks = sqrt(length(global_null_kk)),
                 xlab = paste0("na: ", sum(is.na(global_null_kk))), 
                 ylab = round(t_stat_kk, digits = 5))
            abline(v = t_stat_kk, lwd = 2, col = 'red')
        }
    }
    
}
dev.off()

pdf('test_global_tStat_theta.pdf')
par(mfrow = c(3,1))
for(k in 1:length(buff_val)) {
    
    load(paste0('../Output/origGridInfo/sim_orig_', k,".dat"))
    which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
    sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
    sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1
    
    t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
                      - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
    
    t_stat_max = max(t_stat_orig, na.rm = T)
    t_stat_mean = mean(t_stat_orig, na.rm = T)
    t_stat_median = median(t_stat_orig, na.rm = T)
    
    p_val_df[k, 'max'] = mean(global_t_stat[[k]]$max_t_stat > t_stat_max)
    hist(global_t_stat[[k]]$max_t_stat, main = paste0('max: buff(', buff_val[k], ')'), 
         breaks = sqrt(length(global_t_stat[[k]]$max_t_stat)),
         xlab = paste0("na: ", sum(is.na(global_t_stat[[k]]$max_t_stat))), 
         ylab = round(t_stat_max, digits = 5))
    abline(v = t_stat_max, lwd = 2, col = 'red')
    
    p_val_df[k, 'mean'] = mean(global_t_stat[[k]]$mean_t_stat > t_stat_mean)
    hist(global_t_stat[[k]]$mean_t_stat, main = paste0('mean: buff(', buff_val[k], ')'), 
         breaks = sqrt(length(global_t_stat[[k]]$mean_t_stat)),
         xlab = paste0("na: ", sum(is.na(global_t_stat[[k]]$mean_t_stat))), 
         ylab = round(t_stat_mean, digits = 5))
    abline(v = t_stat_mean, lwd = 2, col = 'red')
    
    p_val_df[k, 'median'] = mean(global_t_stat[[k]]$med_t_stat > t_stat_median)
    hist(global_t_stat[[k]]$med_t_stat, main = paste0('median: buff(', buff_val[k], ')'), 
         breaks = sqrt(length(global_t_stat[[k]]$med_t_stat)),
         xlab = paste0("na: ", sum(is.na(global_t_stat[[k]]$med_t_stat))), 
         ylab = round(t_stat_median, digits = 5))
    abline(v = t_stat_median, lwd = 2, col = 'red')
}
dev.off()
