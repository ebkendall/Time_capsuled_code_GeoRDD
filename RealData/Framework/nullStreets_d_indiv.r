set.seed(2023)

load("../../NegControl/Output_tree/combination/match_count_list.dat")
load("../Data/indexList_MAIN.RData")

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)
match_count = c(match_count_list[[1]][[1]], match_count_list[[2]][[1]],
                match_count_list[[3]][[1]], match_count_list[[4]][[1]],
                match_count_list[[5]][[1]], match_count_list[[6]][[1]],
                match_count_list[[7]][[1]], match_count_list[[8]][[1]])

perc_pval_match = vector(mode = "list", length = length(adjust_val))
p_val_df <- vector(mode = "list", length = length(adjust_val))

int_surface_pval = vector(mode = "list", length = 8)
int_surface_pval_vec = vector(mode = "list", length = 8)

for (k in 1:8) {
    print(k)

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

    t_stat = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
                 - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)
    t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
                      - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)

    t_stat_int_surface = int_surface_info[,1:8]
    t_stat_int_surface_orig = sim_orig$INT_SURFACE[,1:8]

    row_num = 1
    match_vec = match_count_list[[k]]
    print("Matches for old method")
    print(match_vec[[1]])
    print("Matches for adjustment values")
    print(match_vec[[2]])

    perc_pval_match[[k]] = data.frame("num_match" = match_count,
                                    "perc_pval_less_05" = rep(NA, length(match_count)))
    p_val_df[[k]] = matrix(nrow = length(match_count), ncol = nrow(sim_orig$DATA))
    int_surface_pval[[k]] = data.frame("adjust_val" = adjust_val,
                                        "perc_pval_less_05" = rep(NA, length(adjust_val)))

    pval = rep(NA, nrow(sim_orig$DATA))
    pval_int = matrix(NA, nrow = nrow(sim_orig$DATA), ncol = length(adjust_val))

    for (ii in indexList_MAIN) {
        ## find matches
        off_temp = sim_orig$DATA$n_off_1_prec[ii] + sim_orig$DATA$n_off_2_prec[ii]
        ratio_temp = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
                            sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
        # off_temp = sim_orig$DATA$streets1[ii] + sim_orig$DATA$streets2[ii]
        # ratio_temp = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
        #                 sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])

        stat_temp = t_stat_orig[ii]

        dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2))^2/v1) +
                            ((ratio_temp - rat_off)^2 / v2))
        # dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2))^2/v1) +
        #                     ((ratio_temp - combinedMatchingSetupFix2$ratioStreet)^2 / v2))

        w50 = order(dist_temp)[1:match_vec[[1]]]

        null_dist = t_stat[w50]
        pval[ii] = mean(null_dist > stat_temp, na.rm=TRUE)

        for(kk in 1:ncol(pval_int)) {
            w50_k = order(dist_temp)[1:match_vec[[2]][kk]]
            null_dist_int = t_stat_int_surface[w50_k, ]
            pval_int[ii, kk] = mean(null_dist_int[,kk] > t_stat_int_surface_orig[ii,kk], na.rm=TRUE)
        }

        if(sum(is.nan(null_dist)) > 0) print(paste0("WARNING: ", j, ", ", ii))
    }

    perc_pval = mean(pval < 0.05, na.rm=TRUE)
    perc_pval_match[[k]]$perc_pval_less_05[row_num] = perc_pval
    p_val_df[[k]][row_num, ] = pval

    int_surface_pval_vec[[k]] = pval_int

    for(jj in 1:length(adjust_val)) {
        int_surface_pval[[k]]$perc_pval_less_05[jj] = mean(pval_int[,jj] < 0.05, na.rm=TRUE)
    }
    row_num = row_num + 1
}

save(p_val_df, file = "../Output/p_vals_match_rel/p_val_df_new_stat_FINAL.dat")
save(perc_pval_match, file = "../Output/p_vals_match_rel/perc_pval_match_new_stat_FINAL.dat")
save(int_surface_pval, file = paste0("../Output/p_vals_match_rel/int_surface_pval.dat"))
save(int_surface_pval_vec, file = paste0("../Output/p_vals_match_rel/int_surface_pval_vec.dat"))