set.seed(2023)

match_count <- 300
load("../Data/indexList_MAIN.RData")

adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)

perc_pval_match = vector(mode = "list", length = length(adjust_val))
p_val_df <- vector(mode = "list", length = length(adjust_val))
debug_spatial_ind = vector(mode = 'list', length = length(adjust_val))

for (k in 1:length(adjust_val)) {

    print(k)

    load(paste0('../Output/nullGridInfo_randomize/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output/origGridInfo_randomize/sim_orig_', k, '.dat'))

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


    #   v1 = sd(combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2, na.rm=TRUE)^2
    v1 = sd(combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2, na.rm=TRUE)^2

    #   rat_off = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
    rat_off = combinedMatchingSetupFix2$streets1 / combinedMatchingSetupFix2$streets2
    rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
    v2 = sd(rat_off, na.rm=TRUE)^2

    # t_stat = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
    #              - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)
    # t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
    #                   - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
    t_stat = combinedMatchingSetupFix2$spatialDiff
    t_stat_orig = sim_orig$DATA$spatialDiff

    row_num = 1
    perc_pval_match[[k]] = data.frame("num_match" = match_count,
                                    "perc_pval_less_05" = rep(NA, length(match_count)))
    p_val_df[[k]] = matrix(nrow = length(match_count), ncol = nrow(sim_orig$DATA))

    j = match_count
    print(paste0("Match Num: ", j))

    pval = rep(NA, nrow(sim_orig$DATA))

    for (ii in indexList_MAIN) {
    ## find matches
    # off_temp = sim_orig$DATA$n_off_1_prec[ii] + sim_orig$DATA$n_off_2_prec[ii]
    # ratio_temp = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
    #                     sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
    off_temp = sim_orig$DATA$streets1[ii] + sim_orig$DATA$streets2[ii]
    ratio_temp = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                     sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])

    stat_temp = t_stat_orig[ii]

    # dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2))^2/v1) +
    #                     ((ratio_temp - rat_off)^2 / v2))
    dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2))^2/v1) +
                         ((ratio_temp - rat_off)^2 / v2))

    w50 = order(dist_temp)[1:j]

    null_dist = t_stat[w50]
    pval[ii] = mean(null_dist > stat_temp, na.rm=TRUE)

    if(sum(is.nan(null_dist)) > 0) print(paste0("WARNING: ", j, ", ", ii))
    }

    perc_pval = mean(pval < 0.05, na.rm=TRUE)
    perc_pval_match[[k]]$perc_pval_less_05[row_num] = perc_pval
    p_val_df[[k]][row_num, ] = pval
    row_num = row_num + 1
}

save(p_val_df, file = "../Output/p_vals_match_rel/p_val_df_new_stat_FINAL.dat")
save(perc_pval_match, file = "../Output/p_vals_match_rel/perc_pval_match_new_stat_FINAL.dat")
