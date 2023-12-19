set.seed(100)

# match_count <- c(160, 880, 480, 160, 1020, 140, 140, 300) 
match_count = 300
load("../Data/indexList_MAIN.RData")

adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)

perc_pval_match = vector(mode = "list", length = length(adjust_val))
p_val_df <- vector(mode = "list", length = length(adjust_val))

for (k in 1:length(adjust_val)) {
  print(k)
  
  load(paste0('../Output_tree/nullGridInfo/combinedMatchingSetup', k, ".dat"))
  load(paste0('../Output_tree/origGridInfo/sim_orig_', k, '.dat'))

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
  
  v1 = sd(combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2, na.rm=TRUE)^2
  v2 = sd(combinedMatchingSetupFix2$ratioStreet, na.rm=TRUE)^2

  # t_stat_streets = abs(combinedMatchingSetupFix2$count1 / combinedMatchingSetupFix2$streets1
  #                      - combinedMatchingSetupFix2$count2 / combinedMatchingSetupFix2$streets2)
  # t_stat_streets_orig = abs(sim_orig$DATA$count1 / sim_orig$DATA$streets1
  #                           - sim_orig$DATA$count2 / sim_orig$DATA$streets2)
  t_stat_streets = combinedMatchingSetupFix2$spatialDiff
  t_stat_streets_orig = sim_orig$DATA$spatialDiff
  
  row_num = 1
  perc_pval_match[[k]] = data.frame("num_match" = match_count,
                                    "perc_pval_less_05" = rep(NA, length(match_count)))
  p_val_df[[k]] = matrix(nrow = length(match_count), ncol = nrow(sim_orig$DATA))

  j = match_count
  print(paste0("Match Num: ", j))

  pval = rep(NA, nrow(sim_orig$DATA))

  for (ii in indexList_MAIN) {
    ## find matches
    streets_temp = sim_orig$DATA$streets1[ii] + sim_orig$DATA$streets2[ii]
    ratio_temp = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                      sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])
    stat_temp = t_stat_streets_orig[ii]
    
    dist_temp = sqrt(((streets_temp - (combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2))^2/v1) +
                        ((ratio_temp - combinedMatchingSetupFix2$ratioStreet)^2 / v2))
    
    w50 = order(dist_temp)[1:j]
    
    null_dist = t_stat_streets[w50]
    pval[ii] = mean(null_dist > stat_temp)
  }

  perc_pval = mean(pval < 0.05, na.rm=TRUE)
  perc_pval_match[[k]]$perc_pval_less_05[row_num] = perc_pval
  p_val_df[[k]][row_num, ] = pval
  row_num = row_num + 1
}

save(p_val_df, file = paste0("../Output_tree/p_vals_match_rel/p_val_df_FINAL.dat"))
save(perc_pval_match, file = paste0("../Output_tree/p_vals_match_rel/perc_pval_match_FINAL.dat"))

