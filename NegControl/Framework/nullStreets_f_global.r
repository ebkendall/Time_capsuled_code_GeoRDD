# Step Outline -----------------------------------------------------------------
# Step 1: save all the matches in the global_null list for each surface
# Step 2: iterate through and grab the max
# Step 3: get a p-value for each surface and each iterate for each buffer width
# ------------------------------------------------------------------------------

load("../Data/indexList_MAIN.RData")

match_count <- c(160, 880, 480, 160, 1020, 140, 140, 300)

set.seed(2023)

`%notin%` <- Negate(`%in%`)

# Step 1 -----------------------------------------------------------------------

# Run this for the different surface types
global_null = vector(mode = "list", length = 10)

for (k in 3:10) {
 
    print(k)
    n_matches = match_count[k-2]
    print(paste0("Match Num: ", n_matches))
  
    load(paste0('../Output_tree/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output_tree/origGridInfo/sim_orig_', k, '.dat'))

    global_null[[k]] = matrix(nrow = max(indexList_MAIN), ncol = n_matches)

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
    
    # =====================================================================
    
    v1 = sd(combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2, na.rm=TRUE)^2
    v2 = sd(combinedMatchingSetupFix2$ratioStreet, na.rm=TRUE)^2

    t_stat_streets = abs(combinedMatchingSetupFix2$count1 / combinedMatchingSetupFix2$streets1
                         - combinedMatchingSetupFix2$count2 / combinedMatchingSetupFix2$streets2)

    for(ii in indexList_MAIN) {
        ## find matches
        streets_temp = sim_orig$DATA$streets1[ii] + sim_orig$DATA$streets2[ii]
        ratio_temp = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                        sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])
        
        dist_temp = sqrt(((streets_temp - (combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2))^2/v1) +
                            ((ratio_temp - combinedMatchingSetupFix2$ratioStreet)^2 / v2))
      
        w50 = order(dist_temp)[1:n_matches]

        null_dist = t_stat_streets[w50]
        global_null[[k]][ii,] = null_dist
    }

}

# Step 2 -----------------------------------------------------------------------
global_t_stat <- vector(mode = "list", length = 13)

for(k in 3:10) {

    print(k)
    n_matches = match_count[k-2]
    print(paste0("Match Num: ", n_matches))

    global_t_stat[[k]] = data.frame("max_t_stat" = rep(NA, n_matches),
                                    "max_loc" = rep(NA, n_matches))
    
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
        global_t_stat[[k]][rep, 2] = temp_loc[which.max(temp_max)]
    }

}

save(global_t_stat, file = "../Output_tree/combination/global_t_stat_FINAL.dat")

# Step 3 -----------------------------------------------------------------------

p_val_df = rep(NA, 10)

for(k in 3:10) {
  
  load(paste0('../Output_tree/origGridInfo/sim_orig_', k,".dat"))
  t_stat_streets_orig = abs(sim_orig$DATA$count1 / sim_orig$DATA$streets1
                            - sim_orig$DATA$count2 / sim_orig$DATA$streets2)
  
  t_stat = max(t_stat_streets_orig, na.rm = T)

  p_val_df[k] = mean(global_t_stat[[k]]$max_t_stat > t_stat)
}
print(p_val_df)