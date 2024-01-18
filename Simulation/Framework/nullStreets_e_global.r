load("../Data/indexList_MAIN.RData")

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
    
    return(t_stat_df)
}

save_type = c("HotSpot/", "Uniform/", "Random/", "Correlated/")

load("../../NegControl/Output_tree/combination/match_count_list.dat")

# for(trialNum in 1:1000) {
    trialNum = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    set.seed(trialNum)

    `%notin%` <- Negate(`%in%`)

    file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))
    tau = 0.5  
    # Step 1 -----------------------------------------------------------------------

    for (s_name in 1:4) {
        # Run this for the different surface types
        global_null = vector(mode = "list", length = 10)

        load(file_names[s_name])
        gridPointValues = NULL

        if (s_name == 1) {gridPointValues = gridPointValues_hotspot * tau}
        else if (s_name == 2) {gridPointValues = gridPointValues_uniform}
        else if (s_name == 3) {gridPointValues = gridPointValues_cov_r}
        else if (s_name == 4) {gridPointValues = gridPointValues_cov_c_big}
        else {print("Incorrect input to start")}

        for (k in 1:8) {
            # We need a global test for each buffer width
            print(paste0(s_name, " ", k))
            match_vec = match_count_list[[k]]
            n_matches = match_vec[[1]]
            print(paste0("Match Num: ", n_matches))

            load(paste0("../Output_noWater/nullGridInfo/combinedMatchingSetup", k, ".dat"))
            load(paste0('../Output_noWater/origGridInfo/sim_orig_', k,".dat"))

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
            
            combinedMatchingSetupFix2 = combinedMatchingSetupFix
            combinedMatchingSetupFix2$DATA = combinedMatchingSetupFix$DATA[wRatioOk,]
            combinedMatchingSetupFix2$GRID_IND_1 = combinedMatchingSetupFix$GRID_IND_1[wRatioOk]
            combinedMatchingSetupFix2$GRID_IND_2 = combinedMatchingSetupFix$GRID_IND_2[wRatioOk]

            # =====================================================================

            v1 = sd(combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2, na.rm=TRUE)^2
            v2 = sd(combinedMatchingSetupFix2$DATA$ratioArea, na.rm=TRUE)^2

            for(ii in indexList_MAIN) {
                ## find matches
                area_temp = sim_orig$DATA$area1[ii] + sim_orig$DATA$area2[ii]
                ratio_temp = max(sim_orig$DATA$area1[ii] / sim_orig$DATA$area2[ii],
                                sim_orig$DATA$area2[ii] / sim_orig$DATA$area1[ii])
                dist_temp = sqrt(((area_temp - (combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2))^2/v1) +
                                    ((ratio_temp - combinedMatchingSetupFix2$DATA$ratioArea)^2 / v2))
                
                w50 = order(dist_temp)[1:n_matches]
                
                # Calculating the test statistics based on the surface
                tStats_temp = test_stats(gridPointValues, combinedMatchingSetupFix2, w50)
                null_dist = tStats_temp$t_stat

                global_null[[k]][ii,] = null_dist
            }

        }

        save(global_null, file = paste0("../Output_noWater/sim_results/Global/", save_type[s_name],
                                        "global_null_", trialNum, ".dat"))
    }

    # Step 2 -----------------------------------------------------------------------

    for (s_name in 1:4) {

        load(paste0("../Output_noWater/sim_results/Global/", save_type[s_name], "global_null_", trialNum, ".dat"))
        global_t_stat <- vector(mode = "list", length = 10)

        for(k in 1:8) {

            print(k)
            match_vec = match_count_list[[k]]
            n_matches = match_vec[[1]]
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
                global_t_stat[[k]][rep, 1] = mean(temp_max, na.rm = T)
                global_t_stat[[k]][rep, 2] = temp_loc[which.max(temp_max)]
            }

        }

        save(global_t_stat, file = paste0("../Output_noWater/sim_results/Global/", save_type[s_name],
                                        "global_t_stat_", trialNum, ".dat"))
    }

# }


