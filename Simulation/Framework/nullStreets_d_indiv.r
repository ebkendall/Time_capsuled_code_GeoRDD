library("sp")
library("sf")
library("rgeos")
library("raster")

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

    return(t_stat_df)
}

match_count <- c(160, 880, 480, 160, 1020, 140, 140, 300)

for(trialNum in 1:1000) {
    set.seed(trialNum)

    load("../Data/indexList_MAIN.RData")

    file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))

    tau = 0.5  

    perc_pval_match <- vector(mode = "list", length = 4) # Order: HotSpot, Uniform, Random, Correlated
    perc_pval_match[[1]] = perc_pval_match[[2]] = perc_pval_match[[3]] = perc_pval_match[[4]] = vector(mode = "list", length = 10)

    p_val_df <- vector(mode = "list", length = 4) # Order: HotSpot, Uniform, Random, Correlated
    p_val_df[[1]] = p_val_df[[2]] = p_val_df[[3]] = p_val_df[[4]] = vector(mode = "list", length = 10)


    for (k in 3:10) {

    print(k)
    j = match_count[k-2]
    print(paste0("Match Num: ", j))
    
    load(paste0('../Output_noWater/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output_noWater/origGridInfo/sim_orig_', k, '.dat'))

    for (s_name in 1:4) {

        load(file_names[s_name])
        gridPointValues = NULL

        if (s_name == 1) {gridPointValues = gridPointValues_hotspot * tau}
        else if (s_name == 2) {gridPointValues = gridPointValues_uniform}
        else if (s_name == 3) {gridPointValues = gridPointValues_cov_r}
        else if (s_name == 4) {gridPointValues = gridPointValues_cov_c_big}
        else {print("Incorrect input to start")}

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
        
        v1 = sd(combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2, na.rm=TRUE)^2
        v2 = sd(combinedMatchingSetupFix2$DATA$ratioArea, na.rm=TRUE)^2
        
        row_num = 1
        perc_pval_match[[s_name]][[k]] = data.frame("num_match" = j,
                                                    "perc_pval_less_05" = rep(NA, length(j)))
        p_val_df[[s_name]][[k]] = matrix(nrow = length(j), ncol = nrow(sim_orig$DATA))

        pval = rep(NA, nrow(sim_orig$DATA))

        for (ii in indexList_MAIN) {
            ## find matches
            area_temp = sim_orig$DATA$area1[ii] + sim_orig$DATA$area2[ii]
            ratio_temp = max(sim_orig$DATA$area1[ii] / sim_orig$DATA$area2[ii],
                            sim_orig$DATA$area2[ii] / sim_orig$DATA$area1[ii])
            dist_temp = sqrt(((area_temp - (combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2))^2/v1) +
                                ((ratio_temp - combinedMatchingSetupFix2$DATA$ratioArea)^2 / v2))

            w50 = order(dist_temp)[1:j]

            # Calculating the test statistics based on the surface
            tStats_temp = test_stats(gridPointValues, combinedMatchingSetupFix2, w50)
            null_dist = tStats_temp$t_stat

            orig_temp = test_stats_orig(gridPointValues, sim_orig, ii)
            stat_temp = orig_temp[1,1]

            pval[ii] = mean(null_dist > stat_temp)
        }

        perc_pval = mean(pval < 0.05, na.rm=TRUE)
        perc_pval_match[[s_name]][[k]]$perc_pval_less_05[row_num] = perc_pval
        p_val_df[[s_name]][[k]][row_num, ] = pval
        row_num = row_num + 1
    }
    }

    save(p_val_df, file = paste0("../Output_noWater/sim_results/p_vals_match_rel/p_val_df_", trialNum, ".dat"))
    save(perc_pval_match, file = paste0("../Output_noWater/sim_results/p_vals_match_rel/perc_pval_match_", trialNum, ".dat"))

}