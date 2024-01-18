set.seed(2023)

load("../../NegControl/Output_tree/combination/match_count_list.dat")
load("../Data/indexList_MAIN.RData")

for(k in 1:8) {
    debug_invest = vector(mode = 'list', length = max(indexList_MAIN))
    
    kk = 2
    
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
    
    # Matching on offenses ---------------------------------------------------------
    v1 = sd(combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2, na.rm=TRUE)^2
    
    rat_off = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
    rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
    v2 = sd(rat_off, na.rm=TRUE)^2
    
    
    for (ii in indexList_MAIN) {
        debug_invest[[ii]] = vector(mode = 'list', length = 2)
        
        ## find matches
        off_temp = sim_orig$DATA$n_off_1_prec[ii] + sim_orig$DATA$n_off_2_prec[ii]
        ratio_temp = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
                         sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
        
        dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2))^2/v1) +
                             ((ratio_temp - rat_off)^2 / v2))
        
        debug_invest[[ii]][[1]] = matrix(nrow = match_vec[[2]][[kk]], ncol = 2)
        colnames(debug_invest[[ii]][[1]]) = c("ratioCrime", "ratioStreets")
        
        w50_k = order(dist_temp)[1:match_vec[[2]][kk]]
        debug_invest[[ii]][[1]][,"ratioCrime"] = rat_off[w50_k]
        debug_invest[[ii]][[1]][,"ratioStreets"] = combinedMatchingSetupFix2[w50_k, "ratioStreet"]
    }
    
    # Matching on streets ---------------------------------------------------------
    v1 = sd(combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2, na.rm=TRUE)^2
    v2 = sd(combinedMatchingSetupFix2$ratioStreet, na.rm=TRUE)^2
    
    for (ii in indexList_MAIN) {
        ## find matches
        off_temp = sim_orig$DATA$streets1[ii] + sim_orig$DATA$streets2[ii]
        ratio_temp = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                         sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])
        
        dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2))^2/v1) +
                             ((ratio_temp - combinedMatchingSetupFix2$ratioStreet)^2 / v2))
        
        debug_invest[[ii]][[2]] = matrix(nrow = match_vec[[2]][[kk]], ncol = 2)
        colnames(debug_invest[[ii]][[2]]) = c("ratioCrime", "ratioStreets")
        
        w50_k = order(dist_temp)[1:match_vec[[2]][kk]]
        debug_invest[[ii]][[2]][,"ratioCrime"] = rat_off[w50_k]
        debug_invest[[ii]][[2]][,"ratioStreets"] = combinedMatchingSetupFix2[w50_k, "ratioStreet"]
    }
    
    pval_updated = matrix(nrow = 164, ncol = 4)
    
    pdf(paste0("debugRatio", k, ".pdf"))
    par(mfrow=c(2,2))
    for(ii in indexList_MAIN) {
        ratio_street = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                           sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])
        ratio_crime = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
                          sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
        
        pval_updated[ii, 1] = mean(debug_invest[[ii]][[1]][,1] > ratio_crime, na.rm=TRUE)
        pval_updated[ii, 2] = mean(debug_invest[[ii]][[1]][,2] > ratio_street, na.rm=TRUE)
        pval_updated[ii, 3] = mean(debug_invest[[ii]][[2]][,1] > ratio_crime, na.rm=TRUE)
        pval_updated[ii, 4] = mean(debug_invest[[ii]][[2]][,2] > ratio_street, na.rm=TRUE)
    }
    hist(pval_updated[,4], breaks = ceiling(sqrt(nrow(pval_updated))),
         main = "p-vals match streets, ratio streets",
         xlab = paste0("prop < 0.05: ", round(mean(pval_updated[,4] < 0.05, na.rm = T), digits = 3)))
    
    hist(pval_updated[,3], breaks = ceiling(sqrt(nrow(pval_updated))),
         main = "p-vals match streets, ratio crime",
         xlab = paste0("prop < 0.05: ", round(mean(pval_updated[,3] < 0.05, na.rm = T), digits = 3)))
    
    hist(pval_updated[,1], breaks = ceiling(sqrt(nrow(pval_updated))),
         main = "p-vals match crime, ratio crime",
         xlab = paste0("prop < 0.05: ", round(mean(pval_updated[,1] < 0.05, na.rm = T), digits = 3)))
    
    hist(pval_updated[,2], breaks = ceiling(sqrt(nrow(pval_updated))),
         main = "p-vals match crime, ratio streets",
         xlab = paste0("prop < 0.05: ", round(mean(pval_updated[,2] < 0.05, na.rm = T), digits = 3)))
    
    for(ii in indexList_MAIN) {
        
        ratio_street = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                           sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])
        ratio_crime = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
                          sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
        # Match on streets, plot ratio of streets
        p4 = round(median(debug_invest[[ii]][[2]][,2], na.rm=TRUE), digits = 3)
        xlabel = paste0("median: ", p4, 
                        ", observed: ", round(ratio_street, digits = 3),
                        ", pval: ", round(mean(debug_invest[[ii]][[2]][,2] > ratio_street, 
                                               na.rm=TRUE), digits = 3))
        hist(debug_invest[[ii]][[2]][,2], 
             breaks = ceiling(sqrt(nrow(debug_invest[[ii]][[2]]))),
             main = paste0(ii, ", Match streets, ratio streets"),
             xlim = c(min(debug_invest[[ii]][[2]][,2], ratio_street), 
                      max(debug_invest[[ii]][[2]][,2], ratio_street)),
             xlab = xlabel)
        abline(v = ratio_street, col = 'red', lwd = 2)
        
        # Match on streets, plot ratio crime
        p3 = round(median(debug_invest[[ii]][[2]][,1], na.rm=TRUE), digits = 3)
        xlabel = paste0("median: ", p3, 
                        ", observed: ", round(ratio_crime, digits = 3),
                        ", pval: ", round(mean(debug_invest[[ii]][[2]][,1] > ratio_crime, 
                                               na.rm=TRUE), digits = 3))
        hist(debug_invest[[ii]][[2]][,1], 
             breaks = ceiling(sqrt(nrow(debug_invest[[ii]][[2]]))),
             main = paste0(ii, ", Match streets, ratio crime"),
             xlim = c(min(debug_invest[[ii]][[2]][,1], ratio_crime), 
                      max(debug_invest[[ii]][[2]][,1], ratio_crime)),
             xlab = xlabel)
        abline(v = ratio_crime, col = 'red', lwd = 2)
        
        # Match on offenses, plot ratio crime
        p1 = round(median(debug_invest[[ii]][[1]][,1], na.rm=TRUE), digits = 3)
        xlabel = paste0("median: ", p1, 
                        ", observed: ", round(ratio_crime, digits = 3),
                        ", pval: ", round(mean(debug_invest[[ii]][[1]][,1] > ratio_crime,
                                               na.rm=TRUE), digits = 3))
        hist(debug_invest[[ii]][[1]][,1], 
             breaks = ceiling(sqrt(nrow(debug_invest[[ii]][[1]]))),
             main = paste0(ii, ", Match crime, ratio crime"),
             xlim = c(min(debug_invest[[ii]][[1]][,1], ratio_crime), 
                      max(debug_invest[[ii]][[1]][,1], ratio_crime)),
             xlab = xlabel)
        abline(v = ratio_crime, col = 'red', lwd = 2)
        
        # Match on offenses, plot ratio of streets
        p2 = round(median(debug_invest[[ii]][[1]][,2], na.rm=TRUE), digits = 3)
        xlabel = paste0("median: ", p2, 
                        ", observed: ", round(ratio_street, digits = 3),
                        ", pval: ", round(mean(debug_invest[[ii]][[1]][,2] > ratio_street,
                                               na.rm=TRUE), digits = 3))
        hist(debug_invest[[ii]][[1]][,2], 
             breaks = ceiling(sqrt(nrow(debug_invest[[ii]][[1]]))),
             main = paste0(ii, ", Match crime, ratio streets"),
             xlim = c(min(debug_invest[[ii]][[1]][,2], ratio_street), 
                      max(debug_invest[[ii]][[1]][,2], ratio_street)),
             xlab = xlabel)
        abline(v = ratio_street, col = 'red', lwd = 2)
    }
    dev.off()
}

