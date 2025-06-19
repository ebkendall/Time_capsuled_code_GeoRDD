test_stats <- function(gridPointValues, combinedMatchingSetupFix, w50) {
    
    t_stat_df = data.frame("t_stat" = rep(0, length(w50)),
                           "naive_pval" = rep(0, length(w50)))
    
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

options(warn=1)
load("../Data/indexList_MAIN.RData")
n_buff_width = 8

# for(trialNum in 1:1000) {
    trialNum = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    set.seed(trialNum)
    
    file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))
    
    tau = 0.5  
    
    # Theta ------------------------------------------------------------------------
    indiv_results_theta = matrix(NA, nrow = n_buff_width, ncol = 4)
    
    g_mat = matrix(NA, nrow = n_buff_width, ncol = 2)
    colnames(g_mat) = c("max", "mean")
    
    global_results_theta = vector(mode = 'list', length  = 4)
    global_results_theta[[1]] = global_results_theta[[2]] = 
        global_results_theta[[3]] = global_results_theta[[4]] = g_mat
        
    
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
        rat_off = combinedMatchingSetupFix2$DATA$area1 / combinedMatchingSetupFix2$DATA$area2
        rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
        
        for(s_name in 1:4) {
            print(paste0("k ", k, " s ", s_name))
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
            
            n_matches = 2000
            
            # Theta --------------------------------------------------------------------
            indPvalue_theta = rep(NA, length(indexList_MAIN))
            globalPvalue_theta = rep(NA, 2)
            Y_theta = rep(NA, length(indexList_MAIN))
            X_theta = matrix(NA, length(indexList_MAIN), 2)
            
            counter = 1
            for(ii in indexList_MAIN) {
                
                area_temp = sim_orig$DATA$area1[ii] + sim_orig$DATA$area2[ii]
                ratio_temp = max(sim_orig$DATA$area1[ii] / sim_orig$DATA$area2[ii],
                                 sim_orig$DATA$area2[ii] / sim_orig$DATA$area1[ii])
                
                orig_temp = test_stats_orig(gridPointValues, sim_orig, ii)
                Y_theta[counter] = orig_temp[1,1]
                X_theta[counter,] = c(ratio_temp, area_temp)
                
                counter = counter + 1
            }
            
            m_theta = 0.99 # initialize measure of similarity
            store_theta = NULL
            repeat {
                total_match = rep(0, nrow(X_theta))
                store_theta = matrix(NA, nrow = length(indexList_MAIN), ncol = n_matches)
                ind_keep_theta = NULL
                
                for(ii in 1:nrow(X_theta)) {
                    off_temp = X_theta[ii,2]
                    ratio_temp = X_theta[ii,1]
                    
                    w1 = which((combinedMatchingSetupFix2$DATA$area1 + 
                                    combinedMatchingSetupFix2$DATA$area2) > m_theta*off_temp &
                                   (combinedMatchingSetupFix2$DATA$area1 + 
                                        combinedMatchingSetupFix2$DATA$area2) < (1/m_theta)*off_temp)
                    
                    w2 = which(rat_off > m_theta*ratio_temp &
                                   rat_off < (1/m_theta)*ratio_temp)
                    
                    wAll = intersect(w1, w2)
                    
                    total_match[ii] = length(wAll)
                    
                    if (length(wAll) > 10) {
                        tStats_temp = test_stats(gridPointValues, combinedMatchingSetupFix2, wAll)
                        testStatsNULL = tStats_temp$t_stat
                        testStatsNULL = testStatsNULL[which(testStatsNULL > 0)]
                        store_theta[ii,] = sample(testStatsNULL, n_matches, replace=TRUE)
                        ind_keep_theta = c(ind_keep_theta, ii)
                    } 
                }
                
                if(summary(total_match)[2] > 10) {
                    break
                } else {
                    print("Need more matches (theta)!")
                    print(paste0("Current perc = ", m_theta))
                    print(summary(total_match))
                    m_theta = m_theta - 0.01
                }
            }
            
            for (jj in 1:nrow(X_theta)) {
                indPvalue_theta[jj] = mean(store_theta[jj,] > Y_theta[jj])
            }
            
            globalPvalue_theta[1] = mean(apply(store_theta[ind_keep_theta,], 2, max, na.rm=TRUE) > 
                                             max(Y_theta[ind_keep_theta], na.rm=TRUE))
            globalPvalue_theta[2] = mean(apply(store_theta[ind_keep_theta,], 2, mean, na.rm=TRUE) > 
                                             mean(Y_theta[ind_keep_theta], na.rm=TRUE))
            
            indiv_results_theta[k, s_name] = mean(indPvalue_theta < .05, na.rm=TRUE)
            global_results_theta[[s_name]][k,] = globalPvalue_theta
            
        }
    }
    
    print("Theta")
    print(indiv_results_theta)
    print(global_results_theta)
    
    save(indiv_results_theta, file = paste0('../Output_noWater_rev/test_results/indiv_test_', trialNum, '.rda'))
    save(global_results_theta, file = paste0('../Output_noWater_rev/test_results/global_test_', trialNum, '.rda'))
# }

all_indiv = array(NA, dim=c(1000, n_buff_width, 4))

all_global = vector(mode = 'list', length = 4)
all_global[[1]] = all_global[[2]] =
    all_global[[3]] = all_global[[4]] = vector(mode = 'list', length = 8)

for(t in 1:1000) {

    load(paste0('../Output_noWater_rev/test_results/indiv_test_', t, '.rda'))
    all_indiv[t,,] = indiv_results_theta

    load(paste0('../Output_noWater_rev/test_results/global_test_', t, '.rda'))

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
                 round(perc_rej_j[3], 3), ", ",
                 round(perc_rej_j[4], 3), ") & (",
                 round(perc_rej_j[5], 3), ", ",
                 round(perc_rej_j[6], 3), ") & (",
                 round(perc_rej_j[7], 3), ", ",
                 round(perc_rej_j[8], 3), ")"))
}



# hist(Y_theta, breaks = 12, col = 'black', xlim = c(0, max(store_theta, na.rm = T)))
# hist(store_theta[,1], breaks = 12, add = T)
# hist(store_theta[,10], breaks = 12, col = 'blue', add = T)
# hist(store_theta[,100], breaks = 12, col = 'green', add = T)
# hist(store_theta[,1000], breaks = 12, col = 'red', add = T)
# 
# hist(apply(store_theta, 2, max, na.rm=TRUE), main = 'Histogram of sim w/ Max (grey) and Mean (red)',
#      xlim = c(0, max(store_theta, na.rm = T)), xlab = 'Test stats')
# hist(apply(store_theta, 2, mean, na.rm=TRUE), col = 'red', add = T)
# abline(v = max(Y_theta, na.rm=TRUE), lwd = 2, lty = 2)
# abline(v = mean(Y_theta, na.rm=TRUE), col = 'red', lwd = 2, lty = 2)
# 
# library(tidyverse)
# library(gridExtra, quietly = T)
# s_name = c( "Precinct", "Constant", "Random", "Spatial")
# p_max = p_mean = vector(mode = 'list', length = n_buff_width)
# for(i in 1:n_buff_width) {
#     p_max[[i]] = p_mean[[i]] = vector(mode = 'list', length = 4)
#     for(s in 1:4) {
#         adjPVal500 = data.frame("p" = all_global[[s]][[i]][,1])
#         p_max[[i]][[s]] = ggplot(adjPVal500, aes(x=p)) +
#             geom_histogram(color="black", fill="white", bins = sqrt(144)) +
#             xlab("P-Values") +
#             ylab("Frequency") +
#             labs(title = "Corrected p-values (max)",
#                  subtitle = paste0("k = ", i+2, ", ", s_name[s])) +
#             theme(text = element_text(size=8))
# 
#         adjPVal500 = data.frame("p" = all_global[[s]][[i]][,2])
#         p_mean[[i]][[s]] = ggplot(adjPVal500, aes(x=p)) +
#             geom_histogram(color="black", fill="white", bins = sqrt(144)) +
#             xlab("P-Values") +
#             ylab("Frequency") +
#             labs(title = "Corrected p-values (mean)",
#                  subtitle = paste0("k = ", i+2, ", ", s_name[s])) +
#             theme(text = element_text(size=8))
#     }
# }
# pdf("pHist.pdf", onefile = T)
# for(i in 1:n_buff_width) {
#     grid.arrange(p_max[[i]][[2]], p_mean[[i]][[2]], p_max[[i]][[3]], p_mean[[i]][[3]],
#                  p_max[[i]][[4]], p_mean[[i]][[4]], p_max[[i]][[1]], p_mean[[i]][[1]],
#                  nrow = 4, ncol = 2)
# }
# dev.off()

