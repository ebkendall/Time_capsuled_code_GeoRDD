set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
load("../Output_tree_rev/match_count_list.dat")
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

# Theta ------------------------------------------------------------------------
indiv_results_theta = rep(NA, n_buff_width)
global_results_theta = matrix(NA, nrow = n_buff_width, ncol = 3)
colnames(global_results_theta) = c("max", "mean", "median")

# Tau --------------------------------------------------------------------------
indiv_results = matrix(NA, nrow = n_buff_width, ncol = length(adjust_val))
global_results = vector(mode = 'list', length = 3)
global_results[[1]] = global_results[[2]] = global_results[[3]] = 
    matrix(NA, nrow = n_buff_width, ncol = length(adjust_val))
names(global_results) = c("max", "mean", "median")

for(k in 1:n_buff_width) {
    print(k)
    m_theta = match_count_list[k,"theta"]
    m_tau = match_count_list[k,"tau"]
    
    load(paste0('../Output_tree_rev/nullGridInfo/combinedMatchingSetup', k, ".dat"))
    load(paste0('../Output_tree_rev/origGridInfo/origData_', k, '.dat'))
    
    ## Now remove data points where these ratios are much different
    area_ratio = c(na.omit(origData$str_info$area1 / origData$str_info$area2))
    area_ratio[area_ratio < 1] = 1 / area_ratio[area_ratio < 1]
    wMax_a = max(area_ratio)
    wMin_a = min(area_ratio)
    
    street_ratio = c(na.omit(origData$str_info$streets1 / origData$str_info$streets2))
    street_ratio[street_ratio < 1] = 1 / street_ratio[street_ratio < 1]
    wMax_s = max(street_ratio)
    wMin_s = min(street_ratio)
    
    # Which null streets meet these requirements
    wRatioOk = which(combinedMatchingSetupFix$DATA$ratioArea > wMin_a &
                         combinedMatchingSetupFix$DATA$ratioArea < wMax_a & 
                         combinedMatchingSetupFix$DATA$ratioStreet > wMin_s &
                         combinedMatchingSetupFix$DATA$ratioStreet < wMax_s)
    
    combinedMatchingSetupFix2 = combinedMatchingSetupFix$DATA[wRatioOk,]
    int_surface_info = combinedMatchingSetupFix$INT_SURFACE[wRatioOk,]
    
    # Test statistics: Theta
    t_stat = abs(combinedMatchingSetupFix2$count1 / combinedMatchingSetupFix2$streets1
                 - combinedMatchingSetupFix2$count2 / combinedMatchingSetupFix2$streets2)
    t_stat_orig = abs(origData$str_info$count1 / origData$str_info$streets1
                      - origData$str_info$count2 / origData$str_info$streets2)
    
    # Test statistics: Tau (NEW)
    t_stat_int_surface = int_surface_info[,3*(1:8)]
    t_stat_int_surface_orig = origData$str_surf$INT_SURFACE[,3*(1:8)]

    rat_off = combinedMatchingSetupFix2$streets1 / combinedMatchingSetupFix2$streets2
    rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
   
    n_null = nrow(combinedMatchingSetupFix2)
    n_matches = 2000
    
    # Theta --------------------------------------------------------------------
    store_theta = matrix(NA, nrow = length(indexList_MAIN), ncol = n_matches)
    indPvalue_theta = rep(NA, length(indexList_MAIN))
    globalPvalue_theta = rep(NA, 3)
    Y_theta = rep(NA, length(indexList_MAIN))
    X_theta = matrix(NA, length(indexList_MAIN), 2)
    
    counter = 1
    for(ii in indexList_MAIN) {
        
        off_temp = origData$str_info$streets1[ii] + origData$str_info$streets2[ii]
        ratio_temp = max(origData$str_info$streets1[ii] / origData$str_info$streets2[ii],
                         origData$str_info$streets2[ii] / origData$str_info$streets1[ii])
        
        Y_theta[counter] = t_stat_orig[ii]
        X_theta[counter,] = c(ratio_temp, off_temp)
        
        counter = counter + 1
    }
    
    for(ii in 1:nrow(X_theta)) {
        
        off_temp = X_theta[ii,2]
        ratio_temp = X_theta[ii,1]
        
        w1 = which((combinedMatchingSetupFix2$streets1 + 
                        combinedMatchingSetupFix2$streets2) > m_theta*off_temp &
                       (combinedMatchingSetupFix2$streets1 + 
                            combinedMatchingSetupFix2$streets2) < (1/m_theta)*off_temp)
        
        w2 = which(rat_off > m_theta*ratio_temp &
                       rat_off < (1/m_theta)*ratio_temp)
        
        wAll = intersect(w1, w2)
        
        if (length(wAll) > 10) {
            
            testStatsNULL = t_stat[wAll]
            
            testStatsNULL = testStatsNULL[which(testStatsNULL > 0)]
            
            store_theta[ii,] = sample(testStatsNULL, n_matches, replace=TRUE)
        }
        
    }
    
    for (jj in 1:nrow(X_theta)) {
        indPvalue_theta[jj] = mean(store_theta[jj,] > Y_theta[jj])
    }
    
    globalPvalue_theta[1] = mean(apply(store_theta, 2, max, na.rm=TRUE) > 
                                  max(Y_theta, na.rm=TRUE))
    globalPvalue_theta[2] = mean(apply(store_theta, 2, mean, na.rm=TRUE) > 
                                  mean(Y_theta, na.rm=TRUE))
    globalPvalue_theta[3] = mean(apply(store_theta, 2, median, na.rm=TRUE) > 
                                  median(Y_theta, na.rm=TRUE))
    
    indiv_results_theta[k] = mean(indPvalue_theta < .05, na.rm=TRUE)
    global_results_theta[k,] = globalPvalue_theta

    # Tau ----------------------------------------------------------------------
    store = array(NA, dim=c(length(indexList_MAIN), n_matches, length(adjust_val)))
    
    indPvalue = matrix(NA, length(indexList_MAIN), length(adjust_val))
    YtestSave = matrix(NA, length(indexList_MAIN), length(adjust_val))
    globalPvalue = matrix(NA, 3, length(adjust_val))
    
    for (kk in 1:length(adjust_val)) {
        
        Ytest = rep(NA, length(indexList_MAIN))
        Xtest = matrix(NA, length(indexList_MAIN), 2)
        
        counter = 1
        for(ii in indexList_MAIN) {
            
            off_temp = origData$str_info$streets1[ii] + origData$str_info$streets2[ii]
            ratio_temp = max(origData$str_info$streets1[ii] / origData$str_info$streets2[ii],
                            origData$str_info$streets2[ii] / origData$str_info$streets1[ii])
            
            Ytest[counter] = origData$str_surf$INT_SURFACE[ii,kk*3]
            Xtest[counter,] = c(ratio_temp, off_temp)
            
            counter = counter + 1
        }
        
        YtestSave[,kk] = Ytest
        
        for(ii in 1:nrow(Xtest)) {
            
            off_temp = Xtest[ii,2]
            ratio_temp = Xtest[ii,1]
            
            w1 = which((combinedMatchingSetupFix2$streets1 + 
                  combinedMatchingSetupFix2$streets2) > m_tau*off_temp &
                 (combinedMatchingSetupFix2$streets1 + 
                    combinedMatchingSetupFix2$streets2) < (1/m_tau)*off_temp)
    
            w2 = which(rat_off > m_tau*ratio_temp &
                        rat_off < (1/m_tau)*ratio_temp)
            
            wAll = intersect(w1, w2)
            
            if (length(wAll) > 10) {
                
                testStatsNULL = t_stat_int_surface[wAll, kk]
                
                testStatsNULL = testStatsNULL[which(testStatsNULL > 0)]
                
                store[ii,,kk] = sample(testStatsNULL, n_matches, replace=TRUE)
            }
            
        }
        
        for (jj in 1:nrow(Xtest)) {
            indPvalue[jj,kk] = mean(store[jj,,kk] > Ytest[jj])
        }
        
        globalPvalue[1,kk] = mean(apply(store[,,kk], 2, max, na.rm=TRUE) > 
                                    max(Ytest, na.rm=TRUE))
        globalPvalue[2,kk] = mean(apply(store[,,kk], 2, mean, na.rm=TRUE) > 
                                      mean(Ytest, na.rm=TRUE))
        globalPvalue[3,kk] = mean(apply(store[,,kk], 2, median, na.rm=TRUE) > 
                                      median(Ytest, na.rm=TRUE))
    }
    
    indiv_results[k,] = apply(indPvalue < .05, 2, mean, na.rm=TRUE)
    global_results[['max']][k, ] = globalPvalue[1,]
    global_results[['mean']][k, ] = globalPvalue[2,]
    global_results[['median']][k, ] = globalPvalue[3,]
}
print("Theta")
print(indiv_results_theta)
print(global_results_theta)
print("Tau")
print(indiv_results)
print(global_results)
