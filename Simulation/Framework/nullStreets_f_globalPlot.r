load("../Data/indexList_MAIN.RData")

test_stats_orig <- function(gridPointValues, sim_orig) {
    
    t_stat_df = matrix(nrow = 164, ncol = 2)
    colnames(t_stat_df) = c("t_stat", "naive_pval")
    
    for(ii in indexList_MAIN) {
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
        
        t_stat_df[ii, ] = c(abs(count1 / area1 - count2 / area2), pval)
        
    }
    
    return(t_stat_df)
}

save_type = c("HotSpot/", "Uniform/", "Random/", "Correlated/") 

# Step 3 -----------------------------------------------------------------------
p_val_df = vector(mode = "list", length = 8)

for(i in 1:8) {p_val_df[[i]] = data.frame( "HotSpot"     = rep(NA,1000),
                                           "Uniform"     = rep(NA,1000),
                                           "Random"      = rep(NA,1000),
                                           "Correlated"  = rep(NA,1000))}

whichMaxInfo = vector(mode = "list", length = 4)
whichMaxInfo[[1]] = whichMaxInfo[[2]] = whichMaxInfo[[3]] = whichMaxInfo[[4]] = vector(mode = 'list', length = 10)

for (trialNum in 1:1000) {
    set.seed(trialNum)
    print(trialNum)
    
    file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))
    
    tau = 0.5  
    for (s_name in 1:4) {
        
        load(paste0('../Output_noWater/sim_results/Global/', save_type[s_name], 
                    "global_t_stat_", trialNum,".dat"))
        
        load(file_names[s_name])
        gridPointValues = NULL
        
        if (s_name == 1) {gridPointValues = gridPointValues_hotspot * tau}
        else if (s_name == 2) {gridPointValues = gridPointValues_uniform}
        else if (s_name == 3) {gridPointValues = gridPointValues_cov_r}
        else if (s_name == 4) {gridPointValues = gridPointValues_cov_c_big}
        else {print("Incorrect input to start")}
        
        for(k in 1:8) {
            
            load(paste0('../Output_noWater/origGridInfo/sim_orig_', k,".dat"))
            
            t_stat_list = test_stats_orig(gridPointValues, sim_orig)
            
            t_stat = max(t_stat_list[,1], na.rm = T)
            w_max = which.max(na.omit(t_stat_list[,1]))
            
            whichMaxInfo[[s_name]][[k]] = rbind(whichMaxInfo[[s_name]][[k]], 
                                                c(w_max, t_stat))
            
            p_val = mean(global_t_stat[[k]]$max_t_stat > t_stat)
            p_val_df[[k]][trialNum, s_name] = p_val
        }
        
    }
}

save(p_val_df, file = paste0("../Output_noWater/global_p_values.rda"))

display_name = c("Precinct", "Constant", "Random", "Spatial")

print("Global Results")
pdf(paste0("../_visualizations/Plots/global_new_total.pdf"))
par(mfrow=c(2,2))
for (i in 1:8) {
    pval_save = rep(0, 4)
    for(k in 1:4) {
        print(paste0(i, " ", k))
        pval = p_val_df[[i]][,k]
        hist(pval, breaks = sqrt(length(pval)), main = paste0(display_name[k], ": pVal for B", i*100),
             xlab = paste0("Perc. < 0.05 is ",  round(mean(pval < 0.05, na.rm=TRUE), 4)),
             xlim=c(0,1))
        pval_save[k] = round(mean(pval < 0.05, na.rm=TRUE), 4)
    }
    print(paste0((i+2)*100, " & ", pval_save[2], " & ", pval_save[3], " & ",
                                   pval_save[4], " & ", pval_save[1]))
}

for (i in 1:8) {
    for(k in 1:4) {
        plot(table(indexList_MAIN[whichMaxInfo[[k]][[i]][,1]]), 
             ylab = "Freq", main = paste0("Max Obs TStat index: ", 
                                          names(which.max(table(indexList_MAIN[whichMaxInfo[[k]][[i]][,1]]))),
                                          "\n B:", (i+2)*100))
    }
}
dev.off()

# Individual plot
print("Individual Results")

load("../Output_noWater/sim_results/p_vals_match_rel/p_val_df_1.dat")
final_hist = p_val_df

for(j in 1:4) {
  for (k in 1:8) {
    final_hist[[j]][[k]] = final_hist[[j]][[k]][1,]
  }
}

for (i in c(2:1000)) {
    load(paste0("../Output_noWater/sim_results/p_vals_match_rel/p_val_df_", i, ".dat"))
    for(j in 1:4) {
        for(k in 1:8) {
            final_hist[[j]][[k]] = c(final_hist[[j]][[k]], p_val_df[[j]][[k]][1,])
        }
    }
}

folder_type = c("HotSpot", "Uniform", "Random", "Correlated")
display_name = c("Precinct", "Constant", "Random", "Spatial")

pdf("../_visualizations/Plots/indiv_new_adj.pdf")
par(mfrow=c(2,2))
for (i in 1:8) {
    pval_save = rep(0, 4)
    for(k in 1:4) {
        pval = final_hist[[k]][[i]]
        hist(pval, main = paste0(display_name[k], ": pVal for B", (i+2)*100),
            xlab = paste0("Perc. < 0.05 is ",  round(mean(pval < 0.05, na.rm=TRUE), 4)),
            xlim=c(0,1))
        pval_save[k] = round(mean(pval < 0.05, na.rm=TRUE), 4)
    }
    print(paste0((i+2)*100, " & ", pval_save[2], " & ", pval_save[3], " & ",
                                   pval_save[4], " & ", pval_save[1]))
}
dev.off()