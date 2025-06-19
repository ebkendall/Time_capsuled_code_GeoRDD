test_stats <- function(gridPointValues, combinedMatchingSetupFix, w50) {
    
    t_stat_df = data.frame("t_stat" = rep(0, length(w50)),
                           "naive_pval" = rep(0, length(w50)))
    
    rowInd = 1
    
    for (jj in w50) {
        
        if(rowInd %% 1000 == 0) print(rowInd)
        
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

trialNum = 1
set.seed(trialNum)

file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))
surface_type = c("Precinct Effect", "Constant", "Random", "Spatial")

tau = 0.5  

k = 2 # buff_width = 4

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
tot_off = combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2
rat_off = combinedMatchingSetupFix2$DATA$area1 / combinedMatchingSetupFix2$DATA$area2
rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]

# Split this up via the surface
t_stat_plot_list = vector(mode = 'list', length = length(file_names))
var_stat_plot_list = vector(mode = 'list', length = length(file_names))

ratio_breaks = seq(from = quantile(rat_off, probs = 0.10), to = quantile(rat_off, probs = 0.90), 
                   length.out = 15)
sum_breaks   = seq(from = quantile(tot_off, probs = 0.90), to = quantile(tot_off, probs = 0.10), 
                   length.out = 15)
c_name = NULL
r_name = NULL
for(i in 2:length(ratio_breaks)) r_name = c(r_name, paste0(round(ratio_breaks[i-1], digits = 3), 
                                                           "-", round(ratio_breaks[i], digits = 3)))
for(i in 2:length(sum_breaks)) c_name = c(c_name, paste0(round(sum_breaks[i-1], digits = 1),
                                                         "-", round(sum_breaks[i], digits = 1)))

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
    
    w1 = which(tot_off <= head(sum_breaks, 1) & tot_off >= tail(sum_breaks, 1))
    
    w2 = which(rat_off >= head(ratio_breaks, 1) & rat_off <= tail(ratio_breaks, 1))
    
    wAll = intersect(w1, w2)
    
    tStats_temp = test_stats(gridPointValues, combinedMatchingSetupFix2, wAll)
    testStatsNULL = tStats_temp$t_stat
    
    t_stat_new = log(testStatsNULL)
    
    plotting_info = data.frame("t_stat_new" = t_stat_new, "totLength" = tot_off[wAll], "ratioOff" = rat_off[wAll])
    plotting_info = plotting_info[!is.na(plotting_info$totLength), ]
    
    t_stat_plot = matrix(nrow = length(ratio_breaks)-1, ncol = length(sum_breaks)-1)
    var_stat_plot = matrix(nrow = length(ratio_breaks)-1, ncol = length(sum_breaks)-1)
    for(i in 2:length(ratio_breaks)) {
        sub_rat = plotting_info[plotting_info$ratioOff <= ratio_breaks[i] &
                                    plotting_info$ratioOff > ratio_breaks[i-1], , drop = F]
        for(j in 2:length(sum_breaks)) {
            sub_tot = sub_rat[sub_rat$totLength > sum_breaks[j] & sub_rat$totLength <= sum_breaks[j-1], ,drop=F]
            if(nrow(sub_tot) > 0) {
                t_stat_plot[i-1, j-1] = mean(sub_tot$t_stat_new, na.rm = T)
                var_stat_plot[i-1, j-1] = var(sub_tot$t_stat_new, na.rm = T)
            }
        }
    }
    
    colnames(t_stat_plot) = c_name
    rownames(t_stat_plot) = r_name
    
    colnames(var_stat_plot) = c_name
    rownames(var_stat_plot) = r_name
    
    t_stat_plot_list[[s_name]] = t_stat_plot
    var_stat_plot_list[[s_name]] = var_stat_plot
}

library(egg)

# Formatting for ggplot
x <- c_name
y <- r_name
t_stat_plot_list_gg = vector(mode = 'list', length = length(surface_type))
var_stat_plot_list_gg = vector(mode = 'list', length = length(surface_type))
for(k in 1:length(surface_type)) {
    data <- expand.grid(X=x, Y=y)
    data$Z <- c(t(t_stat_plot_list[[k]]))
    t_stat_plot_list_gg[[k]] = ggplot(data, aes(X, Y, fill= Z)) +
        geom_tile() +
        labs(title=paste0("Mean (Surface = ", surface_type[k], ")")) +
        scale_fill_gradient(low = "navy", high = "red", na.value="white") +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none",
              plot.title = element_text(size=10))
    
    data2 <- expand.grid(X=x, Y=y)
    data2$Z <- c(t(var_stat_plot_list[[k]]))
    var_stat_plot_list_gg[[k]] = ggplot(data2, aes(X, Y, fill= Z)) +
        geom_tile() +
        labs(title=paste0("Variance (Surface = ", surface_type[k], ")")) +
        scale_fill_gradient(low = "navy", high = "red", na.value="white") +
        scale_x_discrete(guide = guide_axis(angle = 90)) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none",
              plot.title = element_text(size=10))
}

app3 = ggarrange(t_stat_plot_list_gg[[1]] +
                     theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank()),
                 t_stat_plot_list_gg[[2]] +
                     theme(axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank()),
                 var_stat_plot_list_gg[[1]] +
                     theme(axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank()),
                 var_stat_plot_list_gg[[2]] +
                     theme(axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank()),
                 t_stat_plot_list_gg[[3]],
                 t_stat_plot_list_gg[[4]] +
                     theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank()),
                 var_stat_plot_list_gg[[3]] +
                     theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank()),
                 var_stat_plot_list_gg[[4]] +
                     theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank()),
                 nrow = 2, ncol = 4)


ggsave(filename = "../Plots/simulation_tau_appendix.pdf", plot = app3, width=11, height=5.67)




