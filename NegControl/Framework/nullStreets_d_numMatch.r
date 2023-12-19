set.seed(10)

match_count <- seq(20, 1200, by = 20)
load("../Data/indexList_MAIN.RData")

for (k in 1:8) {

    print(k)

    load(paste0("../Output_tree/nullGridInfo/combinedMatchingSetup", k, ".dat"))
    load(paste0("../Output_tree/origGridInfo/sim_orig_", k,".dat"))

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
    
    t_stat_streets = abs(combinedMatchingSetupFix2$count1 / combinedMatchingSetupFix2$streets1
                         - combinedMatchingSetupFix2$count2 / combinedMatchingSetupFix2$streets2)
    t_stat_streets_orig = abs(sim_orig$DATA$count1 / sim_orig$DATA$streets1
                              - sim_orig$DATA$count2 / sim_orig$DATA$streets2)

    row_num = 1
    perc_pval_match = data.frame("num_match" = match_count,
                                 "perc_pval_less_05" = rep(NA, length(match_count)))
    p_val_df = matrix(nrow = length(match_count), ncol = nrow(sim_orig$DATA))
    
    for(j in match_count) {
      print(j)
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
      perc_pval_match$perc_pval_less_05[row_num] = perc_pval
      p_val_df[row_num, ] = pval
      row_num = row_num + 1
    }
        
    save(p_val_df, file = paste0("../Output_tree/combination/p_val_df_street", k, ".dat"))
    save(perc_pval_match, file = paste0("../Output_tree/combination/perc_pval_match_street", k, ".dat"))
}

# Plotting the results
library(tidyverse, quietly = T)
library(gridExtra, quietly = T)
p = vector(mode = 'list', length = 8)
for(i in 1:12) {
    load(paste0('../Output_tree/combination/perc_pval_match_street', i, '.dat'))
    pval = perc_pval_match[1:60,]
    p[[i]] = ggplot(pval, aes( y=perc_pval_less_05, x=num_match)) + 
        geom_point(color = "red", size = 2) +
        geom_smooth(method = "loess", formula = y ~ x) +
        ggtitle(paste0("Matching's Effect on Type I Error (",i,"00 ft)")) +
        xlab("Number of Resampled Streets") + 
        ylab("Type I Error") +
        ylim(0,max(pval$perc_pval_less_05)) + 
        scale_x_continuous(breaks = pretty(pval$num_match, n = 10)) +
        geom_hline(yintercept=0.05, linetype="dashed", 
                   color = "black", size = 1.5) +
        theme(text = element_text(size=15))
    print(paste0("Min. p-value: ", min(pval$perc_pval_less_05)))
    print(which(pval$perc_pval_less_05 == min(pval$perc_pval_less_05)))
    
}
pdf("../_visualizations/Plots/numMatch_street.pdf", onefile = T)
grid.arrange(p[[1]], p[[2]], ncol = 1, nrow = 2)
grid.arrange(p[[3]], p[[4]], ncol = 1, nrow = 2)
grid.arrange(p[[5]], p[[6]], ncol = 1, nrow = 2)
grid.arrange(p[[7]], p[[8]], ncol = 1, nrow = 2)
grid.arrange(p[[9]], p[[10]], ncol = 1, nrow = 2)
grid.arrange(p[[11]], p[[12]], ncol = 1, nrow = 2)
dev.off()

# 300 is number of matches