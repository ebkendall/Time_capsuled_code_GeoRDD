library(sp)
library(spatstat)

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)
buff_val = 3:10
load("../../NegControl/Output_tree/combination/match_count_list.dat")
load("../Data/gridWithin_prec.rda")
load("../../NegControl/Data/totalStreetBuffInfo_NEW.RData")
load("../Data/ind_prec_df.rda")
load("../Data/nycSub.RData")

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

test_stats_surf <- function(gridPointValues, combinedMatchingSetupFix, 
                            w50, a_v, longStrBroke, gridWithin_prec) {
    
    t_stat_df = NULL
    
    for (jj in w50) {
        
        k = combinedMatchingSetupFix$DATA$precinct[jj]
        i = combinedMatchingSetupFix$DATA$indigo[jj]
        j = combinedMatchingSetupFix$DATA$juliet[jj]
        
        border_line_1_2 = longStrBroke[[k]][[i]][[j]]$shorterStreet
        b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
        b_line_pp = ppp(b_c_1_2[,"x"], b_c_1_2[,"y"],
                        c(min(b_c_1_2[,"x"]), max(b_c_1_2[,"x"])), 
                        c(min(b_c_1_2[,"y"]), max(b_c_1_2[,"y"])))
        b_line_1_2 = psp(b_c_1_2[1:(nrow(b_c_1_2)-1),"x"], 
                         b_c_1_2[1:(nrow(b_c_1_2)-1),"y"],
                         b_c_1_2[2:nrow(b_c_1_2),"x"],
                         b_c_1_2[2:nrow(b_c_1_2),"y"],
                         Window(b_line_pp))
        
        gridVals_ind_1 = combinedMatchingSetupFix$GRID_IND_1[[jj]]
        gridVals_ind_2 = combinedMatchingSetupFix$GRID_IND_2[[jj]]
        
        position1 = which(gridWithin_prec[[k]]$index %in% gridVals_ind_1)
        position2 = which(gridWithin_prec[[k]]$index %in% gridVals_ind_2)

        prec_1_x = gridWithin_prec[[k]]@coords[position1,1]
        prec_1_y = gridWithin_prec[[k]]@coords[position1,2]
        prec_2_x = gridWithin_prec[[k]]@coords[position2,1]
        prec_2_y = gridWithin_prec[[k]]@coords[position2,2]
        
        box_x_min = min(c(prec_1_x, prec_2_x))
        box_x_max = max(c(prec_1_x, prec_2_x))
        box_y_min = min(c(prec_1_y, prec_2_y))
        box_y_max = max(c(prec_1_y, prec_2_y))
        
        pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        
        gridValues1 = gridPointValues[gridVals_ind_1]
        gridValues2 = gridPointValues[gridVals_ind_2]
        
        pp_1_weight = pp_1 %mark% gridValues1
        pp_2_weight = pp_2 %mark% gridValues2
        
        int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks, 
                            sigma = bw.diggle, adjust = a_v,
                            scalekernel = T)
        line_intensity_1 = int_1[b_line_1_2]
        int_line_1 = mean(line_intensity_1)
        
        int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks, 
                            sigma = bw.diggle, adjust = a_v,
                            scalekernel = T)
        line_intensity_2 = int_2[b_line_1_2]
        int_line_2 = mean(line_intensity_2)  
        
        t_stat_df = c(t_stat_df, abs(int_line_1 - int_line_2))
        
        # plot(int_1)
        # plot(pp_1, add = T)
        # plot(b_line_1_2, col = 'green', add = T)
        # 
        # plot(int_2)
        # plot(pp_2, add = T)
        # plot(b_line_1_2, col = 'green', add = T)
    }
    
    return(t_stat_df)
}

test_stats_orig_surf <- function(gridPointValues, buff_ind, sim_orig, ii, a_v,
                                 totalStreetBuffInfo_NEW) {
    
    t_stat_df = NULL
    
    prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[ii])
    prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[ii])
    
    # Making the border line a spatstat object
    border_line_1_2 = totalStreetBuffInfo_NEW[[buff_ind]][[ii]]$centerLine
    b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
    
    # Make sure the line is actually linear
    endpts = which(multiplicity(b_c_1_2) == 1)
    if(length(endpts) > 2) {
        print(paste0(i, " problem!"))
    } 
    if(endpts[1] != 1 | endpts[2] != nrow(b_c_1_2)) {
        # print(paste0(i, " out of order"))
        b_c_1_2_new = rbind(b_c_1_2[endpts[2]:nrow(b_c_1_2), ], b_c_1_2[1:endpts[1], ])
        b_c_1_2 = b_c_1_2_new
    }
    
    b_line_pp = ppp(b_c_1_2[,"x"], b_c_1_2[,"y"],
                    c(min(b_c_1_2[,"x"]), max(b_c_1_2[,"x"])), 
                    c(min(b_c_1_2[,"y"]), max(b_c_1_2[,"y"])))
    b_line_1_2 = psp(b_c_1_2[1:(nrow(b_c_1_2)-1),"x"], 
                     b_c_1_2[1:(nrow(b_c_1_2)-1),"y"],
                     b_c_1_2[2:nrow(b_c_1_2),"x"],
                     b_c_1_2[2:nrow(b_c_1_2),"y"],
                     Window(b_line_pp))
    
    prec1 = nycSub[prec_ind_1, ]
    prec2 = nycSub[prec_ind_2, ]
    
    box_x_min = min(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,1], 
                      prec2@polygons[[1]]@Polygons[[1]]@coords[,1]))
    box_x_max = max(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,1], 
                      prec2@polygons[[1]]@Polygons[[1]]@coords[,1]))
    box_y_min = min(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,2], 
                      prec2@polygons[[1]]@Polygons[[1]]@coords[,2]))
    box_y_max = max(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,2], 
                      prec2@polygons[[1]]@Polygons[[1]]@coords[,2]))
    
    if(length(prec1@polygons[[1]]@Polygons) > 1) {
        for (pj in 1:length(prec1@polygons[[1]]@Polygons)) {
            box_x_min = min(c(box_x_min, prec1@polygons[[1]]@Polygons[[pj]]@coords[,1]))
            box_x_max = max(c(box_x_max, prec1@polygons[[1]]@Polygons[[pj]]@coords[,1]))
            box_y_min = min(c(box_y_min, prec1@polygons[[1]]@Polygons[[pj]]@coords[,2]))
            box_y_max = max(c(box_y_max, prec1@polygons[[1]]@Polygons[[pj]]@coords[,2]))
        }
    }
    
    if(length(prec2@polygons[[1]]@Polygons) > 1) {
        for (pj in 1:length(prec2@polygons[[1]]@Polygons)) {
            box_x_min = min(c(box_x_min, prec2@polygons[[1]]@Polygons[[pj]]@coords[,1]))
            box_x_max = max(c(box_x_max, prec2@polygons[[1]]@Polygons[[pj]]@coords[,1]))
            box_y_min = min(c(box_y_min, prec2@polygons[[1]]@Polygons[[pj]]@coords[,2]))
            box_y_max = max(c(box_y_max, prec2@polygons[[1]]@Polygons[[pj]]@coords[,2]))
        }
    }
    
    prec_1_x = gridWithin_prec[[prec_ind_1]]@coords[,1]
    prec_1_y = gridWithin_prec[[prec_ind_1]]@coords[,2]
    prec_2_x = gridWithin_prec[[prec_ind_2]]@coords[,1]
    prec_2_y = gridWithin_prec[[prec_ind_2]]@coords[,2]
    
    pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
    pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
    
    gridValues1 = gridPointValues[gridWithin_prec[[prec_ind_1]]$index]
    gridValues2 = gridPointValues[gridWithin_prec[[prec_ind_2]]$index]
    
    pp_1_weight = pp_1 %mark% gridValues1
    pp_2_weight = pp_2 %mark% gridValues2
    
    int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks, 
                        sigma = bw.diggle, adjust = a_v,
                        scalekernel = T)
    line_intensity_1 = int_1[b_line_1_2]
    int_line_1 = mean(line_intensity_1)
    
    int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks, 
                        sigma = bw.diggle, adjust = a_v,
                        scalekernel = T)
    line_intensity_2 = int_2[b_line_1_2]
    int_line_2 = mean(line_intensity_2)  
    
    t_stat_df = abs(int_line_1 - int_line_2)
    
    return(t_stat_df)
    
    # plot(int_1)
    # plot(b_line_1_2, col = 'green', add = T)
    # plot(int_2)
    # plot(b_line_1_2, col = 'green', add = T)
    # 
    # l <- solist(int_1, int_2)
    # plot(l, equal.ribbon = TRUE, main = "", 
    #      main.panel = c("prec 1", "prec 2"))
}

# for(trialNum in 1:1000) {
    trialNum = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    set.seed(trialNum)

    load("../Data/indexList_MAIN.RData")

    file_names <- c(paste0("../Data/Surfaces/gridPointValues_hotspot_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_uniform_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_r_", trialNum,".rda"),
                    paste0("../Data/Surfaces/gridPointValues_cov_c_", trialNum,".rda"))

    tau = 0.5  

    # Order: Precinct, Uniform, Random, Correlated
    perc_pval_match <- vector(mode = "list", length = 4)
    perc_pval_match[[1]] = perc_pval_match[[2]] = perc_pval_match[[3]] = 
            perc_pval_match[[4]] = vector(mode = "list", length = length(buff_val))
    # Order: Precinct, Uniform, Random, Correlated
    int_surface_pval = vector(mode = "list", length = 4)
    int_surface_pval[[1]] = int_surface_pval[[2]] = int_surface_pval[[3]] = 
            int_surface_pval[[4]] = vector(mode = "list", length = length(buff_val))
    # Order: Precinct, Uniform, Random, Correlated
    p_val_df <- vector(mode = "list", length = 4)
    p_val_df[[1]] = p_val_df[[2]] = p_val_df[[3]] = 
        p_val_df[[4]] = vector(mode = "list", length = length(buff_val))
    

    for (k in 1:length(buff_val)) {

        print(paste0("--> ", k))
        
        match_vec = match_count_list[[k]]
        j = match_vec[[1]]
        
        buff_ind = k + 2
        
        load(paste0('../Output_noWater/nullGridInfo/combinedMatchingSetup', k, ".dat"))
        load(paste0('../Output_noWater/origGridInfo/sim_orig_', k, '.dat'))
        load(paste0('../../RealData/Data/Street_Seg/streets', buff_ind, '.dat'))
    
        for (s_name in 1:4) {
            print(file_names[s_name])
            
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
            pval_int = matrix(NA, nrow = nrow(sim_orig$DATA), ncol = length(adjust_val))
    
            for (ii in indexList_MAIN) {
                print(paste0("ii:", ii))
                ## find matches
                area_temp = sim_orig$DATA$area1[ii] + sim_orig$DATA$area2[ii]
                ratio_temp = max(sim_orig$DATA$area1[ii] / sim_orig$DATA$area2[ii],
                                sim_orig$DATA$area2[ii] / sim_orig$DATA$area1[ii])
                dist_temp = sqrt(((area_temp - (combinedMatchingSetupFix2$DATA$area1 + combinedMatchingSetupFix2$DATA$area2))^2/v1) +
                                    ((ratio_temp - combinedMatchingSetupFix2$DATA$ratioArea)^2 / v2))
    
                w50 = order(dist_temp)[1:j]
    
                # Calculating the test statistics based on the surface
                print("OG null dist")
                tStats_temp = test_stats(gridPointValues, combinedMatchingSetupFix2, w50)
                null_dist = tStats_temp$t_stat
                print("OG obs test")
                orig_temp = test_stats_orig(gridPointValues, sim_orig, ii)
                stat_temp = orig_temp[1,1]
                
                pval[ii] = mean(null_dist > stat_temp)
                
                # Calculating the intensity surface info
                for(kk in 1:ncol(pval_int)) {
                    print(paste0("kk: ", kk))
                    w50_k = order(dist_temp)[1:match_vec[[2]][kk]]
                    
                    a_v = adjust_val[kk]
                    
                    null_dist = test_stats_surf(gridPointValues, 
                                                  combinedMatchingSetupFix2, 
                                                  w50_k, a_v, longStrBroke,
                                                  gridWithin_prec)
                    stat_temp = test_stats_orig_surf(gridPointValues,buff_ind,
                                                     sim_orig, ii, a_v,
                                                     totalStreetBuffInfo_NEW)
                    
                    pval_int[ii, kk] = mean(null_dist > stat_temp, na.rm=TRUE)
                }
            }
    
            int_surface_pval[[s_name]][[k]] = pval_int
            
            perc_pval = mean(pval < 0.05, na.rm=TRUE)
            perc_pval_match[[s_name]][[k]]$perc_pval_less_05[row_num] = perc_pval
            p_val_df[[s_name]][[k]][row_num, ] = pval
            row_num = row_num + 1
        }
        rm(longStrBroke)
    }

    save(p_val_df, file = paste0("../Output_noWater/sim_results/p_vals_match_rel/p_val_df_", trialNum, ".dat"))
    save(perc_pval_match, file = paste0("../Output_noWater/sim_results/p_vals_match_rel/perc_pval_match_", trialNum, ".dat"))
    save(int_surface_pval, file = paste0("../Output_noWater/sim_results/p_vals_match_rel/int_surface_pval_", trialNum, ".dat"))
# }