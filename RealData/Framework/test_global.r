set.seed(2024)
options(warn=1)
load("../Data/indexList_MAIN.RData")
n_matches = 400
n_buff_width = 8
adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

perc_rejections_indiv = rep(NA, n_buff_width)
perc_rejections_indiv_surf = matrix(nrow = n_buff_width, ncol = length(adjust_val))
rownames(perc_rejections_indiv_surf) = paste0(seq(from = 300, to = 1000, by = 100), "ft")
colnames(perc_rejections_indiv_surf) = adjust_val

p_global = matrix(nrow = n_buff_width, ncol = 3)
colnames(p_global) = c("max", "mean", "med")
rownames(p_global) = paste0(seq(from = 300, to = 1000, by = 100), "ft")
p_global_surf = vector(mode = 'list', length = 3)
p_global_surf[[1]] = p_global_surf[[2]] = 
    p_global_surf[[3]] = matrix(nrow = n_buff_width, ncol = length(adjust_val))

for(k in 1:n_buff_width) {
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
    
    # Null locations
    null_sum = combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2
    null_ratio = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
    null_ratio[null_ratio < 1] = 1 / null_ratio[null_ratio < 1]
    
    sd1 = sd(null_sum, na.rm = T)
    sd2 = sd(null_ratio, na.rm = T)
    # summary(null_sum); print(sd1)
    # summary(null_ratio); print(sd2)
    
    
    # Observed locations
    obs_sum = sim_orig$DATA$n_off_1_prec[indexList_MAIN] + sim_orig$DATA$n_off_2_prec[indexList_MAIN]
    obs_ratio = sim_orig$DATA$n_off_1_prec[indexList_MAIN] / sim_orig$DATA$n_off_2_prec[indexList_MAIN]
    obs_ratio[obs_ratio < 1] = 1 / obs_ratio[obs_ratio < 1]
    
    sd1_obs = sd(obs_sum, na.rm = T)
    sd2_obs = sd(obs_ratio, na.rm = T)
    # summary(obs_sum); print(sd1_obs)
    # summary(obs_ratio); print(sd2_obs)
    
    t_stat = abs(combinedMatchingSetupFix2$n_arr_1 / combinedMatchingSetupFix2$n_off_1
                 - combinedMatchingSetupFix2$n_arr_2 / combinedMatchingSetupFix2$n_off_2)
    t_stat_orig = abs(sim_orig$DATA$n_arr_1_prec / sim_orig$DATA$n_off_1_prec
                      - sim_orig$DATA$n_arr_2_prec / sim_orig$DATA$n_off_2_prec)
    
    t_stat_int_surface = int_surface_info[,1:8]
    t_stat_int_surface_orig = sim_orig$INT_SURFACE[,1:8]
    
    pval = rep(NA, nrow(sim_orig$DATA))
    pval_int = matrix(NA, nrow = nrow(sim_orig$DATA), ncol = length(adjust_val))
    
    global_null = matrix(nrow = 164, ncol = n_matches)
    global_null_int_surf= vector(mode = 'list', length = length(adjust_val))
    global_samp = vector(mode = 'list', length = 164)
    global_samp_int_surf = vector(mode = 'list', length = length(adjust_val))
    for(jj in 1:length(adjust_val)) {
        global_null_int_surf[[jj]] = matrix(nrow = 164, ncol = n_matches)
        global_samp_int_surf[[jj]] = vector(mode = 'list', length = 164)
    }
    
    null_str_position = combinedMatchingSetupFix2[,c("precinct", "indigo", "juliet")]
    matched_null_str_loc = vector(mode = 'list', length = 164)
    
    for(ii in indexList_MAIN) {
        ## find matches
        off_temp = sim_orig$DATA$n_off_1_prec[ii] + sim_orig$DATA$n_off_2_prec[ii]
        ratio_temp = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
                         sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
        
        # Remove relatively extreme values to improve the mahalanobis distance
        remove_extreme = which((null_ratio < 2*ratio_temp) & (null_ratio > (1/2)*ratio_temp))
        null_sum_ii = null_sum[remove_extreme]
        null_ratio_ii = null_ratio[remove_extreme]
        t_stat_ii = t_stat[remove_extreme]
        t_stat_int_surface_ii = t_stat_int_surface[remove_extreme, ]
        v1_ii = sd(null_sum_ii, na.rm = T)^2
        v2_ii = sd(null_ratio_ii, na.rm = T)^2
        null_str_position_ii = null_str_position[remove_extreme, ]
        
        
        dist_temp = sqrt(((off_temp - null_sum_ii)^2/v1_ii) + ((ratio_temp - null_ratio_ii)^2 / v2_ii))
        
        w50 = order(dist_temp)[1:n_matches]
        
        null_dist = t_stat_ii[w50]
        global_null[ii,] = null_dist
        
        ind_order = order(dist_temp)[1:min(5000, length(dist_temp))] # sample from the best 5,000
        t_stat_ind = t_stat_ii[ind_order]
        prob_w50 = exp(-dist_temp[ind_order])
        prob_w50 = prob_w50 / sum(prob_w50)
        prob_dist = cbind(ind_order, t_stat_ind, prob_w50)
        colnames(prob_dist) = c("w50_ind", "w50_tstat", "w50_prob")
        global_samp[[ii]] = prob_dist
        
        stat_temp = t_stat_orig[ii]
        pval[ii] = mean(null_dist > stat_temp, na.rm=TRUE)
        
        matched_null_str_loc[[ii]] = null_str_position_ii[w50, ]
        
        for(kk in 1:length(adjust_val)) {
            null_dist_int = t_stat_int_surface_ii[w50, ]
            global_null_int_surf[[kk]][ii, ] = null_dist_int[,kk]
            
            t_stat_ind_kk = t_stat_int_surface_ii[ind_order, kk]
            prob_dist_kk = cbind(ind_order, t_stat_ind_kk, prob_w50)
            colnames(prob_dist_kk) = c("w50_ind", "w50_tstat", "w50_prob")
            global_samp_int_surf[[kk]][[ii]] = prob_dist_kk
            
            pval_int[ii, kk] = mean(null_dist_int[,kk] > t_stat_int_surface_orig[ii,kk], na.rm=TRUE)
        }
    }
    # Individual test ---------------------------------------------------------
    print(paste0("Buffer: ", k+2, "00 ft Individual -------------------------"))
    print(paste0("Theta: ", round(mean(pval < 0.05, na.rm=TRUE), digits = 4)))
    print("Tau")
    tau_pval = apply(pval_int, 2, function(x){mean(x < 0.05, na.rm=TRUE)})
    print(round(tau_pval, digits = 4))    
    
    perc_rejections_indiv[k] = mean(pval < 0.05, na.rm=TRUE)
    perc_rejections_indiv_surf[k, ] = tau_pval
    
    
    # Global test -------------------------------------------------------------
    print(paste0("Buffer: ", k+2, "00 ft Global -----------------------------"))
    n_reps = 300
    global_t_stat = data.frame("max_t_stat"  = rep(NA, n_reps),
                               "mean_t_stat" = rep(NA, n_reps),
                               "med_t_stat" = rep(NA, n_reps),
                               "max_loc" = rep(NA, n_reps))
    global_t_stat_int_surf = vector(mode='list',length=length(adjust_val))
    for(jj in 1:length(global_t_stat_int_surf)) {
        global_t_stat_int_surf[[jj]] = data.frame("max_t_stat"  = rep(NA, n_reps),
                                                  "mean_t_stat" = rep(NA, n_reps),
                                                  "med_t_stat" = rep(NA, n_reps),
                                                  "max_loc" = rep(NA, n_reps))
    }
    
    # Theta
    for (rep in 1:n_reps) {
        # This is the repetition to get the null distribution
        temp_loc = temp_max = rep(NA, 164)
        for(ii in indexList_MAIN) {
            rand_ind = sample(c(1:n_matches), 1)
            temp_max[ii] = global_null[ii, rand_ind]
            temp_loc[ii] = rand_ind
            # temp_max[ii] = sample(x = global_samp[[ii]][,"w50_tstat"], size = 1,
            #                          prob = global_samp[[ii]][,"w50_prob"])
        }
        global_t_stat[rep, 1] = max(temp_max, na.rm = T)
        global_t_stat[rep, 2] = mean(temp_max, na.rm = T)
        global_t_stat[rep, 3] = median(temp_max, na.rm = T)
        global_t_stat[rep, 4] = temp_loc[which.max(temp_max)]
    }
    
    t_stat_max = max(t_stat_orig, na.rm = T)
    t_stat_mean = mean(t_stat_orig, na.rm = T)
    t_stat_median = median(t_stat_orig, na.rm = T)
    
    p_val_global = rep(NA, 3)
    
    p_val_global[1] = mean(global_t_stat$max_t_stat > t_stat_max)
    p_val_global[2] = mean(global_t_stat$mean_t_stat > t_stat_mean)
    p_val_global[3] = mean(global_t_stat$med_t_stat > t_stat_median)
    print("Theta")
    print(paste0("Max = ", round(p_val_global[1], digits = 4), 
                 ", Mean = ", round(p_val_global[2], digits = 4),
                 ", Med. = ", round(p_val_global[3], digits = 4)))
    
    p_global[k, ] = p_val_global
    
    # Tau
    global_match_loc = vector(mode = 'list', length = length(adjust_val))
    for(av in 1:length(adjust_val)) {
        global_match_loc[[av]] = matrix(0,nrow = n_reps, ncol = ncol(null_str_position))
        
        for (rep in 1:n_reps) {
            # This is the repetition to get the null distribution
            temp_loc = temp_max = rep(NA, 164)
            for(ii in indexList_MAIN) {
                rand_ind = sample(c(1:n_matches), 1)
                temp_max[ii] = global_null_int_surf[[av]][ii, rand_ind]
                temp_loc[ii] = rand_ind
                # temp_max[ii] = sample(x = global_samp_int_surf[[av]][[ii]][,"w50_tstat"], size = 1,
                #                          prob = global_samp_int_surf[[av]][[ii]][,"w50_prob"])
            }
            temp_loc_max = cbind(1:164, temp_loc, temp_max)
            temp_loc_max = temp_loc_max[indexList_MAIN, ]
            temp_loc_max = temp_loc_max[order(abs(temp_loc_max[,3] - median(temp_max, na.rm = T))), ]
            
            med_ii = as.numeric(temp_loc_max[1,1])
            match_ii = as.numeric(temp_loc_max[1,2])
            
            global_match_loc[[av]][rep, ] = as.numeric(matched_null_str_loc[[med_ii]][match_ii, ])
            
            global_t_stat_int_surf[[av]][rep, 1] = max(temp_max, na.rm = T)
            global_t_stat_int_surf[[av]][rep, 2] = mean(temp_max, na.rm = T)
            global_t_stat_int_surf[[av]][rep, 3] = median(temp_max, na.rm = T)
            global_t_stat_int_surf[[av]][rep, 4] = temp_loc[which.max(temp_max)]
        }
    }
    
    p_val_global_surf = matrix(nrow = 3, ncol = length(adjust_val))
    global_match_loc_orig = rep(NA, length(adjust_val))
    
    for(av in 1:length(adjust_val)) {
        for(m in 1:3) {
            if(m == 1) {
                t_stat_kk = max(t_stat_int_surface_orig[,av], na.rm = T)
                global_null_kk = global_t_stat_int_surf[[av]]$max_t_stat
            } else if(m == 2) {
                t_stat_kk = mean(t_stat_int_surface_orig[,av], na.rm = T)
                global_null_kk = global_t_stat_int_surf[[av]]$mean_t_stat
            } else {
                t_stat_kk = median(t_stat_int_surface_orig[,av], na.rm = T) 
                temp_orig = cbind(1:164, t_stat_int_surface_orig[,av])
                temp_orig = temp_orig[indexList_MAIN, ]
                temp_orig = temp_orig[order(abs(temp_orig[,2] - t_stat_kk)), ]
                global_match_loc_orig[av] = as.numeric(temp_orig[1,1])
                    
                global_null_kk = global_t_stat_int_surf[[av]]$med_t_stat
            }
            
            p_val_global_surf[m, av] = mean(global_null_kk > t_stat_kk, na.rm = T)
        }
    }
    
    colnames(p_val_global_surf) = adjust_val
    rownames(p_val_global_surf) = c("max", "mean", "med")
    print("Tau")
    print(round(p_val_global_surf, digits = 4))
    p_global_surf[[1]][k, ] = p_val_global_surf[1,]
    p_global_surf[[2]][k, ] = p_val_global_surf[2,]
    p_global_surf[[3]][k, ] = p_val_global_surf[3,]
}

print("Individual test results ---------------------------------------------")
print("Theta")
print(round(perc_rejections_indiv, digits = 4))
print("Tau")
print(round(perc_rejections_indiv_surf, digits = 4))

print("Global test results -------------------------------------------------")
print("Theta")
print(round(p_global, digits = 4))
print("Tau")
names(p_global_surf) = c("max", "mean", "med")
print(p_global_surf)

curr_results = list(perc_rejections_indiv = perc_rejections_indiv,
                    perc_rejections_indiv_surf = perc_rejections_indiv_surf,
                    p_global = p_global,
                    p_global_surf = p_global_surf)
# save(curr_results, file = 'curr_results_random.rda')


# What do these surfaces look like
load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub
load('../Data/nycSub.RData')
library(spatstat)
library(sp)
library(sf)

load(paste0('../Data/Street_Seg/streets', 4, '.dat')) # longStrBroke

for(a in c(2,4,8)) {
    pdf(paste0("null_str_global_", a, ".pdf"))
    par(mfrow = c(3,2))
    positions_a = global_match_loc[[a]]
    positions_a = positions_a[order(positions_a[,1]), ]
    unique_prec = unique(positions_a[,1])
    
    for(p in unique_prec) {
        print(p)
        load(paste0("../Data/OutputStrInfo_realData/strInfo_", 4, "_", p, ".dat")) # contains the buffer object  
        positions_a_sub = positions_a[positions_a[,1] == p, ,drop = F]
        
        prec_num = nycSub$Precinct[p]
        
        arr_sub = dataArr_sub[dataArr_sub$arrest_precinct == prec_num, ]
        off_sub = dataOff_sub[dataOff_sub$precinct == prec_num, ]
        
        for(ij in 1:nrow(positions_a_sub)) {
            i = positions_a_sub[ij, 2]
            j = positions_a_sub[ij, 3]
            
            poly1 = streetLengthInfo_null[[i]][[j]]$buffer@polygons[[1]]
            poly2 = streetLengthInfo_null[[i]][[j]]$buffer@polygons[[2]]
            
            arr_1 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
                                     poly1@Polygons[[1]]@coords[,1], 
                                     poly1@Polygons[[1]]@coords[,2])
            arr_2 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
                                     poly2@Polygons[[1]]@coords[,1], 
                                     poly2@Polygons[[1]]@coords[,2])
            
            off_1 = point.in.polygon(off_sub$x_coord_cd, off_sub$y_coord_cd,
                                     poly1@Polygons[[1]]@coords[,1], 
                                     poly1@Polygons[[1]]@coords[,2])
            off_2 = point.in.polygon(off_sub$x_coord_cd, off_sub$y_coord_cd,
                                     poly2@Polygons[[1]]@coords[,1], 
                                     poly2@Polygons[[1]]@coords[,2])
            
            arr_1_ind = arr_sub$main_ind[which(arr_1 > 0)]
            arr_2_ind = arr_sub$main_ind[which(arr_2 > 0)]
            
            off_1_ind = off_sub$main_ind[which(off_1 > 0)]
            off_2_ind = off_sub$main_ind[which(off_2 > 0)]
            
            # Making the border line a spatstat object
            border_line_1_2 = longStrBroke[[p]][[i]][[j]]$shorterStreet
            b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
            b_line_pp = ppp(b_c_1_2[,"x"], b_c_1_2[,"y"],
                            c(min(b_c_1_2[,"x"]), max(b_c_1_2[,"x"])), 
                            c(min(b_c_1_2[,"y"]), max(b_c_1_2[,"y"])))
            b_line_1_2 = psp(b_c_1_2[1:(nrow(b_c_1_2)-1),"x"], 
                             b_c_1_2[1:(nrow(b_c_1_2)-1),"y"],
                             b_c_1_2[2:nrow(b_c_1_2),"x"],
                             b_c_1_2[2:nrow(b_c_1_2),"y"],
                             Window(b_line_pp))
            
            # Defining window
            box_x_min = min(c(poly1@Polygons[[1]]@coords[,"x"], 
                              poly2@Polygons[[1]]@coords[,"x"]))
            box_x_max = max(c(poly1@Polygons[[1]]@coords[,"x"], 
                              poly2@Polygons[[1]]@coords[,"x"]))
            box_y_min = min(c(poly1@Polygons[[1]]@coords[,"y"], 
                              poly2@Polygons[[1]]@coords[,"y"]))
            box_y_max = max(c(poly1@Polygons[[1]]@coords[,"y"], 
                              poly2@Polygons[[1]]@coords[,"y"]))
            if(length(poly1@Polygons) > 1) {
                for (pj in 1:length(poly1@Polygons)) {
                    box_x_min = min(c(box_x_min, poly1@Polygons[[pj]]@coords[,"x"]))
                    box_x_max = max(c(box_x_max, poly1@Polygons[[pj]]@coords[,"x"]))
                    box_y_min = min(c(box_y_min, poly1@Polygons[[pj]]@coords[,"y"]))
                    box_y_max = max(c(box_y_max, poly1@Polygons[[pj]]@coords[,"y"]))
                }
            }
            
            if(length(poly2@Polygons) > 1) {
                for (pj in 1:length(poly2@Polygons)) {
                    box_x_min = min(c(box_x_min, poly2@Polygons[[pj]]@coords[,"x"]))
                    box_x_max = max(c(box_x_max, poly2@Polygons[[pj]]@coords[,"x"]))
                    box_y_min = min(c(box_y_min, poly2@Polygons[[pj]]@coords[,"y"]))
                    box_y_max = max(c(box_y_max, poly2@Polygons[[pj]]@coords[,"y"]))
                }
            }
            
            # Calculating the spatial component
            prec_1_x = arr_sub$x_coord_cd[arr_1 > 0]
            prec_1_y = arr_sub$y_coord_cd[arr_1 > 0]
            
            prec_2_x = arr_sub$x_coord_cd[arr_2 > 0]
            prec_2_y = arr_sub$y_coord_cd[arr_2 > 0]
            
            # Random assignment of points in the buffer ------------------------
            poly3 = st_buffer(st_as_sf(border_line_1_2), dist = buff_ind * 100)
            poly3_coord = sf::st_coordinates(poly3)
            # poly3 = gBuffer(border_line_1_2, width=buff_ind * 100)
            p3_1 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
                                    poly3_coord[,"X"], poly3_coord[,"Y"])
            
            prec_3_x_final = arr_sub$x_coord_cd[((arr_1 == 0) & (arr_2 == 0)) & (p3_1 > 0)]
            prec_3_y_final = arr_sub$y_coord_cd[((arr_1 == 0) & (arr_2 == 0)) & (p3_1 > 0)]
            
            assign_p3 = runif(n = length(prec_3_x_final))
            prec_3_x_1 = prec_3_x_final[assign_p3 > 0.5]
            prec_3_y_1 = prec_3_y_final[assign_p3 > 0.5]
            prec_3_x_2 = prec_3_x_final[assign_p3 <= 0.5]
            prec_3_y_2 = prec_3_y_final[assign_p3 <= 0.5]
            
            prec_1_x = c(prec_1_x, prec_3_x_1)
            prec_1_y = c(prec_1_y, prec_3_y_1)
            prec_2_x = c(prec_2_x, prec_3_x_2)
            prec_2_y = c(prec_2_y, prec_3_y_2)
            
            # plot(streetLengthInfo_null[[i]][[j]]$buffer)
            # points(prec_1_x, prec_1_y)
            # points(prec_2_x, prec_2_y, col = 'red')
            # points(prec_3_x_2, prec_3_y_2, col = 'blue')
            # points(prec_3_x_1, prec_3_y_1, col = 'green')
            # ------------------------------------------------------------------
            pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
            pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
            
            # Assigning weights
            single_1 <- !duplicated(pp_1)
            m1 <- multiplicity(pp_1)
            pp_1_weight <- pp_1[single_1] %mark% m1[single_1]
            
            single_2 <- !duplicated(pp_2)
            m2 <- multiplicity(pp_2)
            pp_2_weight <- pp_2[single_2] %mark% m2[single_2]
            
            if(length(pp_1_weight$x) == 1 | length(pp_2_weight$x) == 1) next
            
            if(pp_1_weight$n > 0) {
                if(pp_1_weight$n == 1) {
                    # Cannot do cross validation with one location
                    print(paste0("One data point (side 1, scale ", adjust_val[a], "): ", i))
                    int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks,
                                        adjust = adjust_val[a], scalekernel = T)
                } else {
                    int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks,
                                        sigma = bw.diggle, adjust = adjust_val[a],
                                        scalekernel = T)   
                }
                line_intensity_1 = int_1[b_line_1_2]
                int_line_1 = mean(line_intensity_1) 
            } else {
                int_line_1 = 0
            }
            
            if(pp_2_weight$n > 0) {
                if(pp_2_weight$n == 1) {
                    # Cannot do cross validation with one location
                    if(pp_1_weight$n > 1) {
                        # Take the spatial smoothness of the other side
                        sig_2 = attr(int_1, "sigma") # This has already been scaled
                        int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks, 
                                            sigma = sig_2, scalekernel = T)   
                    } else {
                        print(paste0("One data point (side 1, scale ", adjust_val[a], "): ", i))
                        int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks, 
                                            adjust = adjust_val[a], scalekernel = T)
                    }
                    
                } else {
                    int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks,
                                        sigma = bw.diggle, adjust = adjust_val[a],
                                        scalekernel = T)
                }
                line_intensity_2 = int_2[b_line_1_2]
                int_line_2 = mean(line_intensity_2)
            } else {
                int_line_2 = 0
            }
            
            # Double check, if pp_1_weight$n == 1 and pp_2_weight$n > 1, use sigma from surface 2
            if((pp_1_weight$n == 1) & (pp_2_weight$n > 1)) {
                print(paste0("Informed surface 1 from surface 2"))
                # If the other surface has information, use the same scaled-sigma
                sig_1 = attr(int_2, "sigma") # This has already been scaled
                int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks, 
                                    sigma = sig_1, scalekernel = T)
                line_intensity_1 = int_1[b_line_1_2]
                int_line_1 = mean(line_intensity_1)
            }
            
            plot(int_1, main = paste("Side = 1, Prec = ", p, ", i = ", i, ", j = ", j, '\n',
                                   "adjust = ", adjust_val[a], ", sigma = ", attr(int_1, "sigma")))
            plot(streetLengthInfo_null[[i]][[j]]$buffer, lwd = 2, add = T)
            plot(int_2, main = paste("Side = 2, Prec = ", p, ", i = ", i, ", j = ", j, '\n',
                                     "adjust = ", adjust_val[a], ", sigma = ", attr(int_2, "sigma")))
            plot(streetLengthInfo_null[[i]][[j]]$buffer, lwd = 2, add = T)
            
        }
    }
    
    dev.off()
    
}



