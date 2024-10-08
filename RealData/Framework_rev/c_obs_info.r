source('common_fnc.r')

options(warn=0)

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

same_sig = T

# for (k in 1:length(adjust_val)) {
    k = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    print(k)
    set.seed(k)
    
    buff_ind = k + 2
    
    if(same_sig) {
        orig_str_info = data.frame( "area1" = rep(NA,164), "area2" = rep(NA,164), 
                                    "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                    "n_arr_1" = rep(NA, 164), "n_arr_2" = rep(NA, 164),
                                    "n_off_1" = rep(NA, 164), "n_off_2" = rep(NA, 164),
                                    "naive_pval" = rep(NA, 164),
                                    "n_arr_1_prec" = rep(NA, 164), "n_arr_2_prec" = rep(NA, 164),
                                    "n_off_1_prec" = rep(NA, 164), "n_off_2_prec" = rep(NA, 164),
                                    "naive_pval_prec" = rep(NA, 164),
                                    "arr_sigma" = rep(NA, 164), "off_sigma" = rep(NA, 164), 
                                    "arr_flag" = rep(NA, 164), "off_flag" = rep(NA, 164))
    } else {
        orig_str_info = data.frame( "area1" = rep(NA,164), "area2" = rep(NA,164), 
                                    "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                    "n_arr_1" = rep(NA, 164), "n_arr_2" = rep(NA, 164),
                                    "n_off_1" = rep(NA, 164), "n_off_2" = rep(NA, 164),
                                    "naive_pval" = rep(NA, 164),
                                    "n_arr_1_prec" = rep(NA, 164), "n_arr_2_prec" = rep(NA, 164),
                                    "n_off_1_prec" = rep(NA, 164), "n_off_2_prec" = rep(NA, 164),
                                    "naive_pval_prec" = rep(NA, 164),
                                    "arr_sigma1" = rep(NA, 164), "arr_sigma2" = rep(NA, 164),
                                    "off_sigma1" = rep(NA, 164), "off_sigma2" = rep(NA, 164), 
                                    "arr_flag1" = rep(NA, 164), "arr_flag2" = rep(NA, 164),
                                    "off_flag1" = rep(NA, 164), "off_flag2" = rep(NA, 164))
    }
    
    orig_str_surf = list(INT_SURFACE = data.frame("int_1a" = rep(NA, 164), "int_1b" = rep(NA, 164), "spatialDiff1" = rep(NA, 164),
                                                  "int_2a" = rep(NA, 164), "int_2b" = rep(NA, 164), "spatialDiff2" = rep(NA, 164),
                                                  "int_3a" = rep(NA, 164), "int_3b" = rep(NA, 164), "spatialDiff3" = rep(NA, 164),
                                                  "int_4a" = rep(NA, 164), "int_4b" = rep(NA, 164), "spatialDiff4" = rep(NA, 164),
                                                  "int_5a" = rep(NA, 164), "int_5b" = rep(NA, 164), "spatialDiff5" = rep(NA, 164),
                                                  "int_6a" = rep(NA, 164), "int_6b" = rep(NA, 164), "spatialDiff6" = rep(NA, 164),
                                                  "int_7a" = rep(NA, 164), "int_7b" = rep(NA, 164), "spatialDiff7" = rep(NA, 164),
                                                  "int_8a" = rep(NA, 164), "int_8b" = rep(NA, 164), "spatialDiff8" = rep(NA, 164)),
                         OFF_SURFACE = data.frame("int_1a" = rep(NA, 164), "int_1b" = rep(NA, 164), "spatialDiff1" = rep(NA, 164),
                                                  "int_2a" = rep(NA, 164), "int_2b" = rep(NA, 164), "spatialDiff2" = rep(NA, 164),
                                                  "int_3a" = rep(NA, 164), "int_3b" = rep(NA, 164), "spatialDiff3" = rep(NA, 164),
                                                  "int_4a" = rep(NA, 164), "int_4b" = rep(NA, 164), "spatialDiff4" = rep(NA, 164),
                                                  "int_5a" = rep(NA, 164), "int_5b" = rep(NA, 164), "spatialDiff5" = rep(NA, 164),
                                                  "int_6a" = rep(NA, 164), "int_6b" = rep(NA, 164), "spatialDiff6" = rep(NA, 164),
                                                  "int_7a" = rep(NA, 164), "int_7b" = rep(NA, 164), "spatialDiff7" = rep(NA, 164),
                                                  "int_8a" = rep(NA, 164), "int_8b" = rep(NA, 164), "spatialDiff8" = rep(NA, 164)))
    
    for (i in indexList_MAIN) {
        print(i)
        prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[i])
        prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[i])
        
        poly1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly1
        poly2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly2
        
        poly_ind1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly_ind1
        poly_ind2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly_ind2
        
        area1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$area1
        area2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$area2
        
        s1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength1
        s2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength2
        
        # Making the border line object ---------------------------------------
        border_line_1_2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$centerLine
        b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
        
        endpts = which(multiplicity(b_c_1_2) == 1)
        if(length(endpts) > 2) {
            print(paste0(i, " problem!"))
        } 
        if(endpts[1] != 1 | endpts[2] != nrow(b_c_1_2)) {
            print(paste0(i, " out of order"))
            b_c_1_2_new = rbind(b_c_1_2[endpts[2]:nrow(b_c_1_2), ], b_c_1_2[1:endpts[1], ])
            b_c_1_2 = b_c_1_2_new
            border_line_1_2@lines[[1]]@Lines[[1]]@coords = b_c_1_2
        }
        
        test_b_c = unique(b_c_1_2)
        b_line_pp = ppp(test_b_c[,"x"], test_b_c[,"y"],
                        c(min(test_b_c[,"x"]), max(test_b_c[,"x"])),
                        c(min(test_b_c[,"y"]), max(test_b_c[,"y"])))
        b_line_1_2 = psp(test_b_c[1:(nrow(test_b_c)-1),"x"],
                         test_b_c[1:(nrow(test_b_c)-1),"y"],
                         test_b_c[2:nrow(test_b_c),"x"],
                         test_b_c[2:nrow(test_b_c),"y"],
                         Window(b_line_pp))
        
        # Create a buffer
        poly3 = st_buffer(st_as_sf(border_line_1_2), dist = buff_ind * 100)
        poly3_coord = sf::st_coordinates(poly3)
        
        # Naive p-value info  -------------------------------------------------
        # Arrest Info (total)
        arr_1 = point.in.polygon(dataArr_sub[,"x_coord_cd"], dataArr_sub[,"y_coord_cd"],
                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        arr_2 = point.in.polygon(dataArr_sub[,"x_coord_cd"], dataArr_sub[,"y_coord_cd"],
                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
        # Offense Info (total)
        off_1 = point.in.polygon(dataOff_sub[,"x_coord_cd"], dataOff_sub[,"y_coord_cd"],
                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        off_2 = point.in.polygon(dataOff_sub[,"x_coord_cd"], dataOff_sub[,"y_coord_cd"],
                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
        # Counting the arrests and offenses
        n_arr_1 = sum(arr_1 > 0)
        n_arr_2 = sum(arr_2 > 0)
        n_off_1 = sum(off_1 > 0)
        n_off_2 = sum(off_2 > 0)
        
        count1 = n_arr_1
        count2 = n_arr_2
        n = count1 + count2
        p = 0.5
        pval = NA
        
        if (count1 <= n/2) {
            pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
        } else {
            pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
        }
        
        # Spatial intensity surface construction ------------------------------
        
        # Define window based on the buffer width
        prec1 = poly1
        prec2 = poly2
        
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
        
        poly_box = matrix(c(box_x_min, box_y_max, 
                            box_x_min, box_y_min, 
                            box_x_max, box_y_min,
                            box_x_max, box_y_max), ncol = 2, byrow = T)
        
        # Arrest and Offense Info (precinct specific)
        dataArr_prec_1 = dataArr_sub[dataArr_sub$arrest_precinct == ind_prec_df$prec1[i], ]
        dataArr_prec_2 = dataArr_sub[dataArr_sub$arrest_precinct == ind_prec_df$prec2[i], ]
        dataOff_prec_1_2 = dataOff_sub[dataOff_sub$precinct %in% c(ind_prec_df$prec1[i], ind_prec_df$prec2[i]), ]
        
        # Arrest Data ------------------------------
        arr_1_tot = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
                                     poly3_coord[,"X"], poly3_coord[,"Y"])
        arr_2_tot = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
                                     poly3_coord[,"X"], poly3_coord[,"Y"])
        prec_1_x = dataArr_prec_1$x_coord_cd[arr_1_tot > 0]
        prec_1_y = dataArr_prec_1$y_coord_cd[arr_1_tot > 0]
        prec_2_x = dataArr_prec_2$x_coord_cd[arr_2_tot > 0]
        prec_2_y = dataArr_prec_2$y_coord_cd[arr_2_tot > 0]
        
        poly_arr1 = point.in.polygon(prec_1_x, prec_1_y, poly_box[,1], poly_box[,2])
        poly_arr2 = point.in.polygon(prec_2_x, prec_2_y, poly_box[,1], poly_box[,2])
        
        prec_1_x = prec_1_x[poly_arr1 > 0]
        prec_1_y = prec_1_y[poly_arr1 > 0]
        prec_2_x = prec_2_x[poly_arr2 > 0]
        prec_2_y = prec_2_y[poly_arr2 > 0]
        prec_both_x = c(prec_1_x, prec_2_x)
        prec_both_y = c(prec_1_y, prec_2_y)
        
        pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        pp_both = ppp(prec_both_x, prec_both_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        
        single_1 <- !duplicated(pp_1)
        m1 <- multiplicity(pp_1)
        pp_1_weight <- pp_1[single_1] %mark% m1[single_1]
        
        single_2 <- !duplicated(pp_2)
        m2 <- multiplicity(pp_2)
        pp_2_weight <- pp_2[single_2] %mark% m2[single_2]
        
        single_both <- !duplicated(pp_both)
        m_both <- multiplicity(pp_both)
        pp_both_weight <- pp_both[single_both] %mark% m_both[single_both]
        
        # Offense Data -----------------------------
        off_1_tot = point.in.polygon(dataOff_prec_1_2[,"x_coord_cd"], dataOff_prec_1_2[,"y_coord_cd"],
                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        off_2_tot = point.in.polygon(dataOff_prec_1_2[,"x_coord_cd"], dataOff_prec_1_2[,"y_coord_cd"],
                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
        prec_1_x_off = dataOff_prec_1_2$x_coord_cd[off_1_tot > 0]
        prec_1_y_off = dataOff_prec_1_2$y_coord_cd[off_1_tot > 0]
        prec_2_x_off = dataOff_prec_1_2$x_coord_cd[off_2_tot > 0]
        prec_2_y_off = dataOff_prec_1_2$y_coord_cd[off_2_tot > 0]
        
        # Randomly assign the points on the border
        p3 = point.in.polygon(dataOff_prec_1_2[,"x_coord_cd"], dataOff_prec_1_2[,"y_coord_cd"],
                              poly3_coord[,"X"], poly3_coord[,"Y"])
        prec_3_x_0 = dataOff_prec_1_2$x_coord_cd[(off_1_tot == 0) & (off_2_tot == 0) & (p3>0)]
        prec_3_y_0 = dataOff_prec_1_2$y_coord_cd[(off_1_tot == 0) & (off_2_tot == 0) & (p3>0)]
        
        assign_p3 = runif(n = length(prec_3_x_0))
        prec_3_x_1 = prec_3_x_0[assign_p3 > 0.5]
        prec_3_y_1 = prec_3_y_0[assign_p3 > 0.5]
        prec_3_x_2 = prec_3_x_0[assign_p3 <= 0.5]
        prec_3_y_2 = prec_3_y_0[assign_p3 <= 0.5]
        
        prec_1_x_off = c(prec_1_x_off, prec_3_x_1)
        prec_1_y_off = c(prec_1_y_off, prec_3_y_1)
        prec_2_x_off = c(prec_2_x_off, prec_3_x_2)
        prec_2_y_off = c(prec_2_y_off, prec_3_y_2)
        
        poly_off1 = point.in.polygon(prec_1_x_off, prec_1_y_off, poly_box[,1], poly_box[,2])
        poly_off2 = point.in.polygon(prec_2_x_off, prec_2_y_off, poly_box[,1], poly_box[,2])
        
        prec_1_x_off = prec_1_x_off[poly_off1 > 0]
        prec_1_y_off = prec_1_y_off[poly_off1 > 0]
        prec_2_x_off = prec_2_x_off[poly_off2 > 0]
        prec_2_y_off = prec_2_y_off[poly_off2 > 0]
        prec_both_x_o = c(prec_1_x_off, prec_2_x_off)
        prec_both_y_o = c(prec_1_y_off, prec_2_y_off)
        
        pp_1_o = ppp(prec_1_x_off, prec_1_y_off, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        pp_2_o = ppp(prec_2_x_off, prec_2_y_off, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        pp_both_o = ppp(prec_both_x_o, prec_both_y_o, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        
        single_1_o <- !duplicated(pp_1_o)
        m1_o <- multiplicity(pp_1_o)
        pp_1_weight_o <- pp_1_o[single_1_o] %mark% m1_o[single_1_o]
        
        single_2_o <- !duplicated(pp_2_o)
        m2_o <- multiplicity(pp_2_o)
        pp_2_weight_o <- pp_2_o[single_2_o] %mark% m2_o[single_2_o]
        
        single_both_o <- !duplicated(pp_both_o)
        m_both_o <- multiplicity(pp_both_o)
        pp_both_weight_o <- pp_both_o[single_both_o] %mark% m_both_o[single_both_o]
        
        # par(mfrow=c(1,2))
        # plot(poly_box[,1], poly_box[,2])
        # plot(prec1, add = T)
        # plot(prec2, add = T)
        # points(pp_1_weight$x, pp_1_weight$y, col = 'red')
        # points(pp_2_weight$x, pp_2_weight$y, col = 'green')
        # points(pp_both_weight$x, pp_both_weight$y, col = 'orange')
        # plot(poly_box[,1], poly_box[,2])
        # plot(prec1, add = T)
        # plot(prec2, add = T)
        # points(pp_1_weight_o$x, pp_1_weight_o$y, col = 'red')
        # points(pp_2_weight_o$x, pp_2_weight_o$y, col = 'green')
        # points(pp_both_weight_o$x, pp_both_weight_o$y, col = 'orange')
        
        if(same_sig) {
            arrest_surface_info = intensity_surf_same(pp_1_weight, pp_2_weight, pp_both_weight, 
                                                      b_line_1_2, adjust_val)
            crime_surface_info  = intensity_surf_same(pp_1_weight_o, pp_2_weight_o, pp_both_weight_o, 
                                                      b_line_1_2, adjust_val)
        } else {
            arrest_surface_info = intensity_surf_diff(pp_1_weight, pp_2_weight, 
                                                      b_line_1_2, adjust_val)
            crime_surface_info  = intensity_surf_diff(pp_1_weight_o, pp_2_weight_o, 
                                                      b_line_1_2, adjust_val)
        }
        
        # Naive p-value (precinct specific) info  -----------------------------
        # arr_1_prec_a = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
        #                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
        #                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        # arr_1_prec_b = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
        #                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
        #                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        # arr_2_prec_a = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
        #                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
        #                                 poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        # arr_2_prec_b = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
        #                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
        #                                 poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
        n_arr_1_prec = length(prec_1_x) # sum(arr_1_prec_a > 0) + sum(arr_1_prec_b > 0) 
        n_arr_2_prec = length(prec_2_x) # sum(arr_2_prec_a > 0) + sum(arr_2_prec_b > 0) 
        n_off_1_prec = length(prec_1_x_off) # sum(off_1_tot > 0) 
        n_off_2_prec = length(prec_2_x_off) # sum(off_2_tot > 0)
        
        count1 = n_arr_1_prec
        count2 = n_arr_2_prec
        n = count1 + count2
        p = 0.5
        pval_prec = NA
        
        if (count1 <= n/2) {
            pval_prec = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
        } else {
            pval_prec = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
        }
        # ----------------------------------------------------------------------
        
        orig_str_info[i,] = c(area1, area2, s1, s2, n_arr_1, n_arr_2, n_off_1, n_off_2, pval,
                              n_arr_1_prec, n_arr_2_prec, n_off_1_prec, n_off_2_prec, pval_prec,
                              arrest_surface_info$sig, crime_surface_info$sig, 
                              arrest_surface_info$flag,crime_surface_info$flag)
        
        orig_str_surf$INT_SURFACE[i,] = arrest_surface_info$int_surf_vals
        orig_str_surf$OFF_SURFACE[i,] = crime_surface_info$int_surf_vals
        
    }
    print(warnings())
    
    origData = list(str_info = orig_str_info,
                    str_surf = orig_str_surf)
    save(origData, file=paste0("../Output_rev/origGridInfo/origData_", k,".dat"))
# }
