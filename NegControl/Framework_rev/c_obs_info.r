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
                                    "count1" = rep(NA,164), "count2" = rep(NA,164),
                                    "naive_pval" = rep(NA, 164),
                                    "tree_sigma" = rep(NA, 164), "tree_flag" = rep(NA, 164))
    } else {
        orig_str_info = data.frame( "area1" = rep(NA,164), "area2" = rep(NA,164), 
                                    "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                    "count1" = rep(NA,164), "count2" = rep(NA,164),
                                    "naive_pval" = rep(NA, 164),
                                    "tree_sigma1" = rep(NA, 164), "tree_sigma2" = rep(NA, 164),
                                    "tree_flag1" = rep(NA, 164), "tree_flag2" = rep(NA, 164))
    }
    
    orig_str_surf = list(INT_SURFACE = data.frame("int_1a" = rep(NA, 164), "int_1b" = rep(NA, 164), "spatialDiff1" = rep(NA, 164),
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
        
        p1 = point.in.polygon(treesByPrec[[prec_ind_1]][,1], treesByPrec[[prec_ind_1]][,2],
                              poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                              poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        p2 = point.in.polygon(treesByPrec[[prec_ind_2]][,1], treesByPrec[[prec_ind_2]][,2],
                              poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                              poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
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
        
        # Tree Data ------------------------------------------------------------
        prec_1_x = treesByPrec[[prec_ind_1]][p1 > 0,1]
        prec_1_y = treesByPrec[[prec_ind_1]][p1 > 0,2]
        prec_2_x = treesByPrec[[prec_ind_2]][p2 > 0,1]
        prec_2_y = treesByPrec[[prec_ind_2]][p2 > 0,2]
        
        # Random assignment of center points ----------------------------------
        poly3 = st_buffer(st_as_sf(border_line_1_2), dist = buff_ind * 100)
        poly3_coord = sf::st_coordinates(poly3)
        p3_1 = point.in.polygon(treesByPrec[[prec_ind_1]][,1], treesByPrec[[prec_ind_1]][,2],
                                poly3_coord[,"X"], poly3_coord[,"Y"])
        p3_2 = point.in.polygon(treesByPrec[[prec_ind_2]][,1], treesByPrec[[prec_ind_2]][,2],
                                poly3_coord[,"X"], poly3_coord[,"Y"])
        
        p3_in_1_x = treesByPrec[[prec_ind_1]][(p3_1 > 0) & (p1 == 0), 1]
        p3_in_1_y = treesByPrec[[prec_ind_1]][(p3_1 > 0) & (p1 == 0), 2]
        p3_in_2_x = treesByPrec[[prec_ind_2]][(p3_2 > 0) & (p2 == 0), 1]
        p3_in_2_y = treesByPrec[[prec_ind_2]][(p3_2 > 0) & (p2 == 0), 2]
        
        prec_3_x_final = c(p3_in_1_x, p3_in_2_x)
        prec_3_y_final = c(p3_in_1_y, p3_in_2_y)
        
        assign_p3 = runif(n = length(prec_3_x_final))
        prec_3_x_1 = prec_3_x_final[assign_p3 > 0.5]
        prec_3_y_1 = prec_3_y_final[assign_p3 > 0.5]
        prec_3_x_2 = prec_3_x_final[assign_p3 <= 0.5]
        prec_3_y_2 = prec_3_y_final[assign_p3 <= 0.5]
        
        prec_1_x = c(prec_1_x, prec_3_x_1)
        prec_1_y = c(prec_1_y, prec_3_y_1)
        prec_2_x = c(prec_2_x, prec_3_x_2)
        prec_2_y = c(prec_2_y, prec_3_y_2)
        
        # Trees are within window ----------------------------------------------
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
        
        # par(mfrow=c(1,2))
        # plot(poly_box[,1], poly_box[,2])
        # plot(prec1, add = T)
        # plot(prec2, add = T)
        # points(pp_1_weight$x, pp_1_weight$y, col = 'red')
        # points(pp_2_weight$x, pp_2_weight$y, col = 'green')
        # plot(poly_box[,1], poly_box[,2])
        # plot(prec1, add = T)
        # plot(prec2, add = T)
        # points(pp_both_weight$x, pp_both_weight$y, col = 'orange')
        
        if(same_sig) {
            tree_surface_info = intensity_surf_same(pp_1_weight, pp_2_weight, pp_both_weight, 
                                                    b_line_1_2, adjust_val)
        } else {
            tree_surface_info = intensity_surf_diff(pp_1_weight, pp_2_weight, 
                                                    b_line_1_2, adjust_val)
        }

        # Naive p-value info  -------------------------------------------------
        count1 = length(prec_1_x) # sum(p1 > 0)
        count2 = length(prec_2_x) # sum(p2 > 0)
        
        n = count1 + count2
        p = 0.5
        pval = NA
        
        if (count1 <= n/2) {
            pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
        } else {
            pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
        }
        # ---------------------------------------------------------------------
        
        orig_str_info[i,] = c(area1, area2, s1, s2, count1, count2, pval,
                              tree_surface_info$sig, tree_surface_info$flag)
        
        orig_str_surf$INT_SURFACE[i,] = tree_surface_info$int_surf_vals
        
    }
    
    origData = list(str_info = orig_str_info,
                    str_surf = orig_str_surf)
    save(origData, file=paste0("../Output_tree_rev/origGridInfo/origData_", k,".dat"))
# }
