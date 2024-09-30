source('common_fnc.r')

options(warn=0)

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

same_sig = T

for (index in 1:8) {
    
  # for(k in 1:77) {
    k = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    set.seed(k)

    buff_ind = index + 2

    prec_num = nycSub$Precinct[k]

    print(paste0("Buffer Index: ", buff_ind, ", Precinct: ", k))

    l = 1000 # arbitrary start length

    if(same_sig) {
        null_str_info = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                   "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                   "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                   "count1" = rep(-1,l), "count2" = rep(-1,l), "naive_pval" = rep(-1,l),
                                   "tree_sigma" = rep(-1,l), "tree_flag" = rep(0,l))
    } else {
        null_str_info = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                   "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                   "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                   "count1" = rep(-1,l), "count2" = rep(-1,l), "naive_pval" = rep(-1,l),
                                   "tree_sigma1" = rep(-1,l), "tree_sigma2" = rep(-1,l),
                                   "tree_flag1" = rep(0,l), "tree_flag2" = rep(0,l))
    }

    null_str_surf = list(INT_SURFACE = data.frame("int_1a" = rep(-1,l), "int_1b" = rep(-1,l), "spatialDiff1" = rep(-1,l),
                                                  "int_2a" = rep(-1,l), "int_2b" = rep(-1,l), "spatialDiff2" = rep(-1,l),
                                                  "int_3a" = rep(-1,l), "int_3b" = rep(-1,l), "spatialDiff3" = rep(-1,l),
                                                  "int_4a" = rep(-1,l), "int_4b" = rep(-1,l), "spatialDiff4" = rep(-1,l),
                                                  "int_5a" = rep(-1,l), "int_5b" = rep(-1,l), "spatialDiff5" = rep(-1,l),
                                                  "int_6a" = rep(-1,l), "int_6b" = rep(-1,l), "spatialDiff6" = rep(-1,l),
                                                  "int_7a" = rep(-1,l), "int_7b" = rep(-1,l), "spatialDiff7" = rep(-1,l),
                                                  "int_8a" = rep(-1,l), "int_8b" = rep(-1,l), "spatialDiff8" = rep(-1,l)))
    rowNum = 1

    load(paste0("../../RealData/Data/OutputStrInfo_realData/strInfo_", buff_ind, "_", k, ".dat")) # contains the buffer object
    load(paste0('../../RealData/Data/Street_Seg/streets', buff_ind, '.dat')) # longStrBroke

    print(paste0("Total length: ", length(streetLengthInfo_null)))
    for(i in 1:length(streetLengthInfo_null)) {
      print(paste0("i ", i, ", length ", length(streetLengthInfo_null[[i]])))
      for(j in 1:length(streetLengthInfo_null[[i]])) {
        if(sum(is.na(streetLengthInfo_null[[i]][[j]])) == 0) {
          if(length(streetLengthInfo_null[[i]][[j]]$buffer@polygons) > 1){
            poly1 = streetLengthInfo_null[[i]][[j]]$buffer@polygons[[1]]
            poly2 = streetLengthInfo_null[[i]][[j]]$buffer@polygons[[2]]

            area1 = poly1@area
            area2 = poly2@area
            
            s1 = streetLengthInfo_null[[i]][[j]]$streetLength1
            s2 = streetLengthInfo_null[[i]][[j]]$streetLength2
            
            p1 = point.in.polygon(treesByPrec[[k]]$x, treesByPrec[[k]]$y,
                                  poly1@Polygons[[1]]@coords[,1], 
                                  poly1@Polygons[[1]]@coords[,2])
            p2 = point.in.polygon(treesByPrec[[k]]$x, treesByPrec[[k]]$y,
                                  poly2@Polygons[[1]]@coords[,1], 
                                  poly2@Polygons[[1]]@coords[,2])
            
            
            count1 = sum(p1 > 0)
            count2 = sum(p2 > 0)
            
            if(count1 < 2 | count2 < 2) next
            
            n = count1 + count2
            p = 0.5
            pval = NA
            
            if (count1 <= n/2) {
                pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
            } else {
                pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
            }
            
            # Making the border line a spatstat object
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
            prec_1_x = treesByPrec[[k]]$x[p1 > 0]
            prec_1_y = treesByPrec[[k]]$y[p1 > 0]
            
            prec_2_x = treesByPrec[[k]]$x[p2 > 0]
            prec_2_y = treesByPrec[[k]]$y[p2 > 0]

            # Random assignment of points in the buffer ------------------------
            poly3 = st_buffer(st_as_sf(border_line_1_2), dist = buff_ind * 100)
            poly3_coord = sf::st_coordinates(poly3)
            
            p3_1 = point.in.polygon(treesByPrec[[k]]$x, treesByPrec[[k]]$y,
                                    poly3_coord[,"X"], poly3_coord[,"Y"])
            
            prec_3_x_final = treesByPrec[[k]]$x[((p1 == 0) & (p2 == 0)) & (p3_1 > 0)]
            prec_3_y_final = treesByPrec[[k]]$y[((p1 == 0) & (p2 == 0)) & (p3_1 > 0)]
            
            assign_p3 = runif(n = length(prec_3_x_final))
            prec_3_x_1 = prec_3_x_final[assign_p3 > 0.5]
            prec_3_y_1 = prec_3_y_final[assign_p3 > 0.5]
            prec_3_x_2 = prec_3_x_final[assign_p3 <= 0.5]
            prec_3_y_2 = prec_3_y_final[assign_p3 <= 0.5]
            
            prec_1_x = c(prec_1_x, prec_3_x_1)
            prec_1_y = c(prec_1_y, prec_3_y_1)
            prec_2_x = c(prec_2_x, prec_3_x_2)
            prec_2_y = c(prec_2_y, prec_3_y_2)

            prec_both_x = c(prec_1_x, prec_2_x)
            prec_both_y = c(prec_1_y, prec_2_y)
            # ------------------------------------------------------------------
            
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

            # Skip the null street if:
            # (1) no observed trees, (2) only one tree between both sides
            if((pp_both_weight$n == 0) | (pp_both_weight$n == 1)) {next} 
            
            # par(mfrow = c(1,2))
            # plot(streetLengthInfo_null[[i]][[j]]$buffer)
            # points(pp_1_weight$x, pp_1_weight$y)
            # points(pp_2_weight$x, pp_2_weight$y, col = 'red')
            # plot(streetLengthInfo_null[[i]][[j]]$buffer)
            # points(pp_both_weight$x, pp_both_weight$y, col = 'orange')
            # dev.off()
            
            if(same_sig) {
                tree_surface_info = intensity_surf_same(pp_1_weight, pp_2_weight, pp_both_weight, 
                                                        b_line_1_2, adjust_val)
            } else {
                tree_surface_info = intensity_surf_diff(pp_1_weight, pp_2_weight, 
                                                        b_line_1_2, adjust_val)
            }
            
            
            null_str_info[rowNum,] = c(k, i, j, s1, s2, area1, area2, 
                                       T, count1, count2, pval,
                                       tree_surface_info$sig, 
                                       tree_surface_info$flag)

            null_str_surf$INT_SURFACE[rowNum,] = tree_surface_info$int_surf_vals
            
            rowNum = rowNum + 1
          }
        }
      }
    }

    nullData = list(str_info = null_str_info,
                    str_surf = null_str_surf)
    save(nullData, file=paste0("../Output_tree_rev/nullGridInfo/nullData", 
          index, "_", k,".dat", sep=''))

    rm(longStrBroke)
    rm(streetLengthInfo_null)
  # }
}
