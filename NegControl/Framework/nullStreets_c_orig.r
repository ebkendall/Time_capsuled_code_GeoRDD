# library(sp); library(sf); library(rgeos); library(raster)
library(spatstat)
library(sp)
library(rgeos)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_NEW.RData")
load("../Data/treesByPrec.RData")
load("../Data/streetsByPrec.RData")
Dir = '../Output_tree/origGridInfo/'
print(Dir)

adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)

# for (k in 1:length(adjust_val)) {
  k = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  print(k)
  set.seed(k)
  
  buff_ind = 10    

  sim_orig <- list(DATA = data.frame("area1" = rep(NA,164), "area2" = rep(NA,164), 
                                     "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                     "count1" = rep(NA,164), "count2" = rep(NA,164),
                                     "naive_pval" = rep(NA, 164),
                                     "int_1" = rep(NA, 164), "int_2" = rep(NA, 164), 
                                     "spatialDiff" = rep(NA, 164)))
    
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

        # Collecting Tree counts each precinct
        p1 = point.in.polygon(treesByPrec[[prec_ind_1]][,1], treesByPrec[[prec_ind_1]][,2],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        p2 = point.in.polygon(treesByPrec[[prec_ind_2]][,1], treesByPrec[[prec_ind_2]][,2],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
        count1 = sum(p1 > 0)
        count2 = sum(p2 > 0)
        
        s1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength1
        s2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength2
        
        n = count1 + count2
        p = 0.5
        pval = NA

        if (count1 <= n/2) {
            pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
        } else {
            pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
        }
        
        # Making the border line a spatstat object
        border_line_1_2 = totalStreetBuffInfo_NEW[[10]][[i]]$centerLine
        
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
        
        prec_1_x = treesByPrec[[prec_ind_1]][,1]
        prec_1_y = treesByPrec[[prec_ind_1]][,2]
        prec_2_x = treesByPrec[[prec_ind_2]][,1]
        prec_2_y = treesByPrec[[prec_ind_2]][,2]

        # Random assignment of center points ----------------------------------
        poly3 = gBuffer(border_line_1_2, width=buff_ind * 100)
        p3_1 = point.in.polygon(prec_1_x, prec_1_y,
                                poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
                                poly3@polygons[[1]]@Polygons[[1]]@coords[,2])
        p3_2 = point.in.polygon(prec_2_x, prec_2_y,
                                poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
                                poly3@polygons[[1]]@Polygons[[1]]@coords[,2])

        p3_in_1_x = prec_1_x[(p3_1 > 0) & (p1 == 0)]
        p3_in_1_y = prec_1_y[(p3_1 > 0) & (p1 == 0)]
        p3_in_2_x = prec_2_x[(p3_2 > 0) & (p2 == 0)]
        p3_in_2_y = prec_2_y[(p3_2 > 0) & (p2 == 0)]

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
        # ---------------------------------------------------------------------

        # Removing center points -----------------------------------------------
        # poly3 = gBuffer(border_line_1_2, width=buff_ind * 100)
        # p3_1 = point.in.polygon(prec_1_x, prec_1_y,
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,2])
        # p3_2 = point.in.polygon(prec_2_x, prec_2_y,
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,2])

        # prec_1_x = prec_1_x[(p3_1 == 0) | (p3_1 > 0 & p1 > 0)]
        # prec_1_y = prec_1_y[(p3_1 == 0) | (p3_1 > 0 & p1 > 0)]
        # prec_2_x = prec_2_x[(p3_2 == 0) | (p3_2 > 0 & p2 > 0)]
        # prec_2_y = prec_2_y[(p3_2 == 0) | (p3_2 > 0 & p2 > 0)]
        # ---------------------------------------------------------------------

        # Focus on points only in the box
        poly_box = matrix(c(box_x_min, box_y_max, 
                            box_x_min, box_y_min, 
                            box_x_max, box_y_min,
                            box_x_max, box_y_max), ncol = 2, byrow = T)
        poly_arr1 = point.in.polygon(prec_1_x, prec_1_y, poly_box[,1], poly_box[,2])
        poly_arr2 = point.in.polygon(prec_2_x, prec_2_y, poly_box[,1], poly_box[,2])

        prec_1_x = prec_1_x[poly_arr1 > 0]
        prec_1_y = prec_1_y[poly_arr1 > 0]
        prec_2_x = prec_2_x[poly_arr2 > 0]
        prec_2_y = prec_2_y[poly_arr2 > 0]

        # plot(poly_box[,1], poly_box[,2])
        # plot(prec1, add = T)
        # plot(prec2, add = T)
        # plot(poly1, add = T)
        # plot(poly2, add = T)
        # points(treesByPrec[[prec_ind_1]][,1], treesByPrec[[prec_ind_1]][,2])
        # points(treesByPrec[[prec_ind_2]][,1], treesByPrec[[prec_ind_2]][,2])
        # points(prec_2_x, prec_2_y, col = 'green')
        # points(prec_1_x, prec_1_y, col = 'red')
        # return(0)
        
        pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), 
                                    c(box_y_min, box_y_max))
        pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), 
                                    c(box_y_min, box_y_max))
        
        
        # Calculating the spatial component
        int_1 = density.ppp(pp_1, adjust = adjust_val[k], scalekernel = T)
        line_intensity_1 = int_1[b_line_1_2]
        int_line_1 = mean(line_intensity_1)
        
        int_2 = density.ppp(pp_2, adjust = adjust_val[k], scalekernel = T)
        line_intensity_2 = int_2[b_line_1_2]
        int_line_2 = mean(line_intensity_2)
        
        intensity_diff = abs(int_line_1 - int_line_2)

        sim_orig$DATA[i,] = c(area1, area2, s1, s2, count1, count2, pval,
                            int_line_1, int_line_2, intensity_diff)
        if(count1 == 0 & count2 == 0) print(paste0("No trees: ", k, ", ", i))
    }

  save(sim_orig, file = paste0(Dir, 'sim_orig_', k, '.dat'))
# }
