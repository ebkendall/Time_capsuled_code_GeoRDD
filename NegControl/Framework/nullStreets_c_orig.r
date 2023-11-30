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

for (k in 1:length(adjust_val)) {
  print(k)
  buff_ind = k+1    

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

    # Collecting Tree counts for Precinct 1 across both buffers
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
    
    poly3 = gBuffer(border_line_1_2, width=30)
    p3_1 = point.in.polygon(treesByPrec[[prec_ind_1]][,1], treesByPrec[[prec_ind_1]][,2],
                            poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
                            poly3@polygons[[1]]@Polygons[[1]]@coords[,2])
    p3_2 = point.in.polygon(treesByPrec[[prec_ind_2]][,1], treesByPrec[[prec_ind_2]][,2],
                            poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
                            poly3@polygons[[1]]@Polygons[[1]]@coords[,2])
    # Points
    prec_1_x = treesByPrec[[prec_ind_1]][p3_1 == 0,1]
    prec_1_y = treesByPrec[[prec_ind_1]][p3_1 == 0,2]
    prec_2_x = treesByPrec[[prec_ind_2]][p3_2 == 0,1]
    prec_2_y = treesByPrec[[prec_ind_2]][p3_2 == 0,2]
    
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
    
    pp_1 = ppp(prec_1_x, prec_1_y, 
               c(min(c(prec_1_x, prec_2_x, box_x_min)), max(c(prec_1_x, prec_2_x, box_x_max))), 
               c(min(c(prec_1_y, prec_2_y, box_y_min)), max(c(prec_1_y, prec_2_y, box_y_max))))
    pp_2 = ppp(prec_2_x, prec_2_y, 
               c(min(c(prec_1_x, prec_2_x, box_x_min)), max(c(prec_1_x, prec_2_x, box_x_max))), 
               c(min(c(prec_1_y, prec_2_y, box_y_min)), max(c(prec_1_y, prec_2_y, box_y_max))))
    
    
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
}
