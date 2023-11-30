# library(sp); library(sf); library(rgeos); library(raster)
library(spatstat)
library(sp)

load('../Data/nycSub.RData')
load("../Data/treesByPrec.RData")    # gridWithin_prec treesByPrec
load('../../RealData/Data/Street_Seg/streets10.dat') # longStrBroke

# Iterate through each buffer width
adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)

for (index in 1:length(adjust_val)) {
  
  for(k in 1:77) {
    
    set.seed(k)

    prec_num = nycSub$Precinct[k]

    print(paste0("Adjustment: ", adjust_val[index], ", Precinct: ", k))

    l = 1000 # arbitrary start length

    nullStr_point_data <- list(DATA = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                                 "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                                 "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                                 "count1" = rep(-1,l), "count2" = rep(-1,l), "naive_pval" = rep(-1,l),
                                                 "int_1" = rep(-1,l), "int_2" = rep(-1,l), "spatialDiff" = rep(-1,l)),
                              GRID_IND_1 = vector(mode = 'list', length = l),
                              GRID_IND_2 = vector(mode = 'list', length = l))
    rowNum = 1

    load(paste0("../../RealData/Data/OutputStrInfo_realData/strInfo_", 10, "_", k, ".dat")) # contains the buffer object

    print(paste0("Total length: ", length(streetLengthInfo_null)))
    for(i in 1:length(streetLengthInfo_null)) {
      print(paste0("i ", i))
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
                                  poly1@Polygons[[1]]@coords[,1], poly1@Polygons[[1]]@coords[,2])
            p2 = point.in.polygon(treesByPrec[[k]]$x, treesByPrec[[k]]$y,
                                  poly2@Polygons[[1]]@coords[,1], poly2@Polygons[[1]]@coords[,2])

            count1 = sum(p1 > 0)
            count2 = sum(p2 > 0)

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
            
            pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
            int_1 = density.ppp(pp_1, adjust = adjust_val[index], scalekernel = T)
            line_intensity_1 = int_1[b_line_1_2]
            int_line_1 = mean(line_intensity_1)
            
            pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
            int_2 = density.ppp(pp_2, adjust = adjust_val[index], scalekernel = T)
            line_intensity_2 = int_2[b_line_1_2]
            int_line_2 = mean(line_intensity_2)
            
            intensity_diff = abs(int_line_1 - int_line_2)
            
            
            
            nullStr_point_data$DATA[rowNum,] = c(k, i, j, s1, s2, area1, area2, 
                                                 T, count1, count2, pval,
                                                 int_line_1, int_line_2, intensity_diff)
            
            nullStr_point_data$GRID_IND_1[[rowNum]] = which(p1 > 0)
            nullStr_point_data$GRID_IND_2[[rowNum]] = which(p2 > 0)
            
            rowNum = rowNum + 1
          }
        }
      }
    }

    save(nullStr_point_data, file=paste0("../Output_tree/nullGridInfo/nullData", 
          index, "_", k,".dat", sep=''))
  }

}
