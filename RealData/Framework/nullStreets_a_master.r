library(spatstat)
library(sp)
library(rgeos)

load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub
load('../Data/nycSub.RData')

# Iterate through each buffer width
adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)
buff_ind = 5

load(paste0('../Data/Street_Seg/streets', buff_ind, '.dat')) # longStrBroke

for (index in 1:length(adjust_val)) {
    
  # for(k in 1:77) {
    k = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    set.seed(k)

    prec_num = nycSub$Precinct[k]
    arr_sub = dataArr_sub[dataArr_sub$arrest_precinct == prec_num, ]
    off_sub = dataOff_sub[dataOff_sub$precinct == prec_num, ]

    print(paste0("Adjustment: ", adjust_val[index], ", Precinct: ", k))

    l = 1000 # arbitrary start length

    nullStr_point_data <- list(DATA = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                      "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                      "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                      "n_arr_1" = rep(-1,l), "n_arr_2" = rep(-1,l), "n_off_1" = rep(-1,l),
                                      "n_off_2" = rep(-1,l), "naive_pval" = rep(-1,l), 
                                      "int_1" = rep(-1,l), "int_2" = rep(-1,l), 
                                      "spatialDiff" = rep(-1,l)),
                              ARR_IND_1 = vector(mode = 'list', length = l),
                              ARR_IND_2 = vector(mode = 'list', length = l),
                              OFF_IND_1 = vector(mode = 'list', length = l),
                              OFF_IND_2 = vector(mode = 'list', length = l))
    rowNum = 1

    load(paste0("../Data/OutputStrInfo_realData/strInfo_", buff_ind, "_", k, ".dat")) # contains the buffer object

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
            
            n_arr_1 = sum(arr_1 > 0)
            n_arr_2 = sum(arr_2 > 0)
            n_off_1 = sum(off_1 > 0)
            n_off_2 = sum(off_2 > 0)
            
            arr_1_ind = arr_sub$main_ind[which(arr_1 > 0)]
            arr_2_ind = arr_sub$main_ind[which(arr_2 > 0)]
            
            off_1_ind = off_sub$main_ind[which(off_1 > 0)]
            off_2_ind = off_sub$main_ind[which(off_2 > 0)]
            
            # Naive p-value
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
            prec_1_x = arr_sub$x_coord_cd[arr_1 > 0]
            prec_1_y = arr_sub$y_coord_cd[arr_1 > 0]
            
            prec_2_x = arr_sub$x_coord_cd[arr_2 > 0]
            prec_2_y = arr_sub$y_coord_cd[arr_2 > 0]

            # Random assignment of points in the buffer ------------------------
            # poly3 = gBuffer(border_line_1_2, width=buff_ind * 100)
            # p3_1 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
            #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
            #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,2])

            # prec_3_x_final = arr_sub$x_coord_cd[((arr_1 == 0) & (arr_2 == 0)) & (p3_1 > 0)]
            # prec_3_y_final = arr_sub$y_coord_cd[((arr_1 == 0) & (arr_2 == 0)) & (p3_1 > 0)]
            
            # assign_p3 = runif(n = length(prec_3_x_final))
            # prec_3_x_1 = prec_3_x_final[assign_p3 > 0.5]
            # prec_3_y_1 = prec_3_y_final[assign_p3 > 0.5]
            # prec_3_x_2 = prec_3_x_final[assign_p3 <= 0.5]
            # prec_3_y_2 = prec_3_y_final[assign_p3 <= 0.5]
            
            # prec_1_x = c(prec_1_x, prec_3_x_1)
            # prec_1_y = c(prec_1_y, prec_3_y_1)
            # prec_2_x = c(prec_2_x, prec_3_x_2)
            # prec_2_y = c(prec_2_y, prec_3_y_2)
            
            # plot(streetLengthInfo_null[[i]][[j]]$buffer)
            # points(prec_1_x, prec_1_y)
            # points(prec_2_x, prec_2_y, col = 'red')
            # points(prec_3_x_2, prec_3_y_2, col = 'blue')
            # points(prec_3_x_1, prec_3_y_1, col = 'green')
            # ------------------------------------------------------------------
            
            pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
            int_1 = density.ppp(pp_1, adjust = adjust_val[index], scalekernel = T)
            line_intensity_1 = int_1[b_line_1_2]
            int_line_1 = mean(line_intensity_1)
            
            pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
            int_2 = density.ppp(pp_2, adjust = adjust_val[index], scalekernel = T)
            line_intensity_2 = int_2[b_line_1_2]
            int_line_2 = mean(line_intensity_2)
            
            intensity_diff = abs(int_line_1 - int_line_2)
            
            nullStr_point_data$DATA[rowNum,] = c(k, i, j,
                                                 streetLengthInfo_null[[i]][[j]]$streetLength1,
                                                 streetLengthInfo_null[[i]][[j]]$streetLength2,
                                                 area1, area2, T, 
                                                 n_arr_1, n_arr_2, 
                                                 n_off_1, n_off_2, pval,
                                                 int_line_1, int_line_2, intensity_diff)
            
            nullStr_point_data$ARR_IND_1[[rowNum]] = arr_1_ind
            nullStr_point_data$ARR_IND_2[[rowNum]] = arr_2_ind
            nullStr_point_data$OFF_IND_1[[rowNum]] = off_1_ind
            nullStr_point_data$OFF_IND_2[[rowNum]] = off_2_ind
            
            rowNum = rowNum + 1
          }
        }
      }
    }

    save(nullStr_point_data, file=paste0("../Output/nullGridInfo/nullData", 
          index, "_", k,".dat", sep=''))

  # }
}
