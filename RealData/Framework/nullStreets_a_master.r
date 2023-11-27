library(spatstat)
library(sp)

load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub
load('../Data/nycSub.RData')
load('../Data/Street_Seg/streets10.dat') # longStrBroke

# Iterate through each buffer width
adjust_val = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 4, 6, 10)

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

    load(paste0("../Data/OutputStrInfo_realData/strInfo_10_", k, ".dat")) # contains the buffer object

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
                                  poly1@Polygons[[1]]@coords[,1], poly1@Polygons[[1]]@coords[,2])
            arr_2 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
                                  poly2@Polygons[[1]]@coords[,1], poly2@Polygons[[1]]@coords[,2])
            
            off_1 = point.in.polygon(off_sub$x_coord_cd, off_sub$y_coord_cd,
                                     poly1@Polygons[[1]]@coords[,1], poly1@Polygons[[1]]@coords[,2])
            off_2 = point.in.polygon(off_sub$x_coord_cd, off_sub$y_coord_cd,
                                     poly2@Polygons[[1]]@coords[,1], poly2@Polygons[[1]]@coords[,2])
            
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
            
            # Calculating the spatial component
            prec_1_x = arr_sub$x_coord_cd[arr_1 > 0]
            prec_1_y = arr_sub$y_coord_cd[arr_1 > 0]
            
            prec_2_x = arr_sub$x_coord_cd[arr_2 > 0]
            prec_2_y = arr_sub$y_coord_cd[arr_2 > 0]
            
            pp_1 = ppp(prec_1_x, prec_1_y, 
                        c(min(c(prec_1_x, prec_2_x)), max(c(prec_1_x, prec_2_x))), 
                        c(min(c(prec_1_y, prec_2_y)), max(c(prec_1_y, prec_2_y))))
            pp_2 = ppp(prec_2_x, prec_2_y, 
                        c(min(c(prec_1_x, prec_2_x)), max(c(prec_1_x, prec_2_x))), 
                        c(min(c(prec_1_y, prec_2_y)), max(c(prec_1_y, prec_2_y))))
            
            border_line_1_2 = longStrBroke[[k]][[i]][[j]]$shorterStreet
            
            b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
            b_line_1_2 = psp(b_c_1_2[1:(nrow(b_c_1_2)-1),"x"], 
                             b_c_1_2[1:(nrow(b_c_1_2)-1),"y"],
                             b_c_1_2[2:nrow(b_c_1_2),"x"],
                             b_c_1_2[2:nrow(b_c_1_2),"y"],
                             Window(pp_1))
            
            int_1 = density.ppp(pp_1, adjust = adjust_val[index], scalekernel = T) 
            int_2 = density.ppp(pp_2, adjust = adjust_val[index], scalekernel = T)
            
            line_intensity_1 = int_1[b_line_1_2]
            line_intensity_2 = int_2[b_line_1_2]
            int_line_1 = mean(line_intensity_1)
            int_line_2 = mean(line_intensity_2)
            intensity_diff = abs(int_line_1 - int_line_2)
            
            plot(int_2)
            plot(streetLengthInfo_null[[i]][[j]]$buffer, add = T, border = 'green')
            plot(longStrBroke[[k]][[i]][[j]]$shorterStreet, add = T, col = 'green')
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
