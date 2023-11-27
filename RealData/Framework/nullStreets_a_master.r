# library(sp); library(sf); library(rgeos); library(raster)
library(spatstat)

load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub
load('../Data/nycSub.RData')
load('../Data/monthKey.rda')

# Iterate through each buffer width
for (index in 3:10) {
  
  for(k in 1:77) {
    
    set.seed(k)

    prec_num = nycSub$Precinct[k]
    arr_sub = dataArr_sub[dataArr_sub$arrest_precinct == prec_num, ]
    off_sub = dataOff_sub[dataOff_sub$precinct == prec_num, ]

    print(paste0("Buffer: ", index, ", Precinct: ", k))

    l = 1000 # arbitrary start length

    nullStr_point_data <- list(DATA = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                      "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                      "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                      "n_arr_1" = rep(-1,l), "n_arr_2" = rep(-1,l), "n_off_1" = rep(-1,l),
                                      "n_off_2" = rep(-1,l), "naive_pval" = rep(-1,l)),
                              ARR_IND_1 = vector(mode = 'list', length = l),
                              ARR_IND_2 = vector(mode = 'list', length = l),
                              OFF_IND_1 = vector(mode = 'list', length = l),
                              OFF_IND_2 = vector(mode = 'list', length = l))
    rowNum = 1

    load(paste0("../Data/OutputStrInfo_realData/strInfo_", index, "_", k, ".dat")) # contains the buffer object

    print(paste0("Total length: ", length(streetLengthInfo_null)))
    for(i in 1:length(streetLengthInfo_null)) {
      print(paste0("i ", i))
      for(j in 1:length(streetLengthInfo_null[[i]])) {
        if(!is.na(streetLengthInfo_null[[i]][[j]])) {
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
            
            nullStr_point_data$DATA[rowNum,] = c(k, i, j,
                                                 streetLengthInfo_null[[i]][[j]]$streetLength1,
                                                 streetLengthInfo_null[[i]][[j]]$streetLength2,
                                                 area1, area2, T, n_arr_1, n_arr_2, n_off_1, n_off_2, pval)
            
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

  }
}
