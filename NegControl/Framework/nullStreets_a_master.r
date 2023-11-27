library(sp); library(sf); library(rgeos); library(raster)

load('../Data/nycSub.RData')
load("../Data/treesByPrec.RData")    # gridWithin_prec treesByPrec

# Iterate through each buffer width
for (index in 3:10) {
  
  for(k in 1:77) {
    
    set.seed(k)

    prec_num = nycSub$Precinct[k]

    print(paste0("Buffer: ", index, ", Precinct: ", k))

    l = 1000 # arbitrary start length

    nullStr_point_data <- list(DATA = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                                 "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                                 "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                                 "count1" = rep(-1,l), "count2" = rep(-1,l), "naive_pval" = rep(-1,l)),
                              GRID_IND_1 = vector(mode = 'list', length = l),
                              GRID_IND_2 = vector(mode = 'list', length = l))
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
            
            nullStr_point_data$DATA[rowNum,] = c(k, i, j, s1, s2, area1, area2, T, count1, count2, pval)
            
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
