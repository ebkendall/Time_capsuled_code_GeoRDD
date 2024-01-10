# Iterate through each buffer width

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)

for (index in 1:8) {

    buff_ind = index + 2

    for (k in 1:77) {

        print(paste0("Buffer: ", buff_ind, ", Precinct: ", k))

        l = 1000 # arbitrary start length

        nullStr_point_data <- list(DATA = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                                     "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                                     "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l)),
                                    GRID_IND_1 = vector(mode = 'list', length = l),
                                    GRID_IND_2 = vector(mode = 'list', length = l))
        rowNum = 1

        load(paste0("../Data/OutputStrInfo_noWater/strInfo_", buff_ind, "_", k, ".dat")) # contains the buffer object
        
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

                        nullStr_point_data$DATA[rowNum,] = c(k, i, j,
                                                    streetLengthInfo_null[[i]][[j]]$streetLength1,
                                                    streetLengthInfo_null[[i]][[j]]$streetLength2,
                                                    area1, area2, T)
                        
                        nullStr_point_data$GRID_IND_1[[rowNum]] = streetLengthInfo_null[[i]][[j]]$pointIndex[[1]]
                        nullStr_point_data$GRID_IND_2[[rowNum]] = streetLengthInfo_null[[i]][[j]]$pointIndex[[2]]
                        
                        rowNum = rowNum + 1
                    }
                }
            }
        }

        save(nullStr_point_data, file=paste0("../Output_noWater/nullGridInfo/nullData", 
                index, "_", k,".dat", sep=''))
    }
    #remember to remove the -1

}



