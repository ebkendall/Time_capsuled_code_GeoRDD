for (k in 1:8) {
    combinedMatchingSetup = list()
    
    load(paste0("../Output_rev/nullGridInfo/nullData", k, "_1.dat"))
    
    neg_index = (nullData$str_info$precinct != -1)
    nullData$str_info = nullData$str_info[neg_index, ]
    nullData$str_surf$INT_SURFACE = nullData$str_surf$INT_SURFACE[neg_index, ]
    nullData$str_surf$OFF_SURFACE = nullData$str_surf$OFF_SURFACE[neg_index, ]
    
    combinedMatchingSetup$DATA <- nullData$str_info
    combinedMatchingSetup$INT_SURFACE <- nullData$str_surf$INT_SURFACE
    combinedMatchingSetup$OFF_SURFACE <- nullData$str_surf$OFF_SURFACE

    for(i in 2:77) {
        print(paste0(k, "_", i))
        load(paste0("../Output_rev/nullGridInfo/nullData", k, "_", i,".dat"))

        neg_index = (nullData$str_info$precinct != -1)
        nullData$str_info = nullData$str_info[neg_index, ]
        nullData$str_surf$INT_SURFACE = nullData$str_surf$INT_SURFACE[neg_index, ]
        nullData$str_surf$OFF_SURFACE = nullData$str_surf$OFF_SURFACE[neg_index, ]

        combinedMatchingSetup$DATA = rbind(combinedMatchingSetup$DATA, nullData$str_info)
        combinedMatchingSetup$INT_SURFACE = rbind(combinedMatchingSetup$INT_SURFACE, 
                                                  nullData$str_surf$INT_SURFACE)
        combinedMatchingSetup$OFF_SURFACE = rbind(combinedMatchingSetup$OFF_SURFACE, 
                                                  nullData$str_surf$OFF_SURFACE)
    }

    # Filter out the streets that do not have any streets because those are not relevant
    street1Ind = (combinedMatchingSetup$DATA$streets1 != 0)
    street2Ind = (combinedMatchingSetup$DATA$streets2 != 0)
    streetInd = (street1Ind & street2Ind)
    
    combinedMatchingSetup$DATA = combinedMatchingSetup$DATA[streetInd, ]
    combinedMatchingSetup$INT_SURFACE = combinedMatchingSetup$INT_SURFACE[streetInd, ]
    combinedMatchingSetup$OFF_SURFACE = combinedMatchingSetup$OFF_SURFACE[streetInd, ]
    
    combinedMatchingSetupFix = combinedMatchingSetup

    ## Create ratios of area and streets
    combinedMatchingSetupFix$DATA$ratioArea = combinedMatchingSetupFix$DATA$area1 /
      combinedMatchingSetupFix$DATA$area2
    combinedMatchingSetupFix$DATA$ratioArea[which(combinedMatchingSetupFix$DATA$ratioArea < 1)] =
      1/combinedMatchingSetupFix$DATA$ratioArea[which(combinedMatchingSetupFix$DATA$ratioArea < 1)]

    combinedMatchingSetupFix$DATA$ratioStreet = combinedMatchingSetupFix$DATA$streets1 /
      combinedMatchingSetupFix$DATA$streets2
    combinedMatchingSetupFix$DATA$ratioStreet[which(combinedMatchingSetupFix$DATA$ratioStreet < 1)] =
      1/combinedMatchingSetupFix$DATA$ratioStreet[which(combinedMatchingSetupFix$DATA$ratioStreet < 1)]

    save(combinedMatchingSetupFix, file = paste0("../Output_rev/nullGridInfo/combinedMatchingSetup", k, ".dat"))
}


