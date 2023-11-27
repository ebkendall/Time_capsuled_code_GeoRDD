library(sp); library(sf); library(rgeos); library(raster)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/gridWithin_prec.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_ORIG.RData")

Dir = '../Output_noWater/origGridInfo/'
print(Dir)

for (k in 3:10) {

  sim_orig <- list(DATA = data.frame("area1" = rep(NA,164), "area2" = rep(NA,164), 
                                     "streets1" = rep(NA, 164), "streets2" = rep(NA, 164)),
                   GRID_IND_1 = vector(mode = 'list', length = 164),
                   GRID_IND_2 = vector(mode = 'list', length = 164))
  print(paste0("index ", k, " of 10"))
    
  for (i in indexList_MAIN) {

    # New stuff
    prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[i])
    prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[i])

    nyc_small = nycSub[c(prec_ind_1, prec_ind_2), ]

    t_grid1 = data.frame(gridWithin_prec[[prec_ind_1]]@coords,
                         "prec" = rep(1, nrow(gridWithin_prec[[prec_ind_1]]@coords)),
                         "o_ind" = 1:nrow(gridWithin_prec[[prec_ind_1]]@coords))
    t_grid2 = data.frame(gridWithin_prec[[prec_ind_2]]@coords,
                         "prec" = rep(2, nrow(gridWithin_prec[[prec_ind_2]]@coords)),
                         "o_ind" = 1:nrow(gridWithin_prec[[prec_ind_2]]@coords))

    gridCoords = rbind(t_grid1, t_grid2)
    colnames(gridCoords) = c("x", "y", "prec", "o_ind")

    gridVals_ind_master = c(gridWithin_prec[[prec_ind_1]]$index,
                            gridWithin_prec[[prec_ind_2]]$index)

    tempOverlap = gIntersection(totalStreetBuffInfo_ORIG[[k]][[i]]$buffer, nyc_small, byid = T)
    poly_order = tempOverlap@plotOrder[1:2]

    int_1 = gIntersection(tempOverlap[poly_order[1], ], nycSub[prec_ind_1, ])
    int_2 = gIntersection(tempOverlap[poly_order[1], ], nycSub[prec_ind_2, ])

    if(is.null(int_1)) {
      int_1 = 0
    } else {
      int_1 = int_1@polygons[[1]]@area
    }

    if(is.null(int_2)) {
      int_2 = 0
    } else {
      int_2 = int_2@polygons[[1]]@area
    }

    if(int_1 < int_2) {
      temp = poly_order[1]
      poly_order[1] = poly_order[2]
      poly_order[2] = temp
    }

    poly1 = tempOverlap[poly_order[1], ]
    poly2 = tempOverlap[poly_order[2], ]

    poly_ind1 = poly1@polygons[[1]]@plotOrder[1]
    poly_ind2 = poly2@polygons[[1]]@plotOrder[1]

    area1 = poly1@polygons[[1]]@Polygons[[poly_ind1]]@area
    area2 = poly2@polygons[[1]]@Polygons[[poly_ind2]]@area

    p1 = point.in.polygon(gridCoords[,1], gridCoords[,2],
                          poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                          poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    p2 = point.in.polygon(gridCoords[,1], gridCoords[,2],
                          poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                          poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])

    ind1 <- which(p1 > 0)
    ind2 <- which(p2 > 0)

    s1 = totalStreetBuffInfo_ORIG[[k]][[i]]$streetLength1
    s2 = totalStreetBuffInfo_ORIG[[k]][[i]]$streetLength2

    gridVals_ind_1 = gridVals_ind_master[ind1]
    gridVals_ind_2 = gridVals_ind_master[ind2]


    if(sum(is.na(gridVals_ind_1)) != 0 | sum(is.na(gridVals_ind_2)) != 0) {
      print("NA was found **************************")
    }
    if(sum(gridVals_ind_1 %in% gridVals_ind_2) != 0 | sum(gridVals_ind_2 %in% gridVals_ind_1) != 0) {
      print("Unsucessful division ******************")
    }

    sim_orig$DATA[i,] = c(area1, area2, s1, s2)
    sim_orig$GRID_IND_1[[i]] = gridVals_ind_1
    sim_orig$GRID_IND_2[[i]] = gridVals_ind_2
  }

  save(sim_orig, file = paste0(Dir, 'sim_orig_', k, '.dat'))
}
