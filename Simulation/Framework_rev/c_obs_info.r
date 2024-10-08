library(sp); library(sf); library(rgeos)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/gridWithin_prec.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_ORIG.RData")
load("../Data/totalStreetBuffInfo_NEW.RData")

Dir = '../Output_noWater_rev/origGridInfo/'
print(Dir)

for (k in 1:8) {
    
    sim_orig <- list(DATA = data.frame("area1" = rep(NA,164), "area2" = rep(NA,164), 
                                       "streets1" = rep(NA, 164), "streets2" = rep(NA, 164)),
                     GRID_IND_1 = vector(mode = 'list', length = 164),
                     GRID_IND_2 = vector(mode = 'list', length = 164))
    print(paste0("index ", k, " of 8"))
    buff_ind = k + 2

    for (i in indexList_MAIN) {
        
        # New stuff
        prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[i])
        prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[i])
        
        poly1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly1
        poly2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly2
        
        poly_ind1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly_ind1
        poly_ind2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly_ind2
        
        area1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$area1
        area2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$area2
        
        s1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength1
        s2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength2
        
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
        
        p1 = point.in.polygon(gridCoords[,1], gridCoords[,2],
                              poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                              poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
        p2 = point.in.polygon(gridCoords[,1], gridCoords[,2],
                              poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                              poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
        
        ind1 <- which(p1 > 0)
        ind2 <- which(p2 > 0)
        
        gridVals_ind_1 = gridVals_ind_master[ind1]
        gridVals_ind_2 = gridVals_ind_master[ind2]
        
        # Sanity check
        check1 = point.in.polygon(gridCoords[p1,1], gridCoords[p1,2],
                                  nycSub[prec_ind_1, ]@polygons[[1]]@Polygons[[1]]@coords[,1],
                                  nycSub[prec_ind_1, ]@polygons[[1]]@Polygons[[1]]@coords[,2])
        
        if(sum(check1 > 0) != length(check1)) print(paste0("Check i ", i, " k ", k))
        # plot(poly1, border = 'red')
        # plot(poly2, border = 'blue', add = T)
        # plot(nyc_small[1,], border = 'red', add = T)
        # plot(nyc_small[2,], border = 'blue', add = T)
        # points(gridCoords[ind1, 1], gridCoords[ind1, 2], col = 'red')
        # points(gridCoords[ind2, 1], gridCoords[ind2, 2], col = 'blue')
        
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
