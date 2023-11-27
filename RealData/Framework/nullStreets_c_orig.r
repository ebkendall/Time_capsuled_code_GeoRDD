library(sp); library(sf); library(rgeos); library(raster)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_NEW.RData")
load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub
Dir = '../Output/origGridInfo/'

for (k in 3:10) {
  print(k)
  sim_orig <- list(DATA = data.frame("area1" = rep(NA,164), "area2" = rep(NA,164), 
                                     "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                     "n_arr_1" = rep(NA, 164), "n_arr_2" = rep(NA, 164),
                                     "n_off_1" = rep(NA, 164), "n_off_2" = rep(NA, 164),
                                     "naive_pval" = rep(NA, 164),
                                     "n_arr_1_prec" = rep(NA, 164), "n_arr_2_prec" = rep(NA, 164),
                                     "n_off_1_prec" = rep(NA, 164), "n_off_2_prec" = rep(NA, 164),
                                     "naive_pval_prec" = rep(NA, 164)))
    
  for (i in indexList_MAIN) {
    print(i)
    prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[i])
    prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[i])
    
    poly1 = totalStreetBuffInfo_NEW[[k]][[i]]$poly1
    poly2 = totalStreetBuffInfo_NEW[[k]][[i]]$poly2
    
    poly_ind1 = totalStreetBuffInfo_NEW[[k]][[i]]$poly_ind1
    poly_ind2 = totalStreetBuffInfo_NEW[[k]][[i]]$poly_ind2
    
    area1 = totalStreetBuffInfo_NEW[[k]][[i]]$area1
    area2 = totalStreetBuffInfo_NEW[[k]][[i]]$area2
    
    # Arrest Info (total)
    arr_1 = point.in.polygon(dataArr_sub[,"x_coord_cd"], dataArr_sub[,"y_coord_cd"],
                              poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                              poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    arr_2 = point.in.polygon(dataArr_sub[,"x_coord_cd"], dataArr_sub[,"y_coord_cd"],
                              poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                              poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    
    # Offense Info (total)
    off_1 = point.in.polygon(dataOff_sub[,"x_coord_cd"], dataOff_sub[,"y_coord_cd"],
                             poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                             poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    off_2 = point.in.polygon(dataOff_sub[,"x_coord_cd"], dataOff_sub[,"y_coord_cd"],
                             poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                             poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    
    # Counting the arrests and offenses
    n_arr_1 = sum(arr_1 > 0)
    n_arr_2 = sum(arr_2 > 0)
    n_off_1 = sum(off_1 > 0)
    n_off_2 = sum(off_2 > 0)
    
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
    
    # Arrest and Offense Info (precinct specific)
    dataArr_prec_1 = dataArr_sub[dataArr_sub$arrest_precinct == ind_prec_df$prec1[i], ]
    dataArr_prec_2 = dataArr_sub[dataArr_sub$arrest_precinct == ind_prec_df$prec2[i], ]
    dataOff_prec_1 = dataOff_sub[dataOff_sub$precinct == ind_prec_df$prec1[i], ]
    dataOff_prec_2 = dataOff_sub[dataOff_sub$precinct == ind_prec_df$prec2[i], ]
    
    arr_1_prec_a = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    arr_1_prec_b = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    arr_2_prec_a = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    arr_2_prec_b = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    off_1_prec_a = point.in.polygon(dataOff_prec_1[,"x_coord_cd"], dataOff_prec_1[,"y_coord_cd"],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    off_1_prec_b = point.in.polygon(dataOff_prec_1[,"x_coord_cd"], dataOff_prec_1[,"y_coord_cd"],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    off_2_prec_a = point.in.polygon(dataOff_prec_2[,"x_coord_cd"], dataOff_prec_2[,"y_coord_cd"],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                                    poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    off_2_prec_b = point.in.polygon(dataOff_prec_2[,"x_coord_cd"], dataOff_prec_2[,"y_coord_cd"],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                                    poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    
    # Counting the arrests and offenses (precinct specific)
    n_arr_1_prec = sum(arr_1_prec_a > 0) + sum(arr_1_prec_b > 0)
    n_arr_2_prec = sum(arr_2_prec_a > 0) + sum(arr_2_prec_b > 0)
    n_off_1_prec = sum(off_1_prec_a > 0) + sum(off_1_prec_b > 0)
    n_off_2_prec = sum(off_2_prec_a > 0) + sum(off_2_prec_b > 0)
    
    # Naive p-value (precinct specific)
    count1 = n_arr_1_prec
    count2 = n_arr_2_prec
    n = count1 + count2
    p = 0.5
    pval_prec = NA
    
    if (count1 <= n/2) {
        pval_prec = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
    } else {
        pval_prec = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
    }

    s1 = totalStreetBuffInfo_NEW[[k]][[i]]$streetLength1
    s2 = totalStreetBuffInfo_NEW[[k]][[i]]$streetLength2

    sim_orig$DATA[i,] = c(area1, area2, s1, s2, n_arr_1, n_arr_2, n_off_1, n_off_2, pval,
                          n_arr_1_prec, n_arr_2_prec, n_off_1_prec, n_off_2_prec, pval_prec)
  }

  save(sim_orig, file = paste0(Dir, 'sim_orig_', k, '.dat'))
}
