library(sp); library(sf); library(rgeos); library(raster)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_NEW.RData")
load("../Data/treesByPrec.RData")
load("../Data/streetsByPrec.RData")
Dir = '../Output_tree/origGridInfo/'
print(Dir)

for (k in 3:10) {
  print(k)
  sim_orig <- list(DATA = data.frame("area1" = rep(NA,164), "area2" = rep(NA,164), 
                                     "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                     "count1" = rep(NA,164), "count2" = rep(NA,164),
                                     "naive_pval" = rep(NA, 164)))
    
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

    # Collecting Tree counts for Precinct 1 across both buffers
    p1 = point.in.polygon(treesByPrec[[prec_ind_1]][,1], treesByPrec[[prec_ind_1]][,2],
                          poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                          poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
    p2 = point.in.polygon(treesByPrec[[prec_ind_2]][,1], treesByPrec[[prec_ind_2]][,2],
                          poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                          poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
    
    count1 = sum(p1 > 0)
    count2 = sum(p2 > 0)
    
    s1 = totalStreetBuffInfo_NEW[[k]][[i]]$streetLength1
    s2 = totalStreetBuffInfo_NEW[[k]][[i]]$streetLength2
    
    n = count1 + count2
    p = 0.5
    pval = NA

    if (count1 <= n/2) {
      pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
    } else {
      pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
    }

    sim_orig$DATA[i,] = c(area1, area2, s1, s2, count1, count2, pval)
    if(count1 == 0 & count2 == 0) print(paste0("No trees: ", k, ", ", i))
  }

  save(sim_orig, file = paste0(Dir, 'sim_orig_', k, '.dat'))
}
