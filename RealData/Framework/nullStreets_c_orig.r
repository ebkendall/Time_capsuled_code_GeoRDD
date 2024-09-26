# library(sp); library(sf); library(rgeos); library(raster)
library(spatstat)
library(sp)
library(sf)
# library(rgeos)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_NEW.RData")
load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub
Dir = '../Output/origGridInfo/'

adjust_val = c(0.5, 1, 1.5, 2, 3, 4, 6, 10)
options(warn=1)

# for (k in 1:length(adjust_val)) {
    k = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    print(k)
    set.seed(k)

    buff_ind = k + 2
    
    sim_orig <- list(DATA = data.frame( "area1" = rep(NA,164), "area2" = rep(NA,164), 
                                        "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                        "n_arr_1" = rep(NA, 164), "n_arr_2" = rep(NA, 164),
                                        "n_off_1" = rep(NA, 164), "n_off_2" = rep(NA, 164),
                                        "naive_pval" = rep(NA, 164),
                                        "n_arr_1_prec" = rep(NA, 164), "n_arr_2_prec" = rep(NA, 164),
                                        "n_off_1_prec" = rep(NA, 164), "n_off_2_prec" = rep(NA, 164),
                                        "naive_pval_prec" = rep(NA, 164)),
                    INT_SURFACE = data.frame("spatialDiff1" = rep(NA,164),
                                             "spatialDiff2" = rep(NA,164),
                                             "spatialDiff3" = rep(NA,164),
                                             "spatialDiff4" = rep(NA,164),
                                             "spatialDiff5" = rep(NA,164),
                                             "spatialDiff6" = rep(NA,164),
                                             "spatialDiff7" = rep(NA,164),
                                             "spatialDiff8" = rep(NA,164)))
        
    # pdf("orig_str_global.pdf")
    # par(mfrow = c(2,2))
    # for (i in unique(global_match_loc_orig)) {
    for (i in indexList_MAIN) {
        print(i)
        prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[i])
        prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[i])
        
        poly1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly1
        poly2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly2
        
        poly_ind1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly_ind1
        poly_ind2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly_ind2
        
        area1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$area1
        area2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$area2
        
        # Making the border line object ---------------------------------------
        border_line_1_2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$centerLine
        b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
        
        endpts = which(multiplicity(b_c_1_2) == 1)
        if(length(endpts) > 2) {
            print(paste0(i, " problem!"))
        } 
        if(endpts[1] != 1 | endpts[2] != nrow(b_c_1_2)) {
            print(paste0(i, " out of order"))
            b_c_1_2_new = rbind(b_c_1_2[endpts[2]:nrow(b_c_1_2), ], b_c_1_2[1:endpts[1], ])
            b_c_1_2 = b_c_1_2_new
            border_line_1_2@lines[[1]]@Lines[[1]]@coords = b_c_1_2
        }
        
        test_b_c = unique(b_c_1_2)
        b_line_pp = ppp(test_b_c[,"x"], test_b_c[,"y"],
                        c(min(test_b_c[,"x"]), max(test_b_c[,"x"])),
                        c(min(test_b_c[,"y"]), max(test_b_c[,"y"])))
        b_line_1_2 = psp(test_b_c[1:(nrow(test_b_c)-1),"x"],
                         test_b_c[1:(nrow(test_b_c)-1),"y"],
                         test_b_c[2:nrow(test_b_c),"x"],
                         test_b_c[2:nrow(test_b_c),"y"],
                         Window(b_line_pp))
        # b_line_pp = ppp(b_c_1_2[,"x"], b_c_1_2[,"y"],
        #                 c(min(b_c_1_2[,"x"]), max(b_c_1_2[,"x"])),
        #                 c(min(b_c_1_2[,"y"]), max(b_c_1_2[,"y"])))
        # b_line_1_2 = psp(b_c_1_2[1:(nrow(b_c_1_2)-1),"x"],
        #                     b_c_1_2[1:(nrow(b_c_1_2)-1),"y"],
        #                     b_c_1_2[2:nrow(b_c_1_2),"x"],
        #                     b_c_1_2[2:nrow(b_c_1_2),"y"],
        #                     Window(b_line_pp))
        
        # Create a buffer
        poly3 = st_buffer(st_as_sf(border_line_1_2), dist = buff_ind * 100)
        poly3_coord = sf::st_coordinates(poly3)
        
        # Naive p-value info  -------------------------------------------------
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
        
        # Naive p-value (precinct specific) info  -----------------------------
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

        s1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength1
        s2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$streetLength2

        # Spatial intensity surface construction ------------------------------
        
        # Define window based on the buffer width
        prec1 = poly1
        prec2 = poly2
        
        arr_1_tot = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
                                     poly3_coord[,"X"], poly3_coord[,"Y"])
        arr_2_tot = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
                                     poly3_coord[,"X"], poly3_coord[,"Y"])
        
        prec_1_x = dataArr_prec_1$x_coord_cd[arr_1_tot > 0]
        prec_1_y = dataArr_prec_1$y_coord_cd[arr_1_tot > 0]
        prec_2_x = dataArr_prec_2$x_coord_cd[arr_2_tot > 0]
        prec_2_y = dataArr_prec_2$y_coord_cd[arr_2_tot > 0]
        
        # prec_1_xa = dataArr_prec_1$x_coord_cd[(arr_1_prec_a > 0) | (arr_1_prec_b > 0)]
        # prec_1_ya = dataArr_prec_1$y_coord_cd[(arr_1_prec_a > 0) | (arr_1_prec_b > 0)]
        # prec_2_xa = dataArr_prec_2$x_coord_cd[(arr_2_prec_a > 0) | (arr_2_prec_b > 0)]
        # prec_2_ya = dataArr_prec_2$y_coord_cd[(arr_2_prec_a > 0) | (arr_2_prec_b > 0)]
        
        # prec_1_x = dataArr_sub$x_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec1[i]]
        # prec_1_y = dataArr_sub$y_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec1[i]]
        # prec_2_x = dataArr_sub$x_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec2[i]]
        # prec_2_y = dataArr_sub$y_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec2[i]]

        # # Defining window based on entire precinct
        # prec1 = nycSub[prec_ind_1, ]
        # prec2 = nycSub[prec_ind_2, ]
        
        # prec_1_x = dataArr_sub$x_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec1[i]]
        # prec_1_y = dataArr_sub$y_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec1[i]]
        # prec_2_x = dataArr_sub$x_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec2[i]]
        # prec_2_y = dataArr_sub$y_coord_cd[dataArr_sub$arrest_precinct == ind_prec_df$prec2[i]]

        box_x_min = min(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,1], 
                        prec2@polygons[[1]]@Polygons[[1]]@coords[,1]))
        box_x_max = max(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,1], 
                        prec2@polygons[[1]]@Polygons[[1]]@coords[,1]))
        box_y_min = min(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,2], 
                        prec2@polygons[[1]]@Polygons[[1]]@coords[,2]))
        box_y_max = max(c(prec1@polygons[[1]]@Polygons[[1]]@coords[,2], 
                        prec2@polygons[[1]]@Polygons[[1]]@coords[,2]))
        
        if(length(prec1@polygons[[1]]@Polygons) > 1) {
            for (pj in 1:length(prec1@polygons[[1]]@Polygons)) {
                box_x_min = min(c(box_x_min, prec1@polygons[[1]]@Polygons[[pj]]@coords[,1]))
                box_x_max = max(c(box_x_max, prec1@polygons[[1]]@Polygons[[pj]]@coords[,1]))
                box_y_min = min(c(box_y_min, prec1@polygons[[1]]@Polygons[[pj]]@coords[,2]))
                box_y_max = max(c(box_y_max, prec1@polygons[[1]]@Polygons[[pj]]@coords[,2]))
            }
        }
        
        if(length(prec2@polygons[[1]]@Polygons) > 1) {
            for (pj in 1:length(prec2@polygons[[1]]@Polygons)) {
                box_x_min = min(c(box_x_min, prec2@polygons[[1]]@Polygons[[pj]]@coords[,1]))
                box_x_max = max(c(box_x_max, prec2@polygons[[1]]@Polygons[[pj]]@coords[,1]))
                box_y_min = min(c(box_y_min, prec2@polygons[[1]]@Polygons[[pj]]@coords[,2]))
                box_y_max = max(c(box_y_max, prec2@polygons[[1]]@Polygons[[pj]]@coords[,2]))
            }
        }

        # Removing center points --------------------------------------------------
        # poly3 = gBuffer(border_line_1_2, width= buff_ind * 100)
        # p3_1 = point.in.polygon(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"],
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,2])
        # p3_2 = point.in.polygon(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"],
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,1],
        #                         poly3@polygons[[1]]@Polygons[[1]]@coords[,2])

        # prec_1_x = dataArr_prec_1[(p3_1 == 0) | ((p3_1 > 0 & arr_1_prec_a > 0) | (p3_1 > 0 & arr_1_prec_b > 0)),"x_coord_cd"]
        # prec_1_y = dataArr_prec_1[(p3_1 == 0) | ((p3_1 > 0 & arr_1_prec_a > 0) | (p3_1 > 0 & arr_1_prec_b > 0)),"y_coord_cd"]
        # prec_2_x = dataArr_prec_2[(p3_2 == 0) | ((p3_2 > 0 & arr_2_prec_a > 0) | (p3_2 > 0 & arr_2_prec_b > 0)),"x_coord_cd"]
        # prec_2_y = dataArr_prec_2[(p3_2 == 0) | ((p3_2 > 0 & arr_2_prec_a > 0) | (p3_2 > 0 & arr_2_prec_b > 0)),"y_coord_cd"]
        #  ------------------------------------------------------------------------

        # Focus on points only in the box
        poly_box = matrix(c(box_x_min, box_y_max, 
                            box_x_min, box_y_min, 
                            box_x_max, box_y_min,
                            box_x_max, box_y_max), ncol = 2, byrow = T)
        poly_arr1 = point.in.polygon(prec_1_x, prec_1_y, poly_box[,1], poly_box[,2])
        poly_arr2 = point.in.polygon(prec_2_x, prec_2_y, poly_box[,1], poly_box[,2])

        prec_1_x = prec_1_x[poly_arr1 > 0]
        prec_1_y = prec_1_y[poly_arr1 > 0]
        prec_2_x = prec_2_x[poly_arr2 > 0]
        prec_2_y = prec_2_y[poly_arr2 > 0]
        
        pp_1 = ppp(prec_1_x, prec_1_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))
        pp_2 = ppp(prec_2_x, prec_2_y, c(box_x_min, box_x_max), c(box_y_min, box_y_max))

        # Assigning weights
        single_1 <- !duplicated(pp_1)
        m1 <- multiplicity(pp_1)
        pp_1_weight <- pp_1[single_1] %mark% m1[single_1]

        single_2 <- !duplicated(pp_2)
        m2 <- multiplicity(pp_2)
        pp_2_weight <- pp_2[single_2] %mark% m2[single_2]

        # if(length(pp_1_weight$x) == 1 | length(pp_2_weight$x) == 1) next
        
        # plot(poly_box[,1], poly_box[,2])
        # plot(prec1, add = T)
        # plot(prec2, add = T)
        # plot(poly1, add = T)
        # plot(poly2, add = T)
        # points(dataArr_prec_1[,"x_coord_cd"], dataArr_prec_1[,"y_coord_cd"])
        # points(dataArr_prec_2[,"x_coord_cd"], dataArr_prec_2[,"y_coord_cd"])
        # points(prec_1_x, prec_1_y, col = 'red')
        # points(prec_2_x, prec_2_y, col = 'green')
        # points(pp_1_weight$x, pp_1_weight$y, col = 'red')
        # points(pp_2_weight$x, pp_2_weight$y, col = 'green')
        # # points(prec_1_xa, prec_1_ya, col = 'black')
        # # points(prec_2_xa, prec_2_ya, col = 'blue')
        # # return(0)

        int_surf_vals = NULL
        for(av in 1:length(adjust_val)) {
            if(pp_1_weight$n > 0) {
                if(pp_1_weight$n == 1) {
                    # Cannot do cross validation with one location
                    print(paste0("One data point (side 1, scale ", adjust_val[av], "): ", i))
                    int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks,
                                        adjust = adjust_val[av], scalekernel = T)
                } else {
                    int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks,
                                        sigma = bw.diggle, adjust = adjust_val[av],
                                        scalekernel = T)   
                }
                line_intensity_1 = int_1[b_line_1_2]
                int_line_1 = mean(line_intensity_1) 
            } else {
                int_line_1 = 0
            }

            if(pp_2_weight$n > 0) {
                if(pp_2_weight$n == 1) {
                    # Cannot do cross validation with one location
                    if(pp_1_weight$n > 1) {
                        # Take the spatial smoothness of the other side
                        sig_2 = attr(int_1, "sigma") # This has already been scaled
                        int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks, 
                                            sigma = sig_2, scalekernel = T)   
                    } else {
                        print(paste0("One data point (side 1, scale ", adjust_val[av], "): ", i))
                        int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks, 
                                            adjust = adjust_val[av], scalekernel = T)
                    }
                    
                } else {
                    int_2 = density.ppp(pp_2_weight, weights = pp_2_weight$marks,
                                        sigma = bw.diggle, adjust = adjust_val[av],
                                        scalekernel = T)
                }
                line_intensity_2 = int_2[b_line_1_2]
                int_line_2 = mean(line_intensity_2)
            } else {
                int_line_2 = 0
            }
            
            # Double check, if pp_1_weight$n == 1 and pp_2_weight$n > 1, use sigma from surface 2
            if((pp_1_weight$n == 1) & (pp_2_weight$n > 1)) {
                print(paste0("Informed surface 1 from surface 2"))
                # If the other surface has information, use the same scaled-sigma
                sig_1 = attr(int_2, "sigma") # This has already been scaled
                int_1 = density.ppp(pp_1_weight, weights = pp_1_weight$marks, 
                                    sigma = sig_1, scalekernel = T)
                line_intensity_1 = int_1[b_line_1_2]
                int_line_1 = mean(line_intensity_1)
            }

            intensity_diff = abs(int_line_1 - int_line_2)
            
            int_surf_vals = c(int_surf_vals, intensity_diff)
            
            # plot(int_1, main = paste("Side = 1, ii = ", i, '\n',
            #                          "adjust = ", adjust_val[av], ", sigma = ", attr(int_1, "sigma")))
            # plot(prec1, lwd = 2, add = T)
            # plot(prec2, lwd = 2, add = T)
            # plot(int_2, main = paste("Side = 2, Prec = ", p, ", i = ", i, ", j = ", j, '\n',
            #                          "adjust = ", adjust_val[av], ", sigma = ", attr(int_2, "sigma")))
            # plot(prec1, lwd = 2, add = T)
            # plot(prec2, lwd = 2, add = T)
        }
        
        sim_orig$INT_SURFACE[i,] = int_surf_vals

        sim_orig$DATA[i,] = c(area1, area2, s1, s2, n_arr_1, n_arr_2, n_off_1, n_off_2, pval,
                            n_arr_1_prec, n_arr_2_prec, n_off_1_prec, n_off_2_prec, pval_prec)
        
       # print(warnings())
    }
    
    save(sim_orig, file = paste0(Dir, 'sim_orig_', k, '.dat'))
# }

    
    # pdf("checking_borderline.pdf")
    # par(mfrow=c(2,2))
    # for(i in indexList_MAIN) {
    #     poly1 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly1
    #     poly2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$poly2
    #     
    #     border_line_1_2 = totalStreetBuffInfo_NEW[[buff_ind]][[i]]$centerLine
    #     b_c_1_2 = border_line_1_2@lines[[1]]@Lines[[1]]@coords
    #     
    #     # Make sure the line is actually linear
    #     endpts = which(multiplicity(b_c_1_2) == 1)
    #     if(length(endpts) > 2) {
    #         print(paste0(i, " problem!"))
    #     } 
    #     if(endpts[1] != 1 | endpts[2] != nrow(b_c_1_2)) {
    #         print(paste0(i, " out of order"))
    #         b_c_1_2_new = rbind(b_c_1_2[endpts[2]:nrow(b_c_1_2), ], b_c_1_2[1:endpts[1], ])
    #         b_c_1_2 = b_c_1_2_new
    #     }
    #     
    #     test_b_c = unique(b_c_1_2)
    #     b_line_pp = ppp(test_b_c[,"x"], test_b_c[,"y"],
    #                     c(min(test_b_c[,"x"]), max(test_b_c[,"x"])),
    #                     c(min(test_b_c[,"y"]), max(test_b_c[,"y"])))
    #     b_line_1_2 = psp(test_b_c[1:(nrow(test_b_c)-1),"x"],
    #                      test_b_c[1:(nrow(test_b_c)-1),"y"],
    #                      test_b_c[2:nrow(test_b_c),"x"],
    #                      test_b_c[2:nrow(test_b_c),"y"],
    #                      Window(b_line_pp))
    #     plot(border_line_1_2, lwd = 2)
    #     plot(poly1, add = T)
    #     plot(poly2, add = T)
    #     plot(b_line_1_2, add = T, col = 'red')
    # }
    # dev.off()