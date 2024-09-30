library(spatstat)
library(sp)
library(sf)

load("../Data/nycSub.RData")
load("../Data/ind_prec_df.rda")
load("../Data/indexList_MAIN.RData")
load("../Data/totalStreetBuffInfo_NEW.RData")
load('../Data/dataArr_sub.rda') # dataArr_sub
load('../Data/dataOff_sub.rda') # dataOff_sub


intensity_surf_same = function(pp_1, pp_2, pp_both, b_line_1_2, adjust_val) {

    sigma_both = flag = 0

    if(pp_both$n > 0) {
        # Using the combined points, compute the smoothing parameter
        if(pp_both$n == 1) {
            # Cannot do cross validation with one location
            int_both = density.ppp(pp_both, weights = pp_both$marks,
                                   adjust = 1, scalekernel = T)
            flag = 1
        } else {
            # Try using cross validation
            check_warn = tryCatch(density.ppp(pp_both, weights = pp_both$marks,
                                              sigma = bw.diggle, adjust = 1,
                                              scalekernel = T),
                                  error = function(e) e, 
                                  warning = function(w) w)
            
            # We receive a warning about the cross-validated smoothing value
            if(inherits(check_warn,"warning")) {
                int_both_temp = density.ppp(pp_both, weights = pp_both$marks,
                                            sigma = bw.diggle, adjust = 1,
                                            scalekernel = T)
                hmax_temp = 1.5 * attr(int_both_temp, "sigma")
                repeat{
                    check_warn = tryCatch(density.ppp(pp_both, weights = pp_both$marks,
                                                      sigma = bw.diggle(pp_both, hmax = hmax_temp), 
                                                      adjust = 1, scalekernel = T),
                                          error = function(e) e, 
                                          warning = function(w) w)
                    
                    # For warnings - continue the next iteration with larger hmax
                    if (inherits(check_warn,"warning")) { 
                        hmax_temp = 1.5 * hmax_temp
                        next  
                    }
                    
                    # For errors - cross validation not working
                    if (inherits(check_warn,"error")) {
                        int_both = density.ppp(pp_both, weights = pp_both$marks,
                                               adjust = 1, scalekernel = T)
                        flag = 2
                        break
                    } 
                    
                    # For it working, run again & leave loop
                    int_both = density.ppp(pp_both, weights = pp_both$marks,
                                           sigma = bw.diggle(pp_both, hmax = hmax_temp),
                                           adjust = 1, scalekernel = T)
                    break
                }
            } else {
                # Cross validation works perfectly
                int_both = density.ppp(pp_both, weights = pp_both$marks,
                                       sigma = bw.diggle, adjust = 1,
                                       scalekernel = T)
            }
        }
        
        sigma_both = attr(int_both, "sigma")
    } else {
        sigma_both = -1
        flag = -1
    }

    # Now fit the intensity surfaces using the scaled spatial smoothness
    int_surf_vals = NULL
    for(av in 1:length(adjust_val)) {
        int_line_1 = int_line_2 = 0

        if(pp_1$n > 0) {
            int_1 = density.ppp(pp_1, weights = pp_1$marks,
                                    sigma = sigma_both, adjust = adjust_val[av],
                                    scalekernel = T)

            line_intensity_1 = int_1[b_line_1_2]
            int_line_1 = mean(line_intensity_1) 
        } 

        if(pp_2$n > 0) {
            int_2 = density.ppp(pp_2, weights = pp_2$marks,
                                    sigma = sigma_both, adjust = adjust_val[av],
                                    scalekernel = T)

            line_intensity_2 = int_2[b_line_1_2]
            int_line_2 = mean(line_intensity_2) 
        }

        intensity_diff = abs(int_line_1 - int_line_2)
        int_surf_vals = c(int_surf_vals, int_line_1, int_line_2, intensity_diff)
    }

    intensity_surface_info = list(sig = sigma_both,
                                  flag = flag,
                                  int_surf_vals = int_surf_vals)

    return(intensity_surface_info)
}

intensity_surf_diff = function(pp_1, pp_2, b_line_1_2, adjust_val) {
    
    sigma_1 = sigma_2 = flag_1 = flag_2 = 0
    
    # Using the combined points, compute the smoothing parameter
    if(pp_1$n > 0) {
        if(pp_1$n == 1) {
            # Cannot do cross validation with one location
            int_1 = density.ppp(pp_1, weights = pp_1$marks,
                                adjust = 1, scalekernel = T)
            flag_1 = 1
        } else {
            # Try using cross validation
            check_warn = tryCatch(density.ppp(pp_1, weights = pp_1$marks,
                                              sigma = bw.diggle, adjust = 1,
                                              scalekernel = T),
                                  error = function(e) e, 
                                  warning = function(w) w)
            
            # We receive a warning about the cross-validated smoothing value
            if(inherits(check_warn,"warning")) {
                int_1_temp = density.ppp(pp_1, weights = pp_1$marks,
                                         sigma = bw.diggle, adjust = 1,
                                         scalekernel = T)
                hmax_temp = 1.5 * attr(int_1_temp, "sigma")
                repeat{
                    check_warn = tryCatch(density.ppp(pp_1, weights = pp_1$marks,
                                                      sigma = bw.diggle(pp_1, hmax = hmax_temp), 
                                                      adjust = 1, scalekernel = T),
                                          error = function(e) e, 
                                          warning = function(w) w)
                    
                    # For warnings - continue the next iteration with larger hmax
                    if (inherits(check_warn,"warning")) { 
                        hmax_temp = 1.5 * hmax_temp
                        next  
                    }
                    
                    # For errors - cross validation not working
                    if (inherits(check_warn,"error")) {
                        int_1 = density.ppp(pp_1, weights = pp_1$marks,
                                            adjust = 1, scalekernel = T)
                        flag_1 = 2
                        break
                    } 
                    
                    # For it working, run again & leave loop
                    int_1 = density.ppp(pp_1, weights = pp_1$marks,
                                        sigma = bw.diggle(pp_1, hmax = hmax_temp),
                                        adjust = 1, scalekernel = T)
                    break
                }
            } else {
                # Cross validation works perfectly
                int_1 = density.ppp(pp_1, weights = pp_1$marks,
                                    sigma = bw.diggle, adjust = 1,
                                    scalekernel = T)
            }
        }
        
        sigma_1 = attr(int_1, "sigma")
    } else {
        sigma_1 = -1
        flag_1 = -1
    }
    
    if(pp_2$n > 0) {
        if(pp_2$n == 1) {
            # Cannot do cross validation with one location
            int_2 = density.ppp(pp_2, weights = pp_2$marks,
                                adjust = 1, scalekernel = T)
            flag_2 = 1
        } else {
            # Try using cross validation
            check_warn = tryCatch(density.ppp(pp_2, weights = pp_2$marks,
                                              sigma = bw.diggle, adjust = 1,
                                              scalekernel = T),
                                  error = function(e) e, 
                                  warning = function(w) w)
            
            # We receive a warning about the cross-validated smoothing value
            if(inherits(check_warn,"warning")) {
                int_2_temp = density.ppp(pp_2, weights = pp_2$marks,
                                         sigma = bw.diggle, adjust = 1,
                                         scalekernel = T)
                hmax_temp = 1.5 * attr(int_2_temp, "sigma")
                repeat{
                    check_warn = tryCatch(density.ppp(pp_2, weights = pp_2$marks,
                                                      sigma = bw.diggle(pp_2, hmax = hmax_temp), 
                                                      adjust = 1, scalekernel = T),
                                          error = function(e) e, 
                                          warning = function(w) w)
                    
                    # For warnings - continue the next iteration with larger hmax
                    if (inherits(check_warn,"warning")) { 
                        hmax_temp = 1.5 * hmax_temp
                        next  
                    }
                    
                    # For errors - cross validation not working
                    if (inherits(check_warn,"error")) {
                        int_2 = density.ppp(pp_2, weights = pp_2$marks,
                                            adjust = 1, scalekernel = T)
                        flag_2 = 2
                        break
                    } 
                    
                    # For it working, run again & leave loop
                    int_2 = density.ppp(pp_2, weights = pp_2$marks,
                                        sigma = bw.diggle(pp_2, hmax = hmax_temp),
                                        adjust = 1, scalekernel = T)
                    break
                }
            } else {
                # Cross validation works perfectly
                int_2 = density.ppp(pp_2, weights = pp_2$marks,
                                    sigma = bw.diggle, adjust = 1,
                                    scalekernel = T)
            }
        }
        
        sigma_2 = attr(int_2, "sigma")
    } else {
        sigma_2 = -1
        flag_2 = -1
    }
    
    # Now fit the intensity surfaces using the scaled spatial smoothness
    int_surf_vals = NULL
    for(av in 1:length(adjust_val)) {
        int_line_1 = int_line_2 = 0

        if(pp_1$n > 0) {
            int_1 = density.ppp(pp_1, weights = pp_1$marks,
                                sigma = sigma_1, adjust = adjust_val[av],
                                scalekernel = T)

            line_intensity_1 = int_1[b_line_1_2]
            int_line_1 = mean(line_intensity_1)
        }

        if(pp_2$n > 0) {
            int_2 = density.ppp(pp_2, weights = pp_2$marks,
                                sigma = sigma_2, adjust = adjust_val[av],
                                scalekernel = T)

            line_intensity_2 = int_2[b_line_1_2]
            int_line_2 = mean(line_intensity_2)
        }

        intensity_diff = abs(int_line_1 - int_line_2)
        int_surf_vals = c(int_surf_vals, int_line_1, int_line_2, intensity_diff)
    }
    
    intensity_surface_info = list(sig = c(sigma_1, sigma_2),
                                  flag = c(flag_1, flag_2),
                                  int_surf_vals = int_surf_vals)
    
    return(intensity_surface_info)
}

