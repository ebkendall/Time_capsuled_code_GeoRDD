library(spatstat)

load('RealData/Data/nycSub.RData')
load('RealData/Data/dataArr_sub.rda')
load('RealData/Data/totalStreetBuffInfo_NEW.RData')
load('RealData/Data/ind_prec_df.rda')
load('RealData/Data/indexList_MAIN.RData')


###############################################################################
######## Checking the distribution of the "matched" test statistics ###########
###############################################################################
load('RealData/Output/nullGridInfo_randomize/combinedMatchingSetup8.dat')
load('RealData/Output/origGridInfo_randomize/sim_orig_8.dat')

area_ratio = c(na.omit(sim_orig$DATA$area1 / sim_orig$DATA$area2))
area_ratio[area_ratio < 1] = 1 / area_ratio[area_ratio < 1]
wMax_a = max(area_ratio)
wMin_a = min(area_ratio)

street_ratio = c(na.omit(sim_orig$DATA$streets1 / sim_orig$DATA$streets2))
street_ratio[street_ratio < 1] = 1 / street_ratio[street_ratio < 1]
wMax_s = max(street_ratio)
wMin_s = min(street_ratio)

wRatioOk = which(combinedMatchingSetupFix$DATA$ratioArea > wMin_a &
                     combinedMatchingSetupFix$DATA$ratioArea < wMax_a & 
                     combinedMatchingSetupFix$DATA$ratioStreet > wMin_s &
                     combinedMatchingSetupFix$DATA$ratioStreet < wMax_s)

combinedMatchingSetupFix2 = combinedMatchingSetupFix$DATA[wRatioOk,]

# Wherever there is a 0 for the offense count, everything gets scaled by 1
which_zeros = which(combinedMatchingSetupFix2$n_off_1 == 0 | combinedMatchingSetupFix2$n_off_2 == 0)
combinedMatchingSetupFix2$n_arr_1[which_zeros] = combinedMatchingSetupFix2$n_arr_1[which_zeros] + 1 
combinedMatchingSetupFix2$n_arr_2[which_zeros] = combinedMatchingSetupFix2$n_arr_2[which_zeros] + 1 
combinedMatchingSetupFix2$n_off_1[which_zeros] = combinedMatchingSetupFix2$n_off_1[which_zeros] + 1 
combinedMatchingSetupFix2$n_off_2[which_zeros] = combinedMatchingSetupFix2$n_off_2[which_zeros] + 1

which_zeros_orig = which(sim_orig$DATA$n_off_1_prec == 0 | sim_orig$DATA$n_off_2_prec == 0)
sim_orig$DATA$n_arr_1_prec[which_zeros_orig] = sim_orig$DATA$n_arr_1_prec[which_zeros_orig] + 1
sim_orig$DATA$n_arr_2_prec[which_zeros_orig] = sim_orig$DATA$n_arr_2_prec[which_zeros_orig] + 1
sim_orig$DATA$n_off_1_prec[which_zeros_orig] = sim_orig$DATA$n_off_1_prec[which_zeros_orig] + 1
sim_orig$DATA$n_off_2_prec[which_zeros_orig] = sim_orig$DATA$n_off_2_prec[which_zeros_orig] + 1

# Looking into things ---------------------------------------------------------
temp <- sim_orig$DATA[which_zeros_orig,]
prec_zero_pairs <- ind_prec_df[which_zeros_orig, ]
dataArr_prec_1 = dataArr_sub[dataArr_sub$arrest_precinct == prec_zero_pairs[2,1], ]
dataArr_prec_2 = dataArr_sub[dataArr_sub$arrest_precinct == prec_zero_pairs[2,2], ]
plot(nycSub[nycSub$Precinct %in% prec_zero_pairs[2,], ])
points(dataArr_prec_1$x_coord_cd, dataArr_prec_1$y_coord_cd)
points(dataArr_prec_2$x_coord_cd, dataArr_prec_2$y_coord_cd, col = 'red')
#  ----------------------------------------------------------------------------

# v1 = sd(combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2, na.rm=TRUE)^2
v1 = sd(combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2, na.rm=TRUE)^2

# rat_off = combinedMatchingSetupFix2$n_off_1 / combinedMatchingSetupFix2$n_off_2
rat_off = combinedMatchingSetupFix2$streets1 / combinedMatchingSetupFix2$streets2
rat_off[rat_off < 1] = 1 / rat_off[rat_off < 1]
v2 = sd(rat_off, na.rm=TRUE)^2

t_stat = combinedMatchingSetupFix2$spatialDiff
t_stat_orig = sim_orig$DATA$spatialDiff

matched_ind = vector(mode = 'list', length = max(indexList_MAIN))
freq_ind = NULL

for (ii in indexList_MAIN) {
    ## find matches
    # off_temp = sim_orig$DATA$n_off_1_prec[ii] + sim_orig$DATA$n_off_2_prec[ii]
    # ratio_temp = max(sim_orig$DATA$n_off_1_prec[ii] / sim_orig$DATA$n_off_2_prec[ii],
    #                  sim_orig$DATA$n_off_2_prec[ii] / sim_orig$DATA$n_off_1_prec[ii])
    off_temp = sim_orig$DATA$streets1[ii] + sim_orig$DATA$streets2[ii]
    ratio_temp = max(sim_orig$DATA$streets1[ii] / sim_orig$DATA$streets2[ii],
                     sim_orig$DATA$streets2[ii] / sim_orig$DATA$streets1[ii])
    
    stat_temp = t_stat_orig[ii]
    
    # dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$n_off_1 + combinedMatchingSetupFix2$n_off_2))^2/v1) +
    #                      ((ratio_temp - rat_off)^2 / v2))
    dist_temp = sqrt(((off_temp - (combinedMatchingSetupFix2$streets1 + combinedMatchingSetupFix2$streets2))^2/v1) +
                         ((ratio_temp - rat_off)^2 / v2))
    
    w50 = order(dist_temp)[1:300]
    
    matched_ind[[ii]] = data.frame("ind" = w50, "tstat" = t_stat[w50])
    freq_ind = c(freq_ind, w50)
}

pdf('tstat_invest2.pdf')
par(mfrow=c(2,2))
tot_p_val = NULL
for(ii in indexList_MAIN) {
    p_val = mean(matched_ind[[ii]]$tstat > t_stat_orig[ii], na.rm=TRUE)
    hist(matched_ind[[ii]]$tstat, breaks = 100, 
         xlim = c(0, max(matched_ind[[ii]]$tstat, t_stat_orig[ii])),
         main = paste0('Border: ', ind_prec_df[ii, 1], ", ", ind_prec_df[ii,2]),
         xlab = paste0("P-value: ", p_val))
    abline(v = t_stat_orig[ii], col = 'red', lwd = 2)
    tot_p_val = c(tot_p_val, p_val)
}
dev.off()

hist(tot_p_val)


###############################################################################
############## Comparing different cross validation strategies ################
###############################################################################
i = 109
prec_ind = 67
prec_i = nycSub[prec_ind,]

box_x_min = min(c(prec_i@polygons[[1]]@Polygons[[1]]@coords[,1]))
box_x_max = max(c(prec_i@polygons[[1]]@Polygons[[1]]@coords[,1]))
box_y_min = min(c(prec_i@polygons[[1]]@Polygons[[1]]@coords[,2]))
box_y_max = max(c(prec_i@polygons[[1]]@Polygons[[1]]@coords[,2]))

if(length(prec_i@polygons[[1]]@Polygons) > 1) {
    for (pj in 1:length(prec_i@polygons[[1]]@Polygons)) {
        box_x_min = min(c(box_x_min, prec_i@polygons[[1]]@Polygons[[pj]]@coords[,1]))
        box_x_max = max(c(box_x_max, prec_i@polygons[[1]]@Polygons[[pj]]@coords[,1]))
        box_y_min = min(c(box_y_min, prec_i@polygons[[1]]@Polygons[[pj]]@coords[,2]))
        box_y_max = max(c(box_y_max, prec_i@polygons[[1]]@Polygons[[pj]]@coords[,2]))
    }
}

poly_box = matrix(c(box_x_min, box_y_max, 
                    box_x_min, box_y_min, 
                    box_x_max, box_y_min,
                    box_x_max, box_y_max), ncol = 2, byrow = T)

prec_i_x = dataArr_sub$x_coord_cd[dataArr_sub$arrest_precinct == i]
prec_i_y = dataArr_sub$y_coord_cd[dataArr_sub$arrest_precinct == i]
poly_arrest = point.in.polygon(prec_i_x, prec_i_y, poly_box[,1], poly_box[,2])

pp_i = ppp(prec_i_x[poly_arrest > 0], prec_i_y[poly_arrest > 0], 
           c(box_x_min, box_x_max), 
           c(box_y_min, box_y_max))

plot(nycSub[prec_ind,])
points(prec_i_x, prec_i_y)
plot(pp_i, col = 'red', add = T)

b_line = psp(box_x_min, box_y_min, 
             box_x_max, box_y_max,
             Window(pp_i))

# Default spatial window
int_1 = density.ppp(pp_i) 
int_line_1 = mean(int_1[b_line]); print(int_line_1)
# BW Diggle
sigma2 = bw.diggle(pp_i)
int_2 = density.ppp(pp_i, sigma = sigma2)
int_line_2 = mean(int_2[b_line]); print(int_line_2)
# BW CvL
sigma3 = bw.CvL(pp_i)
int_3 = density.ppp(pp_i, sigma = sigma3)
int_line_3 = mean(int_3[b_line]); print(int_line_3)
# BW Scott
sigma4 = bw.scott(pp_i)
int_4 = density.ppp(pp_i, sigma = sigma4)
int_line_4 = mean(int_4[b_line]); print(int_line_4)
# BW ppl
sigma5 = bw.ppl(pp_i)
int_5 = density.ppp(pp_i, sigma = sigma5)
int_line_5 = mean(int_5[b_line]); print(int_line_5)

sim_pp = runifpoint(1000, win = Window(pp_i))
sim_pp_cords = cbind(sim_pp$x, sim_pp$y)
sim_pp_cords = sim_pp_cords[order(sim_pp$x), ]

sub_sim_pp = ppp(sim_pp_cords[1:200, 1], sim_pp_cords[1:200, 2], window = Window(pp_i))
plot(sub_sim_pp, add = T, col = 'red')

int_sim = density.ppp(sim_pp)
plot(int_sim)
plot(sim_pp, add = T)

sim_pp2 = superimpose.ppplist(sim_pp, sub_sim_pp, sub_sim_pp, sub_sim_pp, sub_sim_pp, 
                              sub_sim_pp, sub_sim_pp, sub_sim_pp, sub_sim_pp, sub_sim_pp)
int_sim2 = density.ppp(sim_pp2)
plot(int_sim2)
