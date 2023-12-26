load('../Data/dataArr.RData')
load('../Data/ind_prec_df.rda')
load('../Data/indexList_MAIN.RData')
load('../Data/nycSub.RData')
load('../Data/totalStreetBuffInfo_NEW.RData')

library(rgeos, quietly = T)
library(ggmap, quietly = T)
library(tidyverse, quietly = T)

# Year choice
dataArrF = dataArrF[dataArrF$year == 2014, ]

# Figure 1 ---------------------------------------------------------------------
png(filename = "Plots/prec77arrLoc1.png", width = 2000, height = 1000,
    units = "px", pointsize = 12, bg = "white", res = NA)

prec_num = 49

plot(nycSub[prec_num, ], lty = 2)
points(dataArrF$x_coord_cd[dataArrF$arrest_precinct == nycSub$Precinct[prec_num]],
       dataArrF$y_coord_cd[dataArrF$arrest_precinct == nycSub$Precinct[prec_num]],
       col = 'blue', cex = 2.5)
dev.off()

png(filename = "Plots/prec77arrLoc2.png", width = 2000, height = 1000,
    units = "px", pointsize = 12, bg = "white", res = NA)
dataArrF_fig1 = dataArrF[dataArrF$arrest_precinct == nycSub$Precinct[prec_num], ]

qmplot(lon, lat, data = dataArrF_fig1, alpha = I(.9), size = I(8), colour = 'darkred', darken = .1, source = "osm") +
  theme(legend.position = "none")

dev.off()

# Figure 2 ---------------------------------------------------------------------

# Borders ind1, ind2
ind1 = 98; ind2 = 104
ind_choice = ind_prec_df[c(ind1, ind2), ]
arr_choice = dataArrF[dataArrF$arrest_precinct %in% c(ind_choice[1,], ind_choice[2, ]), ]

png(filename = "Plots/arrestsAndBuffer.png", width = 2000, height = 1000,
    units = "px", pointsize = 12, bg = "white", res = NA)
par(mfrow = c(1,2))

# ind1 first
nyc_small = nycSub[c(which(nycSub$Precinct == ind_choice[1,1]),
                     which(nycSub$Precinct == ind_choice[1,2])), ]

plot(totalStreetBuffInfo_NEW[[5]][[ind1]]$poly2, lwd = 3)
plot(totalStreetBuffInfo_NEW[[5]][[ind1]]$poly1, add = T, lwd = 3)
plot(nyc_small, lty = 2, add = T)
points(arr_choice$x_coord_cd[arr_choice$arrest_precinct == ind_choice[1,1]],
       arr_choice$y_coord_cd[arr_choice$arrest_precinct == ind_choice[1,1]],
       col = 'red', cex = 2.5)
points(arr_choice$x_coord_cd[arr_choice$arrest_precinct == ind_choice[1,2]],
       arr_choice$y_coord_cd[arr_choice$arrest_precinct == ind_choice[1,2]],
       col = 'blue', cex = 2.5)

# ind2 second
nyc_small = nycSub[c(which(nycSub$Precinct == ind_choice[2,1]),
                     which(nycSub$Precinct == ind_choice[2,2])), ]

plot(totalStreetBuffInfo_NEW[[5]][[ind2]]$poly1, lwd = 3)
plot(totalStreetBuffInfo_NEW[[5]][[ind2]]$poly2, add = T, lwd = 3)
plot(nyc_small, lty = 2, add = T)
points(arr_choice$x_coord_cd[arr_choice$arrest_precinct == ind_choice[2,1]],
       arr_choice$y_coord_cd[arr_choice$arrest_precinct == ind_choice[2,1]],
       col = 'red', cex = 2.5)
points(arr_choice$x_coord_cd[arr_choice$arrest_precinct == ind_choice[2,2]],
       arr_choice$y_coord_cd[arr_choice$arrest_precinct == ind_choice[2,2]],
       col = 'blue', cex = 2.5)

dev.off()

# Figure 6 (streets) -----------------------------------------------------------
load('../Data/streetsByPrec.RData')
load('../Data/OutputStrInfo_realData/strInfo_5_53.dat')

png(filename = "Plots/precStreetsAndBuff.png", width = 2000, height = 1000,
    units = "px", pointsize = 12, bg = "white", res = NA)

areas = NULL
for(i in 1:length(streetLengthInfo_null)) {
  if (!is.na(streetLengthInfo_null[[i]][[1]])) {
    areas = c(areas, streetLengthInfo_null[[i]][[1]]$buffer@polygons[[1]]@area)
  } else {
    areas = c(areas, NA)
  }
}
areas_ind = order(areas, decreasing = T)
plot(nycSub[53, ], lwd = 2.5)
plot(streetLengthInfo_null[[areas_ind[1]]][[1]]$buffer, border = "blue", lwd = 2, add = T)
temp = nycSub[53, ]
plot(streetsByPrec[[53]], add = T, col = "grey", lwd = 1.5)

borderBuff = gBuffer(temp, width = -500)
newSubStreets = gIntersection(streetsByPrec[[53]], borderBuff)
plot(newSubStreets, add = T, col = 'red', lwd = 2)
plot(streetLengthInfo_null[[areas_ind[1]]][[1]]$buffer,  add = T, border = "blue", lwd = 2)
plot(streetLengthInfo_null[[areas_ind[12]]][[1]]$buffer, add = T, border = "blue", lwd = 2)
plot(streetLengthInfo_null[[areas_ind[60]]][[1]]$buffer, add = T, border = "blue", lwd = 2)

dev.off()
