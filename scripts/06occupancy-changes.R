################################################################################
# Calculating metrics of occupancy change/range shift from 2004 to 2023
# Amanda Zak
# August 2023
################################################################################

# load data
pred <- read.csv("output/PredictedOccupancy2023.csv")
pred04 <- read.csv("output/PredictedOccupancy2004.csv")
pred23_all <- read.csv("output/PredOcc2004model2023covs.csv")
pred23_sf <- read.csv("output/PredOcc2004model2023sf.csv")


# Cell-by-cell change from 2004 to 2023
# only for cells in both 2023 and 2004
pred_overlap <- pred[which(pred$OBJECTID %in% pred04$OBJECTID),]
pred04_overlap <- pred04[which(pred04$OBJECTID %in% pred$OBJECTID),]
pred_overlap$Change <- (pred_overlap$OccReal-pred04_overlap$OccReal)


################################################################################

# Calculate average occupancy prob and estimated occupied cells for each prediction

# 2023
mean(pred$OccReal) # 12.1%
sum(pred$OccReal)*(360*360)*(1e-6) # 2,659 cells

# 2004
mean(pred04$OccReal) # 44.6%
sum(pred04$OccReal) # 82,555 occ cells

# 2023 w/ 2004 model, sf only
mean(pred23_sf$OccReal) # 40.2%
sum(pred23_sf$OccReal) # 68,220 occ cells

# 2023 w/ 2004 model, all
mean(pred23_all$OccReal) # 58.7%
sum(pred23_all$OccReal) # 99,492 occ cells
 

###############################################################################
# Calculating area of >60% occ prob per decade, and proportion of total study area

# 2004
hop04 <- length(which(pred04$OccReal >= 0.6))*(360*360)*(1e-6)
# 4335 sq km
length(which(pred04$OccReal >= 0.6))/nrow(pred04)
# 18.1%

# 2023
hop23 <- length(which(pred$OccReal >= 0.6))*(360*360)*(1e-6)
# 916 sq km
length(which(pred$OccReal >= 0.6))/nrow(pred)
# 4.17%

# change in area
(hop23-hop04)/hop04
# 78.9% loss

# change in proportion area
(0.0417023-0.1805732)/0.1805732
# 76.9% loss

# 2023 prediction climate only
hop23_sf <- length(which(pred23_sf$OccReal >= 0.6))*(360*360)*(1e-6)
# 2971 sq km 
length(which(pred23_sf$OccReal >= 0.6))/nrow(pred)
# 13.5% 
(hop23_sf-hop04)/hop04
# 31.5% loss
(0.1352228-0.1805732)/0.1805732
# 25.1% loss

# 2023 prediction all
hop23_all <- length(which(pred23_all$OccReal >= 0.6))*(360*360)*(1e-6)
# 10,690 sq km
length(which(pred23_all$OccReal >= 0.6))/nrow(pred)
# 48.7%
(hop23_all-hop04)/hop04
# 146.6% gain
(0.4865131-0.1805732)/0.1805732
# 169.4% increase


##############################################################################
# Do the same values for 2004 and 2023 yield different occ probs?

# Values:
snow <- 140
forest <- 0.8
es <- 0.1

# 2023
a <- -2.1594519 # intercept
sfb <- 0.8779362 # snowfall beta
fb <- 0.7056529 # forest in buffer beta
esb <- 1.6499307 # es of forest in buffer beta

# standardize values
sfStnd <- (snow-stnd[1,1])/stnd[2,1]
fbStnd <- (forest-stnd[1,7])/stnd[2,7]
esbStnd <- (es-stnd[1,10])/stnd[2,10]

# calc occ prob and back-transform
Occ <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd
OccReal <- exp(Occ)/(1+exp(Occ))
OccReal

# 2004
a <- -0.0985304 # intercept
sfb <- 0.4419514  # snowfall beta
fb <- -0.0156494 # forest in buffer beta
esb <- 0.7554851 # es of forest in buffer beta

# standardize values
sfStnd <- (snow-stnd04[1,1])/stnd04[2,1]
fbStnd <- (forest-stnd04[1,5])/stnd04[2,5]
esbStnd <- (es-stnd04[1,8])/stnd04[2,8]

# calc occ prob and back-transform
Occ <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd
OccReal <- exp(Occ)/(1+exp(Occ))
OccReal