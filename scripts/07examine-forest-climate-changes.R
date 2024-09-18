################################################################################
# This code examines the changes in forest cover and climate from 2004 to 2023 and
# their connection with changes in occupancy
# Amanda Zak
# February 2024
################################################################################

# load predicted occupancy
pred04 <- read.csv("output/PredictedOccupancy2004.csv")
pred <- read.csv("output/PredictedOccupancy2023.csv")

# Cells/area with high probability of occupancy
high_prob_04 <- length(which(pred04$OccReal >= 0.60))
high_prob_23 <- length(which(pred$OccReal >= 0.60))
area_04 <- high_prob_04*129600*0.000001 # convert to sq km
area_23 <- high_prob_23*129600*0.000001 # convert to sq km

# decline in high occ prob area from 2004 to 2023
(area_04-area_23)/area_04

# How much of the high-prob area was lost due to forest loss?
lost <- pred04[which(!pred04$OBJECTID %in% pred$OBJECTID),] # forested cells lost
gained <- pred[which(!pred$OBJECTID %in% pred04$OBJECTID),] # forested cells gained
net_loss <- nrow(lost)-nrow(gained) # 15,714 forested cells lost
area_lost <- net_loss*129600*0.000001
net_loss/nrow(pred04) # or 8.5% of the 2004 study area lost
high_prob_loss <- (length(which(lost$OccReal >= 0.60)) - length(which(gained$OccReal >= 0.60)))
# 5,929 cells lost that were high-occ-prob
high_prob_loss/high_prob_04
# 17.7% of high occ prob areas lost due to removal from study
hist(lost$OccReal)


###############################################################################
# Change in e-s cover when looking at just the kept cells

# import grids
grid <- readRDS("output/grid_keep.RData")
grid04 <- readRDS("output/grid04_keep.RData")
# 2023 sum of e-s cover
es23 <- sum(grid$PercES)*129600*0.000001
# 2004 sum of e-s cover
es04 <- sum(grid04$PercES)*129600*0.000001
# Percentages of total area in each yr
tot23 <- (nrow(grid)*129600*0.000001)
tot04 <- (nrow(grid04)*129600*0.000001)
es23/tot23
es04/tot04
# 0.0113 to 0.0253

###############################################################################
# Do any parts of the study area experience a no-analog climate in 2023?

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

b <- min(grid$sf2) # Set the minimum for the breakpoints
e <- max(grid04$sf1) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 12) # Make a neat vector for the breakpoints

hg1 <- hist(grid04$sf1, breaks = ax, plot = FALSE) # Save first histogram data
hg2 <- hist(grid$sf2, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg1, col = c1, main = "Histograms of Avg Annual Snowfall \n in 2004 (blue) and 2023 (pink)") # Plot 1st histogram using a transparent color
plot(hg2, col = c2, add = TRUE) # Add 2nd histogram using different color

# find max snowfall in 2004
maxSf_04 <- max(grid04$sf1)
# find % of 2023 grid cells above that max
grid[which(grid$sf2 > maxSf_04),]
# 0% - no cold-extreme novel climate

# find min snowfall in 2004
minSf_04 <- min(grid04$sf1)
# find % of 2023 grid cells below that max
nrow(grid[which(grid$sf2 < minSf_04),])/nrow(grid)
# 52 grid cells, or 0.0307% of the study area is warm-extreme novel climate
# so basically no novel climate, just shifting averages
