################################################################################
# Identifying the >90% forested cells and removing nonforest cells from analysis
# This code saves subsets of the 2023 and 2004 grid cells with >90% forest cover
# Amanda Zak
# February 2024
################################################################################

# Import .csv files of area of forest cover of all grid cells for each year
for23 <- read.csv("data/Forested_Area_2023.csv")
for04 <- read.csv("data/Forested_Area_2004.csv")

# calculate % forested cover of grid cell
for23$percFor <- (for23$FORES_1)/129600
for04$percFor <- (for04$FORES_1)/129600

# extract >=90% forested cells
for23_keep <- subset(for23,for23$percFor >= 0.90)
for04_keep <- subset(for04,for04$percFor >= 0.90)

# save these subsets to refer to later
saveRDS(for23_keep,"output/for23_keep.RData")
saveRDS(for04_keep,"output/for04_keep.RData")

# load grid cell .csv files
grid <- read.csv("data/StudyAreaGrid.csv")
grid04 <- read.csv("data/StudyAreaGrid2004.csv")

# cut grid cells down to kept ones
grid04_keep <- subset(grid04,grid04$OBJECTID %in% for04_keep$OBJECTID_1)
grid_keep <- subset(grid,grid$OBJECTID %in% for23_keep$OBJECTID_1)
saveRDS(grid_keep,"output/grid_keep.RData")
saveRDS(grid04_keep,"output/grid04_keep.RData")
