################################################################################
# 2004 data capture history script using newest land covariates
# Amanda Zak
# September 2023
################################################################################
library("stringr")
library("dplyr")

setwd("C:/Users/alz5215/OneDrive - The Pennsylvania State University/Documents/Research/2004data")

### First need to join coordinate and transect tables

# import transect coordinates
coords <- read.csv("Snowshoe hare database with elevation and all 3 species.csv", skipNul = TRUE)

# import old inp
inp <- read.csv("2005haredata.inp.csv", skipNul = TRUE)

# need to match up site and site ID somehow
coords2 <- na.omit(distinct(coords[,c("SiteID","SiteID1","Latitude1","Longitude1")]))
coords2$SiteID2 <- paste0(coords2$SiteID,ifelse(coords2$SiteID1 == "t","t","r"))
inp$SiteID2 <- str_sub_all(inp$Site,3)
inp$SiteID2 <- str_sub_all(inp$SiteID2,end=-4)
inp$SiteID2 <- as.character(inp$SiteID2)
coordsJoin <- inner_join(coords2,inp,by="SiteID2",relationship="many-to-one")
# 2 duplicated rows - 81r (83,84) and 147r
# get average coordinates for both
coordsJoin[83,3] <- (coordsJoin[83,3] + coordsJoin[84,3])/2
coordsJoin[83,4] <- (coordsJoin[83,4] + coordsJoin[84,4])/2
coordsJoin <- coordsJoin[-84,]
coordsJoin[138,3] <- (coordsJoin[138,3] + coordsJoin[139,3])/2
coordsJoin[138,4] <- (coordsJoin[138,4] + coordsJoin[139,4])/2
coordsJoin <- coordsJoin[-139,]

# save
coordsOut <- coordsJoin[,c(1,5,3,4)]
saveRDS(coordsJoin,"C:/Users/alz5215/OneDrive - The Pennsylvania State University/Documents/Research/2004data/coordsJoin.RData")
write.csv(coordsOut,"C:/Users/alz5215/OneDrive - The Pennsylvania State University/Documents/Research/2004data/coordsOut.csv", row.names = FALSE)

### Import land cover variables and add snowfall

# Import 2004 buffer covariates
trans <- read.csv("Transects2004_Covariates.csv")

# calculate ESfB
trans$PercESpF4500 <- trans$PercES4500/trans$PercF4500

# join to capture history and pare down to essential columns
captHist <- inner_join(coordsJoin,trans,by="SiteID2")
captHist <- captHist[,c(6:9,33:40)]


### Standardize variables and save values
for (i in 5:12) {
  var <- colnames(captHist)[i]
  eval(parse(text=paste0("captHist$",var,"Stnd <- (captHist$",var," - mean(captHist$",
                         var,"))/sd(captHist$",var,")")))
}

stnd04 <- matrix(data = NA, ncol = 8, nrow = 2)
colnames(stnd04) <- colnames(captHist)[5:12]
for (i in 1:ncol(stnd04)) {
  eval(parse(text=paste0("stnd04[1,",i,"] <- mean(captHist$",
                         colnames(captHist)[i+4],")")))
  eval(parse(text=paste0("stnd04[2,",i,"] <- sd(captHist$",colnames(captHist)[i+4],")")))
}

write.csv(stnd04,"Newstandardization2004.csv",row.names = FALSE)


### Write inp file
fn <- "New2004capthist.inp"
captHist$capthist <- paste0(captHist$Site," ",captHist$Obs1,captHist$Obs2,captHist$Obs3,
                            " ","1"," ",captHist$PercESStnd," ",captHist$PercES4500Stnd," ",
                            captHist$PercESpF4500Stnd," ",captHist$PercF4500Stnd," ",
                            captHist$sf1Stnd,";")

write(paste("/* 2004 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Proportion of early-successional cover within the survey unit",
            "\n   Proportion early-successional cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Average snowfall during the period 1994-2004, in cm",
            "*/"),
      file=fn)
write(captHist$capthist, fn, append=T)




