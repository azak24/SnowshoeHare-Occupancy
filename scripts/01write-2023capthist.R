###############################################################################
# New capt hist script that uses the most recent 2004 and 2023 land covariates
# Amanda Zak
# September 2023
###############################################################################
library(dplyr)
library(tidyr)

# upload DNA results
dna <- read.csv("C:/Users/alz5215/OneDrive - The Pennsylvania State University/Documents/Research/Year 1 Data/DNA Results/DNAforR.csv")

# remove the data flagged "remove"
dna <- subset(dna, dna$Note != "Remove")

# list of completed surveys
visits <- read.csv("VisitsForR.csv")

# add detection covariates to visits
snowcover <- read.csv("snowcover.csv")
snowcover$ID <- paste0(snowcover$transectid,"-",snowcover$survey_visit_no)
snowcover$scYN <- ifelse(snowcover$snow_condition == "0 - 25%",0,1)
snowcover$scPerc1 <- ifelse(snowcover$snow_condition == "25 - 75%",1,0)
snowcover$scPerc2 <- ifelse(snowcover$snow_condition == "75 - 100%",1,0)
visits <- left_join(visits,snowcover,by="ID",keep=FALSE)

# Create snowshoe hare capt hist
sh <- matrix(data = NA, nrow = length(unique(visits$Transect)), ncol = 10)
colnames(sh) <- c("Transect","V1","V2","V3","sc1","sc2","sc3","scPerc1","scPerc2","scPerc3")
sh[,"Transect"] <- c(unique(visits$Transect))
sh <- as.data.frame(sh)
for (i in 1:nrow(sh)) { # for each transect
  for (j in 1:3) { # for each visit
    if (j %in% visits[which(visits$Transect == sh$Transect[i]),"Visit"]) { # if the transect was visited
      # if sh detected, matrix gets a 1
      # if not detected, matrix gets a 0
      sh[i,j+1] <- ifelse("Snowshoe Hare" %in% dna[which(dna$Transect == sh$Transect[i] & dna$Visit == j),"Species"],
                          "1","0")
      # fill in snow cover covs
      sh[i,j+4] <- visits[which(visits$Transect == sh$Transect[i] & visits$Visit == j),9]
      sh[i,j+7] <- paste0(visits[which(visits$Transect == sh$Transect[i] & visits$Visit == j),10]," ",
                          visits[which(visits$Transect == sh$Transect[i] & visits$Visit == j),11])
    } else {
      sh[i,j+1] <- "." # if not visited, gets placeholder
      sh[i,j+4] <- "9"
      sh[i,j+7] <- "9 9"
    }
  }
}


# Create cottontail capt hist (eastern and appalachian combined)
ec <- matrix(data = NA, nrow = length(unique(visits$Transect)), ncol = 10)
colnames(ec) <- c("Transect","V1","V2","V3","sc1","sc2","sc3","scPerc1","scPerc2","scPerc3")
ec[,"Transect"] <- c(unique(visits$Transect))
ec <- as.data.frame(ec)
for (i in 1:nrow(ec)) { # for each transect
  for (j in 1:3) { # for each visit
    if (j %in% visits[which(visits$Transect == ec$Transect[i]),"Visit"]) { # if the transect was visited
      # if ec or ac detected, matrix gets a 1
      # if not detected, matrix gets a 0
      ec[i,j+1] <- ifelse("Eastern Cottontail" %in% dna[which(dna$Transect == ec$Transect[i] & dna$Visit == j),"Species"],
                          "1",ifelse("Appalachian Cottontail" %in% dna[which(dna$Transect == ec$Transect[i] & dna$Visit == j),"Species"],
                                     "1","0"))
      # fill in snow cover
      ec[i,j+4] <- visits[which(visits$Transect == ec$Transect[i] & visits$Visit == j),9]
      ec[i,j+7] <- paste0(visits[which(visits$Transect == ec$Transect[i] & visits$Visit == j),10]," ",
                          visits[which(visits$Transect == ec$Transect[i] & visits$Visit == j),11])
    } else {
      ec[i,j+1] <- "."
      ec[i,j+4] <- "9"
      ec[i,j+7] <- "9 9"
    }
  }
}

# check totals
nrow(sh[which(sh$V1 == 1 | sh$V2 == 1 | sh$V3 == 1),]) # 30
nrow(ec[which(ec$V1 == 1 | ec$V2 == 1 | ec$V3 == 1),]) # 41

# Multi-species capture history for SH and EC
multi <- matrix(data = NA, nrow = length(unique(visits$Transect)), ncol = 10)
colnames(multi) <- c("Transect","V1","V2","V3","sc1","sc2","sc3","scPerc1","scPerc2","scPerc3")
multi[,"Transect"] <- c(unique(visits$Transect))
multi <- as.data.frame(multi)
for (i in 1:nrow(multi)) {
  for (j in 1:3) {
    multi[i,j+4] <- sh[i,j+4]
    multi[i,j+7] <- sh[i,j+7]
    if (sh[i,j+1] == "0" & ec[i,j+1] == "0") {
      multi[i,j+1] <- "00"
    } else if (sh[i,j+1] == "1" & ec[i,j+1] == "0") {
      multi[i,j+1] <- "01"
    } else if (sh[i,j+1] == "0" & ec[i,j+1] == "1") {
      multi[i,j+1] <- "02"
    } else if (sh[i,j+1] == "1" & ec[i,j+1] == "1") {
      multi[i,j+1] <- "03"
    } else if (sh[i,j+1] == "." & ec[i,j+1] == ".") {
      multi[i,j+1] <- "."
    }
  }
}


# Import transect covariates and coordinates
coords <- read.csv("TransectCoords.csv")
grid <- read.csv("NewStudyAreaGrid.csv")

# calculate percent of forested(landscape) that is e-s(landscape)
grid$PercESpF4500 <- grid$PercES4500/grid$PercF4500

# Join coords then covariates from grid
sh <- inner_join(sh,coords,by=join_by("Transect"=="TransectPoints.TransectID"))
sh <- sh[,c(1:10,13,14)]
colnames(sh)[11:12] <- c("Longitude","Latitude")
saveRDS(sh,"sh.RData")
multi <- inner_join(multi,coords,by=join_by("Transect"=="TransectPoints.TransectID"))
multi <- multi[,c(1:10,13,14)]
colnames(multi)[11:12] <- c("Longitude","Latitude")
ec <- inner_join(ec,coords,by=join_by("Transect"=="TransectPoints.TransectID"))
ec <- ec[,c(1:10,13,14)]
colnames(ec)[11:12] <- c("Longitude","Latitude")
ch_sh <- inner_join(sh,grid,by=join_by("Latitude","Longitude"))
ch_multi <- inner_join(multi,grid,by=join_by("Latitude","Longitude"))
ch_ec <- inner_join(ec,grid,by=join_by("Latitude","Longitude"))
saveRDS(ch_sh,"ch_sh.RData")

# Only keep needed columns
ch_sh <- ch_sh[,c(1:12,17:20,23:28)]
ch_multi <- ch_multi[,c(1:12,17:20,23:28)]
ch_ec <- ch_ec[,c(1:12,17:20,23:28)]

# standardize covariates
for (i in 13:22) {
  var <- colnames(ch_sh)[i]
  eval(parse(text=paste0("ch_sh$",var,"Stnd <- (ch_sh$",var," - mean(ch_sh$",
                         var,"))/sd(ch_sh$",var,")")))
}
saveRDS(ch_sh,"ch_sh_stnd.RData")
for (i in 13:22) {
  var <- colnames(ch_multi)[i]
  eval(parse(text=paste0("ch_multi$",var,"Stnd <- (ch_multi$",var," - mean(ch_multi$",
                         var,"))/sd(ch_multi$",var,")")))
}
for (i in 13:22) {
  var <- colnames(ch_ec)[i]
  eval(parse(text=paste0("ch_ec$",var,"Stnd <- (ch_ec$",var," - mean(ch_ec$",
                         var,"))/sd(ch_ec$",var,")")))
}
# save standardized values for back-transforming
stnd <- matrix(data = NA, ncol = 10, nrow = 2)
colnames(stnd) <- colnames(ch_sh)[13:22]
for (i in 1:ncol(stnd)) {
  eval(parse(text=paste0("stnd[1,",i,"] <- mean(ch_sh$",
                         colnames(ch_sh)[i+12],")")))
  eval(parse(text=paste0("stnd[2,",i,"] <- sd(ch_sh$",colnames(ch_sh)[i+12],")")))
}

write.csv(stnd,"Newstandardization.csv",row.names=FALSE)



# Write final sh, ec, and multi files into text files for MARK

fn <- "New_sh_only.inp"
ch_sh$tag <- paste("/*", ch_sh$Transect, "*/")
ch_sh$capthist <- paste0(ch_sh$tag," ",ch_sh$V1,ch_sh$V2,ch_sh$V3," ","1"," ",
                         ch_sh$sc1," ",ch_sh$sc2," ",ch_sh$sc3," ",ch_sh$scPerc1," ",
                         ch_sh$scPerc2," ",ch_sh$scPerc3," ",ch_sh$PercESStnd," ",
                         ch_sh$PercentConStnd," ",ch_sh$PercF4500Stnd," ",
                         ch_sh$PercES4500Stnd," ",
                         ch_sh$PercESpF4500Stnd," ",ch_sh$PercHD4500Stnd,
                         " ",ch_sh$sf2Stnd," ",ch_sh$sd2Stnd," ",
                         ch_sh$dd2Stnd," ",ch_sh$bf2Stnd,";")

write(paste("/* 2023 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Snow cover present (0 = <=25% ground cover, 1 = >25% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Percent snow cover present (00 = <25% ground cover, 10 = 25 - 75% ground cover, 01 = >75% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Proportion early-successional cover within sample unit",
            "\n   Proportion conifer cover within sample unit",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion early-successional cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion human-developed cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Average snowfall during the period 2013-2023, in cm",
            "\n   Average number of days per year with snow cover during the period 2013-2023",
            "\n   Average degree days (sum of maximum daily temperatures below 0 Celsius) during the period 2013-2023, in Celsius",
            "\n   Average number of days per year with a maximum temperature below 0 Celsius during the period 2013-2023",
            "*/"),
      file=fn)
write(ch_sh$capthist, fn, append=T)


fn <- "New_cottontail.inp"
ch_ec$tag <- paste("/*", ch_ec$Transect, "*/")
ch_ec$capthist <- paste0(ch_ec$tag," ",ch_ec$V1,ch_ec$V2,ch_ec$V3," ","1"," ",
                         ch_ec$sc1," ",ch_ec$sc2," ",ch_ec$sc3," ",ch_ec$scPerc1," ",
                         ch_ec$scPerc2," ",ch_ec$scPerc3," ",ch_ec$PercESStnd," ",
                         ch_ec$PercentConStnd," ",ch_ec$PercF4500Stnd," ",
                         ch_ec$PercES4500Stnd," ",
                         ch_ec$PercESpF4500Stnd," ",ch_ec$PercHD4500Stnd,
                         " ",ch_ec$sf2Stnd," ",ch_ec$sd2Stnd," ",
                         ch_ec$dd2Stnd," ",ch_ec$bf2Stnd,";")

write(paste("/* 2023 Cottontail Rabbit Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Snow cover present (0 = <=25% ground cover, 1 = >25% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Percent snow cover present (00 = <25% ground cover, 10 = 25 - 75% ground cover, 01 = >75% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Proportion early-successional cover within sample unit",
            "\n   Proportion conifer cover within sample unit",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion early-successional cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion human-developed cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Average snowfall during the period 2013-2023, in cm",
            "\n   Average number of days per year with snow cover during the period 2013-2023",
            "\n   Average degree days (sum of maximum daily temperatures below 0 Celsius) during the period 2013-2023, in Celsius",
            "\n   Average number of days per year with a maximum temperature below 0 Celsius during the period 2013-2023",
            "*/"),
      file=fn)
write(ch_ec$capthist, fn, append=T)

# snowshoe hare climate only
fn <- "New_sh_climate.inp"
ch_sh$tag <- paste("/*", ch_sh$Transect, "*/")
ch_sh$capthist <- paste0(ch_sh$tag," ",ch_sh$V1,ch_sh$V2,ch_sh$V3," ","1"," ",
                         ch_sh$sc1," ",ch_sh$sc2," ",ch_sh$sc3," ",ch_sh$scPerc1," ",
                         ch_sh$scPerc2," ",ch_sh$scPerc3," ",
                         " ",ch_sh$sf2Stnd," ",ch_sh$sd2Stnd," ",
                         ch_sh$dd2Stnd," ",ch_sh$bf2Stnd,";")

write(paste("/* 2023 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Snow cover present (0 = <=25% ground cover, 1 = >25% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Percent snow cover present (00 = <25% ground cover, 10 = 25 - 75% ground cover, 01 = >75% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Average snowfall during the period 2013-2023, in cm",
            "\n   Average number of days per year with snow cover during the period 2013-2023",
            "\n   Average degree days (sum of maximum daily temperatures below 0 Celsius) during the period 2013-2023, in Celsius",
            "\n   Average number of days per year with a maximum temperature below 0 Celsius during the period 2013-2023",
            "*/"),
      file=fn)
write(ch_sh$capthist, fn, append=T)

# snowshoe hare detection only
fn <- "New_sh_detection.inp"
ch_sh$tag <- paste("/*", ch_sh$Transect, "*/")
ch_sh$capthist <- paste0(ch_sh$tag," ",ch_sh$V1,ch_sh$V2,ch_sh$V3," ","1"," ",
                         ch_sh$sc1," ",ch_sh$sc2," ",ch_sh$sc3," ",ch_sh$scPerc1," ",
                         ch_sh$scPerc2," ",ch_sh$scPerc3,";")

write(paste("/* 2023 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Snow cover present (0 = <=25% ground cover, 1 = >25% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Percent snow cover present (00 = <25% ground cover, 10 = 25 - 75% ground cover, 01 = >75% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "*/"),
      file=fn)
write(ch_sh$capthist, fn, append=T)

################################################################################
# take subset with >90% forested in buffer
# then model occupancy with es and esB
# can compare effects 

sh90 <- subset(ch_sh, ch_sh$PercF4500 >= 0.9)

# restandardize
# standardize covariates
for (i in 13:23) {
  var <- colnames(sh90)[i]
  eval(parse(text=paste0("sh90$",var,"Stnd <- (sh90$",var," - mean(sh90$",
                         var,"))/sd(sh90$",var,")")))
}
# save standardized values for back-transforming
stnd90 <- matrix(data = NA, ncol = 11, nrow = 2)
colnames(stnd90) <- colnames(ch_sh)[13:23]
for (i in 1:ncol(stnd90)) {
  eval(parse(text=paste0("stnd90[1,",i,"] <- mean(ch_sh$",
                         colnames(ch_sh)[i+12],")")))
  eval(parse(text=paste0("stnd90[2,",i,"] <- sd(ch_sh$",colnames(ch_sh)[i+12],")")))
}

write.csv(stnd90,"standardization90.csv",row.names=FALSE)

fn <- "sh90.inp"
sh90$tag <- paste("/*", sh90$Transect, "*/")
sh90$capthist <- paste0(sh90$tag," ",sh90$V1,sh90$V2,sh90$V3," ","1"," ",
                        sh90$sc1," ",sh90$sc2," ",sh90$sc3," ",sh90$scPerc1," ",
                        sh90$scPerc2," ",sh90$scPerc3," ",sh90$PercESStnd," ",
                        sh90$PercentConStnd," ",sh90$PercES4500Stnd," ",
                        sh90$PercESpF4500Stnd," ",sh90$PercHD4500Stnd,
                        " ",sh90$PercF4500Stnd," ",sh90$sf2Stnd," ",sh90$sd2Stnd," ",
                        sh90$dd2Stnd," ",sh90$bf2Stnd,";")

write(paste("/* 2023 Snowshoe Hare Capture History Northern PA",
            "\n Only transects that have > 90% forested cover in the buffer",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates: ",
            "\n   Snow cover present (0 = <=25% ground cover, 1 = >25% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Percent snow cover present (00 = <25% ground cover, 10 = 25 - 75% ground cover, 01 = >75% ground cover) for each capture occasion ('9' for occasions not surveyed)",
            "\n   Proportion early-successional cover within sample unit",
            "\n   Proportion conifer cover within sample unit",
            "\n   Proportion early-successional cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion human-developed cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Average snowfall during the period 2013-2023, in cm",
            "\n   Average number of days per year with snow cover during the period 2013-2023",
            "\n   Average degree days (sum of maximum daily temperatures below 0 Celsius) during the period 2013-2023, in Celsius",
            "\n   Average number of days per years with a maximum temperature below 0 Celsius during the period 2013-2023",
            "*/"),
      file=fn)
write(sh90$capthist, fn, append=T)

###########################################################################
# Code for generating shorter testing files, fewer covariates

gridESdist <- read.csv("C:/Users/alz5215/OneDrive - The Pennsylvania State University/Documents/Research/GIS/ES_Distances/StudyAreaGrid_ESdist.csv")

ch_sh$ESdist <- c(NA)
# add ES dist values to ch_sh based on coords
for (i in 1:nrow(ch_sh)) {
  coords <- ch_sh[i,c("Latitude","Longitude")]
  dist <- gridESdist[which(gridESdist$Latitude == coords[[1]]),"NEAR_DIST"]
  ch_sh[i,"ESdist"] <- dist
}

# standardize
meanDist <- mean(ch_sh$ESdist)
sdDist <- sd(ch_sh$ESdist)
ch_sh$ESdistStnd <- (ch_sh$ESdist - mean(ch_sh$ESdist))/sd(ch_sh$ESdist)


fn <- "sh_ESdist.inp"
ch_sh$tag <- paste("/*", ch_sh$Transect, "*/")
ch_sh$capthist <- paste0(ch_sh$tag," ",ch_sh$V1,ch_sh$V2,ch_sh$V3," ","1"," ",
                         ch_sh$PercF4500Stnd," ",ch_sh$PercESpF4500Stnd," ",
                         ch_sh$ESdistStnd," ",ch_sh$sf2Stnd,";")

write(paste("/* 2023 Snowshoe Hare Capture History Northern PA",
            "\n  Transect ID",
            "\n  3 capture occasions ('.' for occasions not surveyed)",
            "\n  Covariates (standardized): ",
            "\n   Proportion forested cover within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Proportion of forested cover that is early-successional within a square-shaped buffer around the sample unit with a side length of 9 km",
            "\n   Distance to nearest patch of early-successional cover >5 ha in area present in 2004",
            "\n   Average snowfall during the period 2013-2023, in cm",
            "*/"),
      file=fn)
write(ch_sh$capthist, fn, append=T)
