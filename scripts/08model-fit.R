################################################################################
# Assessing occupancy model fit
# Amanda Zak
# March 2024
################################################################################
library(dplyr)
library(stringr)
library(tidyr)
library(RMark)
library(truncnorm)

# Using the Mackenzie-Bailey approach (Mackenzie et al. 2017 book, page 157)
# With modification - all revisits will be assigned to "1.." capt hist for calculating stats
# Because the "1.." will be biased low due to them being re-surveyed and moved to other cohorts

# Load 2023 capture history
ch_sh <- readRDS("output/ch_sh.RData")
# need to reassign transect 287 from .0. to 0..
ch_sh[which(ch_sh$Transect == 287),c(2,3)] <- c("0",".")
# create single capt hist column
ch_sh$ch <- paste(ch_sh$V1,ch_sh$V2,ch_sh$V3,sep="")

# reassign the "01." capthists to "10." so that they correctly collapse down to "occupied" later
for (row in 1:nrow(ch_sh)) {
  if (ch_sh$ch[row] == "01.") {
    ch_sh$ch[row] <- "10."
  }
}

# number of "0" sites
#num0 <- sum(startsWith(ch_sh$ch,"0"))
# number of "1" sites
#num1 <- sum(startsWith(ch_sh$ch,"1"))

# list of possible capture histories, broken into cohorts
cohort <- c("1","2","3")
c1 <- c("0..","1..")
c2 <- c("00.","11.","10.","01.")
c3 <- c("111","100","010","001","110","101","011","000")
all_c <- c(c1,c2,c3)

# Find how many instances of each observed capthist
rows <- length(all_c) # number of unique capt hists possible
ohc <- as.data.frame(matrix(nrow = rows, ncol = 3))
ohc[,1] <- c(1,1,2,2,2,2,3,3,3,3,3,3,3,3) # fill in cohort numbers
ohc[,2] <- all_c # fill in capt hists
for (x in 1:rows) {
  num_obs <- nrow(ch_sh[which(ch_sh$ch == ohc[x,2]),])
  ohc[x,3] <- num_obs
}


# Get covariate values for each transect
stnd <- read.csv("data/standardization2023.csv",row.names = NULL)
# standardize
for (i in c(17,25,28)) {
  var <- colnames(ch_sh)[i]
  eval(parse(text=paste0("ch_sh$",var,"Stnd <- (ch_sh$",var," - mean(ch_sh$",
                         var,"))/sd(ch_sh$",var,")")))
}

# Function to calculate occupancy probability *for 2023*, back-transformed
a <- -2.1594519
sfb <- 0.8779362
fb <- 0.7056529
esb <- 1.6499307
calcOccProb <- function(sfStnd,fbStnd,esbStnd) { #function takes standardized values
  logOcc <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd # calculate occupancy probability on logit scale
  realOcc <- exp(logOcc)/(1+exp(logOcc)) # back-transforms to real probability
  return(realOcc)
}

# Detection probability
p <- 0.723

# Create a data frame with the probability formulas for each capt hist
ch <- c("0..","1..","00.","11.","10.","01.","000","100","010","001","110","101","011","111")
probs <- c("(1-psi) + (psi*(1-p))","(psi*p)",
           "(1-psi) + (psi*(1-p)*(1-p))","(psi*p*p)","(psi*p*(1-p))","(psi*(1-p)*p)",
           "(1-psi) + (psi*(1-p)*(1-p)*(1-p))","(psi*p*(1-p)*(1-p))","(psi*(1-p)*p*(1-p))",
           "(psi*(1-p)*(1-p)*p)","(psi*p*p*(1-p))","(psi*p*(1-p)*p)","(psi*(1-p)*p*p)",
           "(psi*p*p*p)")
possibleCH <- as.data.frame(ch) %>% cbind(probs)

# calculate Eh for each capthist/cohort (expected number of units w/ a specific capt hist)

# empty array to hold ehc
ehc <- as.data.frame(matrix(nrow = rows, ncol = 3))

# cycle through each cohort
for (i in 1:length(cohort)) {
  
  # get the sites that are in the cohort
  sc <- eval(parse(text=paste("ch_sh[which(ch_sh$ch %in% c",i,"),]",sep=""))) 
  
  # get sf, f, and esf values for the sc sites (standardized)
  sf <- sc$sf2Stnd
  f <- sc$PercF4500Stnd
  esf <- sc$PercESpF4500Stnd
  
  # get a vector of occupancy probabilities for each site in sc
  psi <- calcOccProb(sf,f,esf)
  
  # set up an empty vector to hold expected number of sites for all capt hists in the cohort
  eh <- c()
  
  # cycle through each capthist
  for (j in 1:eval(parse(text=paste("length(c",i,")",sep="")))) { 
    
    # get capt hist text string
    eval(parse(text=paste("string <- c",i,"[",j,"]",sep="")))
    
    # find the matching probability string
    Pr <- possibleCH[which(possibleCH$ch == string),2]
    
    # evaluate the probability string across all sites in the cohort
    prob <- eval(parse(text = paste(Pr)))
    
    # calculate the expected number of sites with that capt hist
    eh[j] <- sum(prob)
    
    # save these values
    iteration <- ifelse(i < 3,ifelse(i== 1,i*j,i+j),(i*2)+j)
    ehc[iteration,1] <- i
    ehc[iteration,2] <- string
    ehc[iteration,3] <- sum(prob)
    
  }
}

# adjustments:
# update ohc into only 0s and 1s
ohc0 <- sum(ohc$V3[which(startsWith(ohc$V2,"0"))])
ohc1 <- sum(ohc$V3[which(startsWith(ohc$V2,"1"))])
# update ehc into only 0s and 1s
ehc0 <- sum(ehc$V3[which(startsWith(ehc$V2,"0"))])
ehc1 <- sum(ehc$V3[which(startsWith(ehc$V2,"1"))])

# Calculate chi-square test statistic JUST ON SINGLE-OCCASION SITES
ohc_vec <- c(ohc0,ohc1)
ehc_vec <- c(ehc0,ehc1)
test <- ((ohc_vec-ehc_vec)^2)/ehc_vec
test_stat <- sum(test)

# Save separate from results of bootstrap
chi_sq_obs <- test_stat
saveRDS(test_stat,"output/chi_sq_obs.RData")


### Bootstrapping

# Subset ch_sh for this step
ch_boot_orig <- ch_sh[,c(1,29,30,31,32)]

# Add column with number of occasions
ch_boot_orig$occasions <- nchar(gsub(".","",ch_boot_orig$ch,fixed=TRUE))

# Empty list to hold chi-square statistics from each iteration
chi_sq_boot <- c()

nrep <- 1000
try(
  for (iter in 1:nrep) {
    
    # Make a copy of ch_boot_orig to use in the loop, will rewrite for every iteration
    ch_boot <- ch_boot_orig
    
    # Reinitialize p at the start of ever iteration
    p <- 0.723
    
    # 1. Generate capture history
    
    # generate occupancy state for the cell based on psi
    ch_boot$psi <- calcOccProb(ch_boot$sf2Stnd,ch_boot$PercF4500Stnd,ch_boot$PercESpF4500Stnd)
    rand_num <- runif(nrow(ch_boot))
    ch_boot$occ_state <- ifelse(rand_num <= ch_boot$psi,1,0)
    
    # add a capthist column
    ch_boot$capthist <- NA
    
    # go through each site and generate detection histories for each occasion
    for (j in 1:nrow(ch_boot)) {
      
      # empty list to hold outcomes for each occasion
      ch_list <- c()
      
      # for first occasion, sample all sites:
      if (ch_boot$occ_state[j] == 0) {
        
        # if site is not occupied, will get 0..
        ch_list[1] <- "0.."
        
      } else {
        
        # if site is occupied, need to generate detection
        rand_num <- runif(1) # generate random number
        det <- ifelse(rand_num <= p,"1","0..")
        ch_list[1] <- det
      }
      
      # If site had a 1, it will get a second occasion
        # With a 50% chance it will get a third occasion
      
      if (ch_list[1] == "1") {
        
        # Second occasion
        rand_num <- runif(1)
        det <- as.character(ifelse(rand_num <= p,1,0))
        ch_list[2] <- det
        
        # Third occasion
        rand_num <- runif(1)
        if (rand_num < 0.5) { # 50% chance of being surveyed a 3rd time
          rand_num <- runif(1)
          det <- as.character(ifelse(rand_num <= p,1,0))
          ch_list[3] <- det
        } else {
          ch_list[3] <- "."
        }
      }
      
      # Write the ch_list to the capthist
      ch_boot$capthist[j] <- paste0(ch_list,collapse="")
    }
    
    # 2. Fit occupancy model and save model coefficients
    
    ch_boot_mark <- ch_boot[,c(9,3,4,5)]
    ch_boot_mark$group <- rep(1)
    colnames(ch_boot_mark) <- c("ch","sf","forest","es","group")
    
    output <- mark(ch_boot_mark,model="Occupancy",
                   model.parameters=list(Psi=list(formula=~sf+forest+es),
                                         p=list(formula=~1)), 
                   output = T, invisible = T, delete = T)
    
    # Occupancy coefficients
    a_boot <- output$results$beta[2,1]
    sfb_boot <- output$results$beta[3,1]
    fb_boot <- output$results$beta[4,1]
    esb_boot <- output$results$beta[5,1]
    calcOccProb_boot <- function(sfStnd,fbStnd,esbStnd) { #function takes standardized values
      logOcc <- a_boot + sfb_boot*sfStnd + fb_boot*fbStnd + esb_boot*esbStnd # calculate occupancy probability on logit scale
      realOcc <- exp(logOcc)/(1+exp(logOcc)) # back-transforms to real probability
      return(realOcc)
    }
    
    # Detection probability
    p_log <- output$results$beta[1,1]
    p <- exp(p_log)/(1+exp(p_log))
    
    # 3. Calculate expected capture histories
    
    # observed capture histories
    rows <- length(all_c) # number of unique capt hists possible
    ohc_boot <- as.data.frame(matrix(nrow = rows, ncol = 3))
    ohc_boot[,1] <- c(1,1,2,2,2,2,3,3,3,3,3,3,3,3) # fill in cohort numbers
    ohc_boot[,2] <- all_c # fill in capt hists
    for (x in 1:rows) {
      num_obs <- nrow(ch_boot[which(ch_boot$capthist == ohc_boot[x,2]),])
      ohc_boot[x,3] <- num_obs
    }
    
    ehc_boot <- as.data.frame(matrix(nrow = rows, ncol = 3))
    
    for (i in 1:length(cohort)) {
      
      # get the sites that are in the cohort
      sc <- eval(parse(text=paste("ch_boot[which(ch_boot$capthist %in% c",i,"),]",sep=""))) 
      
      # get sf, f, and esf values for the sc sites (standardized)
      sf <- sc$sf2Stnd
      f <- sc$PercF4500Stnd
      esf <- sc$PercESpF4500Stnd
      
      # get a vector of occupancy probabilities for each site in sc
      psi <- calcOccProb_boot(sf,f,esf)
      
      # cycle through each capthist
      for (j in 1:eval(parse(text=paste("length(c",i,")",sep="")))) { 
        
        # get capt hist text string
        eval(parse(text=paste("string <- c",i,"[",j,"]",sep="")))
        
        # find the matching probability string
        Pr <- possibleCH[which(possibleCH$ch == string),2]
        
        # evaluate the probability string for each site
        prob <- eval(parse(text = paste(Pr)))
        
        # save these values
        iteration <- ifelse(i < 3,ifelse(i== 1,i*j,i+j),(i*2)+j)
        ehc_boot[iteration,1] <- i
        ehc_boot[iteration,2] <- string
        ehc_boot[iteration,3] <- sum(prob)
        
      }
    }
    
    # Adjustments:
    # update ohc into only 0s and 1s
    ohc0 <- sum(ohc_boot$V3[which(startsWith(ohc_boot$V2,"0"))])
    ohc1 <- sum(ohc_boot$V3[which(startsWith(ohc_boot$V2,"1"))])
    # update ehc into only 0s and 1s
    ehc0 <- sum(ehc_boot$V3[which(startsWith(ehc_boot$V2,"0"))])
    ehc1 <- sum(ehc_boot$V3[which(startsWith(ehc_boot$V2,"1"))])
    
    # 4. Calculate chi-square test statistic JUST ON SINGLE-OCCASION SITES
    ohc_vec <- c(ohc0,ohc1)
    ehc_vec <- c(ehc0,ehc1)
    test <- ((ohc_vec-ehc_vec)^2)/ehc_vec
    test_stat_boot <- sum(test)
    
    
    # 5. Save to list
    chi_sq_boot[iter] <- test_stat_boot
  }
)

saveRDS(chi_sq_boot, "output/chi_sq_boot.RData")

# get mean and sd of chi_sq_boot
chi_mn <- mean(chi_sq_boot)
chi_sd <- sd(chi_sq_boot)
# plot a truncated normal dist over it to see if it matches
x <- seq(0,2,0.01)
hist(chi_sq_boot,freq=F,breaks=200)
abline(v = chi_sq_obs, col = "red")

# p-value will be proportion of chi-sq-stats above the chi-sq-obs
length(which(chi_sq_boot > chi_sq_obs))/length(chi_sq_boot)
# 0.704 for 1000 iter

################################################################################
# 2004 - modified code ----------------------------------------------------

capthist04 <- as_tibble(read.delim("data/2004capthist.inp"))[9:246,] %>% 
  separate_wider_delim(1, delim = " ", too_many = "debug",
                       names = c("id","id2","ch","group","es","esb","esfb","fb","sf")) %>% 
  separate_wider_delim(1, delim = "*", names = c("id1","id")) %>% 
  separate_wider_delim(10, delim = ";", names = c("sf","discard"))
capthist04 <- capthist04[,c(2,4:10)]

# list of possible capture histories, broken into cohorts
cohort <- c("1","2","3")
c1 <- c("0..","1..")
c2 <- c("00.","11.","10.","01.")
c3 <- c("111","100","010","001","110","101","011","000")
all_c <- c(c1,c2,c3)

# Find how many instances of each observed capthist
rows <- length(all_c) # number of unique capt hists possible
ohc <- matrix(nrow = rows, ncol = 3)
ohc[,1] <- c(1,1,2,2,2,2,3,3,3,3,3,3,3,3) # fill in cohort numbers
ohc[,2] <- all_c # fill in capt hists
for (x in 1:rows) {
  num_obs <- nrow(capthist04[which(capthist04$ch == ohc[x,2]),])
  ohc[x,3] <- num_obs
}

# Function to calculate occupancy probability *for 2004*, back-transformed
a <- -0.0985304
sfb <- 0.4419514
fb <- -0.0156494
esb <- 0.7554851
calcOccProb <- function(sfStnd,fbStnd,esbStnd) { #function takes standardized values
  logOcc <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd # calculate occupancy probability on logit scale
  realOcc <- exp(logOcc)/(1+exp(logOcc)) # back-transforms to real probability
  return(realOcc)
}

# Detection probability
p <- 0.609

# Create a data frame with the probability formulas for each capt hist
ch <- c("0..","1..","00.","11.","10.","01.","000","100","010","001","110","101","011","111")
probs <- c("(1-psi) + (psi*(1-p))","(psi*p)",
           "(1-psi) + (psi*(1-p)*(1-p))","(psi*p*p)","(psi*p*(1-p))","(psi*(1-p)*p)",
           "(1-psi) + (psi*(1-p)*(1-p)*(1-p))","(psi*p*(1-p)*(1-p))","(psi*(1-p)*p*(1-p))",
           "(psi*(1-p)*(1-p)*p)","(psi*p*p*(1-p))","(psi*p*(1-p)*p)","(psi*(1-p)*p*p)",
           "(psi*p*p*p)")
possibleCH <- as.data.frame(ch) %>% cbind(probs)

# calculate Eh for each capthist/cohort (expected number of units w/ a specific capt hist)

# empty array to hold ehc
ehc <- matrix(nrow = rows, ncol = 3)

# cycle through each cohort
for (i in 1:length(cohort)) {
  
  # get the sites that are in the cohort
  sc <- eval(parse(text=paste("capthist04[which(capthist04$ch %in% c",i,"),]",sep=""))) 
  
  # get sf, f, and esf values for the sc sites (standardized)
  sf <- as.numeric(sc$sf)
  f <- as.numeric(sc$fb)
  esf <- as.numeric(sc$esfb)
  
  # get a vector of occupancy probabilities for each site in sc
  psi <- calcOccProb(sf,f,esf)
  
  # set up an empty vector to hold each eh
  eh <- c()
  
  # cycle through each capthist
  for (j in 1:eval(parse(text=paste("length(c",i,")",sep="")))) { 
    
    # get capt hist text string
    eval(parse(text=paste("string <- c",i,"[",j,"]",sep="")))
    
    # find the matching probability string
    Pr <- possibleCH[which(possibleCH$ch == string),2]
    
    # evaluate the probability string for each site
    prob <- eval(parse(text = paste(Pr)))
    
    # calculate the expected number of sites with that capt hist
    eh[j] <- sum(prob)
    
    # save these values
    iteration <- ifelse(i < 3,ifelse(i== 1,i*j,i+j),(i*2)+j)
    ehc[iteration,1] <- i
    ehc[iteration,2] <- string
    ehc[iteration,3] <- sum(prob)
    
  }
}

# adjustments:
ohc <- as.data.frame(ohc)
ohc$V3 <- as.numeric(ohc$V3)
ehc <- as.data.frame(ehc)
ehc$V3 <- as.numeric(ehc$V3)
# update ohc into only 0s and 1s
ohc0 <- sum(ohc$V3[which(startsWith(ohc$V2,"0"))])
ohc1 <- sum(ohc$V3[which(startsWith(ohc$V2,"1"))])
# update ehc into only 0s and 1s
ehc0 <- sum(ehc$V3[which(startsWith(ehc$V2,"0"))])
ehc1 <- sum(ehc$V3[which(startsWith(ehc$V2,"1"))])

# Calculate chi-square test statistic JUST ON SINGLE-OCCASION SITES
ohc_vec <- c(ohc0,ohc1)
ehc_vec <- c(ehc0,ehc1)
test <- ((ohc_vec-ehc_vec)^2)/ehc_vec
test_stat <- sum(test)

# Save separate from results of bootstrap
chi_sq_obs04 <- test_stat
saveRDS(test_stat,"output/chi_sq_obs04.RData")

### Bootstrapping


# Add column with number of occasions
capthist04$occasions <- nchar(gsub(".","",capthist04$ch,fixed=TRUE))

# Empty list to hold chi-square statistics from each iteration
chi_sq_boot04 <- c()

nrep <- 1000
for (iter in 1:nrep) {
  
  # Make a copy of capthist04 to use in the loop, will rewrite for every iteration
  ch_boot <- capthist04
  
  # Reinitialize p at the start of ever iteration
  p <- 0.609
  
  # 1. Generate capture history
  
  # generate occupancy state for the cell based on psi
  ch_boot$psi <- calcOccProb(as.numeric(ch_boot$sf),as.numeric(ch_boot$fb),as.numeric(ch_boot$esfb))
  ch_boot$rand_num <- runif(nrow(ch_boot))
  ch_boot$occ_state <- ifelse(ch_boot$rand_num <= ch_boot$psi,1,0)
  
  # add a capthist column
  ch_boot$capthist <- NA
  
  # go through each site and generate detection histories for each occasion
  for (j in 1:nrow(ch_boot)) {
    
    # empty list to hold outcomes for each occasion
    ch_list <- c()
    
    # for first occasion, sample all sites:
    if (ch_boot$occ_state[j] == 0) {
      
      # if site is not occupied, will get 0..
      ch_list[1] <- "0.."
      
    } else {
      
      # if site is occupied, need to generate detection
      rand_num <- runif(1) # generate random number
      det <- ifelse(rand_num <= p,"1","0..")
      ch_list[1] <- det
    }
    
    # If site had a 1, it will get a second occasion
    # With a 50% chance it will get a third occasion
    
    if (ch_list[1] == "1") {
      
      # Second occasion
      rand_num <- runif(1)
      det <- as.character(ifelse(rand_num <= p,1,0))
      ch_list[2] <- det
      
      # Third occasion
      rand_num <- runif(1)
      if (rand_num < 0.5) { # 50% chance of being surveyed a 3rd time
        rand_num <- runif(1)
        det <- as.character(ifelse(rand_num <= p,1,0))
        ch_list[3] <- det
      } else {
        ch_list[3] <- "."
      }
    }
    
    # Write the ch_list to the capthist
    ch_boot$capthist[j] <- paste0(ch_list,collapse="")
  }
  
  # 2. Fit occupancy model and save model coefficients
  
  ch_boot_mark <- ch_boot[,c(2,8,7,6,3)]
  colnames(ch_boot_mark) <- c("ch","sf","forest","es","group")
  ch_boot_mark[,2:4] <- lapply(ch_boot_mark[,2:4],as.numeric)
  
  output <- mark(ch_boot_mark,model="Occupancy",
                 model.parameters=list(Psi=list(formula=~sf+forest+es),
                                       p=list(formula=~1)), 
                 output = F, invisible = T, delete = T)
  
  # Occupancy coefficients
  a_boot <- output$results$beta[2,1]
  sfb_boot <- output$results$beta[3,1]
  fb_boot <- output$results$beta[4,1]
  esb_boot <- output$results$beta[5,1]
  calcOccProb <- function(sfStnd,fbStnd,esbStnd) { #function takes standardized values
    logOcc <- a_boot + sfb_boot*sfStnd + fb_boot*fbStnd + esb_boot*esbStnd # calculate occupancy probability on logit scale
    realOcc <- exp(logOcc)/(1+exp(logOcc)) # back-transforms to real probability
    return(realOcc)
  }
  
  # Detection probability
  p_log <- output$results$beta[1,1]
  p <- exp(p_log)/(1+exp(p_log))
  
  # 3. Calculate expected capture histories
  
  # observed capture histories
  rows <- length(all_c) # number of unique capt hists possible
  ohc_boot <- matrix(nrow = rows, ncol = 3)
  ohc_boot[,1] <- c(1,1,2,2,2,2,3,3,3,3,3,3,3,3) # fill in cohort numbers
  #ohc_boot[,2] <- gsub(".","",all_c,fixed=TRUE) # fill in capt hists
  ohc_boot[,2] <- all_c
  for (x in 1:rows) {
    num_obs <- nrow(ch_boot[which(ch_boot$capthist == ohc_boot[x,2]),])
    ohc_boot[x,3] <- num_obs
  }
  
  ehc_boot <- matrix(nrow = rows, ncol = 3)
  
  for (i in 1:length(cohort)) {
    
    # get the sites that are in the cohort
    sc <- eval(parse(text=paste("ch_boot[which(ch_boot$capthist %in% c",i,"),]",sep=""))) 
    
    # get sf, f, and esf values for the sc sites (standardized)
    sf <- as.numeric(sc$sf)
    f <- as.numeric(sc$fb)
    esf <- as.numeric(sc$esfb)
    
    # get a vector of occupancy probabilities for each site in sc
    psi <- calcOccProb(sf,f,esf)
    
    # cycle through each capthist
    for (j in 1:eval(parse(text=paste("length(c",i,")",sep="")))) { 
      
      # get capt hist text string
      eval(parse(text=paste("string <- c",i,"[",j,"]",sep="")))
      
      # find the matching probability string
      Pr <- possibleCH[which(possibleCH$ch == string),2]
      
      # evaluate the probability string for each site
      prob <- eval(parse(text = paste(Pr)))
      
      # save these values
      iteration <- ifelse(i < 3,ifelse(i== 1,i*j,i+j),(i*2)+j)
      ehc_boot[iteration,1] <- i
      ehc_boot[iteration,2] <- string
      ehc_boot[iteration,3] <- sum(prob)
      
    }
  }

  # Adjustments:
  ohc_boot <- as.data.frame(ohc_boot)
  ohc_boot$V3 <- as.numeric(ohc_boot$V3)
  ehc_boot <- as.data.frame(ehc_boot)
  ehc_boot$V3 <- as.numeric(ehc_boot$V3)
  # update ohc into only 0s and 1s
  ohc0 <- sum(ohc_boot$V3[which(startsWith(ohc_boot$V2,"0"))])
  ohc1 <- sum(ohc_boot$V3[which(startsWith(ohc_boot$V2,"1"))])
  # update ehc into only 0s and 1s
  ehc0 <- sum(ehc_boot$V3[which(startsWith(ehc_boot$V2,"0"))])
  ehc1 <- sum(ehc_boot$V3[which(startsWith(ehc_boot$V2,"1"))])
  
  # 4. Calculate chi-square test statistic JUST ON SINGLE-OCCASION SITES
  ohc_vec <- c(ohc0,ohc1)
  ehc_vec <- c(ehc0,ehc1)
  test <- ((ohc_vec-ehc_vec)^2)/ehc_vec
  test_stat_boot <- sum(test)
  
  # 5. Save to list
  chi_sq_boot04[iter] <- test_stat_boot
}

saveRDS(chi_sq_boot04, "output/chi_sq_boot04.RData")

# get mean and sd of chi_sq_boot
chi_mn <- mean(chi_sq_boot04)
chi_sd <- sd(chi_sq_boot04)
# plot a truncated normal dist over it to see if it matches
x <- seq(0,20,0.1)
hist(chi_sq_boot04,freq=F,breaks=50)
abline(v = chi_sq_obs04, col = "red")

# p-value will be proportion of chi-sq-stats above the chi-sq-obs
length(which(chi_sq_boot04 > chi_sq_obs04))/length(chi_sq_boot04)
# 0.943 1000

###############################################################################
# 2004 model fit to 2023 data ---------------------------------------------

# Calculate OHC with 2023 capt hist

# Load 2023 capture history
ch_sh <- readRDS("output/ch_sh.RData")
# need to reassign transect 287 from .0. to 0..
ch_sh[which(ch_sh$Transect == 287),c(2,3)] <- c("0",".")
# create single capt hist column
ch_sh$ch <- paste(ch_sh$V1,ch_sh$V2,ch_sh$V3,sep="")

# reassign the "01." capthists to "10." so that they correctly collapse down to "occupied" later
for (row in 1:nrow(ch_sh)) {
  if (ch_sh$ch[row] == "01.") {
    ch_sh$ch[row] <- "10."
  }
}

# list of possible capture histories, broken into cohorts
cohort <- c("1","2","3")
c1 <- c("0..","1..")
c2 <- c("00.","11.","10.","01.")
c3 <- c("111","100","010","001","110","101","011","000")
all_c <- c(c1,c2,c3)

# Find how many instances of each observed capthist
rows <- length(all_c) # number of unique capt hists possible
ohc <- as.data.frame(matrix(nrow = rows, ncol = 3))
ohc[,1] <- c(1,1,2,2,2,2,3,3,3,3,3,3,3,3) # fill in cohort numbers
ohc[,2] <- all_c # fill in capt hists
for (x in 1:rows) {
  num_obs <- nrow(ch_sh[which(ch_sh$ch == ohc[x,2]),])
  ohc[x,3] <- num_obs
}

# Get covariate values for each transect
# This time need to standardize 2023 values onto **2004 SCALE**
stnd04 <- read.csv("data/standardization2004.csv",row.names = NULL)
# standardize
ch_sh$sf2Stnd <- (ch_sh$sf2 - stnd04$sf1[1])/(stnd04$sf1[2])
ch_sh$PercF4500Stnd <- (ch_sh$PercF4500 - stnd04$PercF4500[1])/(stnd04$PercF4500[2])
ch_sh$PercESpF4500Stnd <- (ch_sh$PercESpF4500 - stnd04$PercESpF4500[1])/(stnd04$PercESpF4500[2])


# Calculate EHC from 2004 model

# Function to calculate occupancy probability *for 2004*, back-transformed
a <- -0.0985304
sfb <- 0.4419514
fb <- -0.0156494
esb <- 0.7554851
calcOccProb <- function(sfStnd,fbStnd,esbStnd) { #function takes standardized values
  logOcc <- a + sfb*sfStnd + fb*fbStnd + esb*esbStnd # calculate occupancy probability on logit scale
  realOcc <- exp(logOcc)/(1+exp(logOcc)) # back-transforms to real probability
  return(realOcc)
}

# Detection probability
p <- 0.609

# Create a data frame with the probability formulas for each capt hist
ch <- c("0..","1..","00.","11.","10.","01.","000","100","010","001","110","101","011","111")
probs <- c("(1-psi) + (psi*(1-p))","(psi*p)",
           "(1-psi) + (psi*(1-p)*(1-p))","(psi*p*p)","(psi*p*(1-p))","(psi*(1-p)*p)",
           "(1-psi) + (psi*(1-p)*(1-p)*(1-p))","(psi*p*(1-p)*(1-p))","(psi*(1-p)*p*(1-p))",
           "(psi*(1-p)*(1-p)*p)","(psi*p*p*(1-p))","(psi*p*(1-p)*p)","(psi*(1-p)*p*p)",
           "(psi*p*p*p)")
possibleCH <- as.data.frame(ch) %>% cbind(probs)

# calculate Eh for each capthist/cohort (expected number of units w/ a specific capt hist)

# empty array to hold ehc
ehc <- matrix(nrow = rows, ncol = 3)

# cycle through each cohort
for (i in 1:length(cohort)) {
  
  # get the sites that are in the cohort
  sc <- eval(parse(text=paste("ch_sh[which(ch_sh$ch %in% c",i,"),]",sep="")))
  
  # get sf, f, and esf values for the sc sites (standardized)
  sf <- as.numeric(sc$sf2Stnd)
  f <- as.numeric(sc$PercF4500Stnd)
  esf <- as.numeric(sc$PercESpF4500Stnd)
  
  # get a vector of occupancy probabilities for each site in sc
  psi <- calcOccProb(sf,f,esf)
  
  # set up an empty vector to hold each eh
  eh <- c()
  
  # cycle through each capthist
  for (j in 1:eval(parse(text=paste("length(c",i,")",sep="")))) { 
    
    # get capt hist text string
    eval(parse(text=paste("string <- c",i,"[",j,"]",sep="")))
    
    # find the matching probability string
    Pr <- possibleCH[which(possibleCH$ch == string),2]
    
    # evaluate the probability string for each site
    prob <- eval(parse(text = paste(Pr)))
    
    # calculate the expected number of sites with that capt hist
    eh[j] <- sum(prob)
    
    # save these values
    iteration <- ifelse(i < 3,ifelse(i== 1,i*j,i+j),(i*2)+j)
    ehc[iteration,1] <- i
    ehc[iteration,2] <- string
    ehc[iteration,3] <- sum(prob)
    
  }
}

# adjustments:
ohc <- as.data.frame(ohc)
ohc$V3 <- as.numeric(ohc$V3)
ehc <- as.data.frame(ehc)
ehc$V3 <- as.numeric(ehc$V3)
# update ohc into only 0s and 1s
ohc0 <- sum(ohc$V3[which(startsWith(ohc$V2,"0"))])
ohc1 <- sum(ohc$V3[which(startsWith(ohc$V2,"1"))])
# update ehc into only 0s and 1s
ehc0 <- sum(ehc$V3[which(startsWith(ehc$V2,"0"))])
ehc1 <- sum(ehc$V3[which(startsWith(ehc$V2,"1"))])

# Calculate chi-square test statistic JUST ON SINGLE-OCCASION SITES
ohc_vec <- c(ohc0,ohc1)
ehc_vec <- c(ehc0,ehc1)
test <- ((ohc_vec-ehc_vec)^2)/ehc_vec
test_stat <- sum(test)

chi_sq_obs04_23 <- test_stat
saveRDS(chi_sq_obs04_23, "output/chi_sq_obs04_23.RData")

# Simulate populations under 2004 model

# Add column with number of occasions
ch_sh$occasions <- nchar(gsub(".","",ch_sh$ch,fixed=TRUE))

# Empty list to hold chi-square statistics from each iteration
chi_sq_boot04_23 <- c()

nrep <- 1000
for (iter in 1:nrep) {
  
  # Make a copy of capthist04 to use in the loop, will rewrite for every iteration
  ch_boot <- ch_sh
  
  # Reinitialize p at the start of every iteration
  p <- 0.609

  # 1. Generate capture history
  
  # generate occupancy state for the cell based on psi
  ch_boot$psi <- calcOccProb(as.numeric(ch_boot$sf2Stnd),as.numeric(ch_boot$PercF4500Stnd),as.numeric(ch_boot$PercESpF4500Stnd))
  ch_boot$rand_num <- runif(nrow(ch_boot))
  ch_boot$occ_state <- ifelse(ch_boot$rand_num <= ch_boot$psi,1,0)
  
  # add a capthist column
  ch_boot$capthist <- NA
  
  # go through each site and generate detection histories for each occasion
  for (j in 1:nrow(ch_boot)) {
    
    # empty list to hold outcomes for each occasion
    ch_list <- c()
    
    # for first occasion, sample all sites:
    if (ch_boot$occ_state[j] == 0) {
      
      # if site is not occupied, will get 0..
      ch_list[1] <- "0.."
      
    } else {
      
      # if site is occupied, need to generate detection
      rand_num <- runif(1) # generate random number
      det <- ifelse(rand_num <= p,"1","0..")
      ch_list[1] <- det
    }
    
    # If site had a 1, it will get a second occasion
    # With a 50% chance it will get a third occasion
    
    if (ch_list[1] == "1") {
      
      # Second occasion
      rand_num <- runif(1)
      det <- as.character(ifelse(rand_num <= p,1,0))
      ch_list[2] <- det
      
      # Third occasion
      rand_num <- runif(1)
      if (rand_num < 0.5) { # 50% chance of being surveyed a 3rd time
        rand_num <- runif(1)
        det <- as.character(ifelse(rand_num <= p,1,0))
        ch_list[3] <- det
      } else {
        ch_list[3] <- "."
      }
    }
    
    # Write the ch_list to the capthist
    ch_boot$capthist[j] <- paste0(ch_list,collapse="")
  }
  
  # 2. Fit occupancy model and save model coefficients
  
  ch_boot_mark <- ch_boot[,c("capthist","sf2Stnd","PercF4500Stnd","PercESpF4500Stnd")]
  ch_boot_mark$group <- rep(1)
  colnames(ch_boot_mark) <- c("ch","sf","forest","es","group")
  ch_boot_mark[,2:4] <- lapply(ch_boot_mark[,2:4],as.numeric)
  
  output <- mark(ch_boot_mark,model="Occupancy",
                 model.parameters=list(Psi=list(formula=~sf+forest+es),
                                       p=list(formula=~1)), 
                 output = F, invisible = T, delete = T,
                 retry = 3) # will retry up to 3 times for inestimable parameters
  
  # Occupancy coefficients
  a_boot <- output$results$beta[2,1]
  sfb_boot <- output$results$beta[3,1]
  fb_boot <- output$results$beta[4,1]
  esb_boot <- output$results$beta[5,1]
  calcOccProb_boot <- function(sfStnd,fbStnd,esbStnd) { #function takes standardized values
    logOcc <- a_boot + sfb_boot*sfStnd + fb_boot*fbStnd + esb_boot*esbStnd # calculate occupancy probability on logit scale
    realOcc <- exp(logOcc)/(1+exp(logOcc)) # back-transforms to real probability
    return(realOcc)
  }
  
  # Detection probability
  p_log <- output$results$beta[1,1]
  p <- exp(p_log)/(1+exp(p_log))
  
  # 3. Calculate expected capture histories
  
  # observed capture histories
  rows <- length(all_c) # number of unique capt hists possible
  ohc_boot <- matrix(nrow = rows, ncol = 3)
  ohc_boot[,1] <- c(1,1,2,2,2,2,3,3,3,3,3,3,3,3) # fill in cohort numbers
  ohc_boot[,2] <- all_c
  for (x in 1:rows) {
    num_obs <- nrow(ch_boot[which(ch_boot$capthist == ohc_boot[x,2]),])
    ohc_boot[x,3] <- num_obs
  }
  
  ehc_boot <- matrix(nrow = rows, ncol = 3)
  
  for (i in 1:length(cohort)) {
    
    # get the sites that are in the cohort
    sc <- eval(parse(text=paste("ch_boot[which(ch_boot$capthist %in% c",i,"),]",sep="")))
    
    # get sf, f, and esf values for the sc sites (standardized)
    sf <- as.numeric(sc$sf2Stnd)
    f <- as.numeric(sc$PercF4500Stnd)
    esf <- as.numeric(sc$PercESpF4500Stnd)
    
    # get a vector of occupancy probabilities for each site in sc
    psi <- calcOccProb_boot(sf,f,esf)
    
    # cycle through each capthist
    for (j in 1:eval(parse(text=paste("length(c",i,")",sep="")))) { 
      
      # get capt hist text string
      eval(parse(text=paste("string <- c",i,"[",j,"]",sep="")))
      
      # find the matching probability string
      Pr <- possibleCH[which(possibleCH$ch == string),2]
      
      # evaluate the probability string for each site
      prob <- eval(parse(text = paste(Pr)))
      
      # save these values
      iteration <- ifelse(i < 3,ifelse(i== 1,i*j,i+j),(i*2)+j)
      ehc_boot[iteration,1] <- i
      ehc_boot[iteration,2] <- string
      ehc_boot[iteration,3] <- sum(prob)
      ehc_boot[iteration,3] <- sum(prob)
      
      ### the following code is to force all 1-occasion sites to be 0..
      ### which is not necessary - keeping code in for reference
      # if (string == "0..") { # if 0.., expected will be equal to number of sites in cohort
      #   ehc_boot[iteration,3] <- nrow(sc)
      # } else if (string == "1..") { # if 1.., expected will need to be 0
      #   ehc_boot[iteration,3] <- 0
      # } else { # for all other capt hists, estimate accordingly
      #   # calculate the expected number of sites with that capt hist
      #   ehc_boot[iteration,3] <- sum(prob)
      # }

    }
  }
  
  # Adjustments:
  ohc_boot <- as.data.frame(ohc_boot)
  ohc_boot$V3 <- as.numeric(ohc_boot$V3)
  ehc_boot <- as.data.frame(ehc_boot)
  ehc_boot$V3 <- as.numeric(ehc_boot$V3)
  # update ohc into only 0s and 1s
  ohc0 <- sum(ohc_boot$V3[which(startsWith(ohc_boot$V2,"0"))])
  ohc1 <- sum(ohc_boot$V3[which(startsWith(ohc_boot$V2,"1"))])
  # update ehc into only 0s and 1s
  ehc0 <- sum(ehc_boot$V3[which(startsWith(ehc_boot$V2,"0"))])
  ehc1 <- sum(ehc_boot$V3[which(startsWith(ehc_boot$V2,"1"))])
  
  # 4. Calculate chi-square test statistic JUST ON SINGLE-OCCASION SITES
  ohc_vec <- c(ohc0,ohc1)
  ehc_vec <- c(ehc0,ehc1)
  test <- ((ohc_vec-ehc_vec)^2)/ehc_vec
  test_stat_boot <- sum(test)
  
  # 5. Save to list
  chi_sq_boot04_23[iter] <- test_stat_boot
}

saveRDS(chi_sq_boot04_23, "output/chi_sq_boot04_23.RData")


# get mean and sd of chi_sq_boot
chi_mn <- mean(chi_sq_boot04_23)
chi_sd <- sd(chi_sq_boot04_23)
# plot a truncated normal dist over it to see if it matches
x <- seq(0,20,0.1)
hist(chi_sq_boot04_23,freq=F,breaks=50)
abline(v = chi_sq_obs04_23, col = "red")

length(which(chi_sq_boot04_23 > chi_sq_obs04_23))/length(chi_sq_boot04_23)
# p-value of 0.00 for 1000




###############################################################################
# old code - using unmarked package ---------------------------------------
# Using function built in to unmarked package
# will not work for my study design

library(unmarked)
library(AICcmodavg)

# Load 2023 capture history
ch_sh <- readRDS("ch_sh.RData")
# need to reassign transect 287 from .0. to 0..
ch_sh[which(ch_sh$Transect == 287),c(2,3)] <- c("0",".")
ch_sh$ch <- paste(ch_sh$V1,ch_sh$V2,ch_sh$V3,sep="")

# Get covariate values for each transect
stnd <- read.csv("Newstandardization.csv",row.names = NULL)
# standardize
for (i in c(17,25,28)) {
  var <- colnames(ch_sh)[i]
  eval(parse(text=paste0("ch_sh$",var,"Stnd <- (ch_sh$",var," - mean(ch_sh$",
                         var,"))/sd(ch_sh$",var,")")))
}

# Build unmarked capt hist
y <- ch_sh[,2:4]
siteCovs <- ch_sh[,c(17,25,28)]
obsCovs <- list(sc = ch_sh[,5:7])

# Replace "." with NA
y[y == "."] <- NA
y$V1 <- as.numeric(y$V1)
y$V2 <- as.numeric(y$V2)
y$V3 <- as.numeric(y$V3)
obsCovs$sc[obsCovs$sc == "9"] <- NA

capthist <- unmarkedFrameOccu(y, siteCovs, obsCovs = NULL)
siteCovs(capthist) <- scale(siteCovs(capthist))
summary(capthist)


occModel <- occu(~1 ~1 + sf2 + PercF4500 + PercESpF4500, capthist)

chisq <- function(fm) {
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

pb <- parboot(occModel, statistic=chisq, nsim=100, parallel=FALSE)

# or 

##compute observed chi-square
obs <- mb.chisq(occModel)
##compute observed chi-square, assess significance, and estimate c-hat
obs.boot <- mb.gof.test(occModel, nsim = 1000)
obs.boot

