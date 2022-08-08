############################################################################
# This script is for the re-analysis of large-mammals in North Anatolia.
# Capreolus capreolus: single-species dynamic occupancy models.
# Data is from 11/2007 to 11/2018. Data prep file: data_prep_052021.R

# We decided to do seasons as:
# Active/Summer season: May-June-July-Aug-Sept-Oct
# Inactive/Winter season: Dec-Jan-Feb-March
# we dropped April and November data to allow for transition periods
# so data needs to be equal from 2007 December to 2018 October: 24096 - 24226
# This makes 11 Winter and 11 Summer. From 2008 Winter until 2018 Summer.

# GLOBAL MODEL AFTER J OF MAMMALOGY REVISION
# initial occupancy (psi1) = habitat + popden + habitat:popden + re(area)
# colonization (gamma) = season + popden + elevation + re(area)
# desertion (epsilon) = season + popden + elevation + re(area)
# detection (p) = season + habitat + season:habitat + re(area)

# added areas as random effect
# standardized all continous covariates by subtracting the mean and dividing by sd

# Date: 12/07/2022
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data
###########################################################################

## 1.1. House keeping ----
##########################

rm(list = ls())

## 1.2. Loading libraries ----
##############################

load.libraries <- function(){
  library(dplyr)
  library(jagsUI)
  library(ggplot2)
  library(ggthemes)
  library(MCMCvis)
  library(rphylopic)
  library(wesanderson)
}

load.libraries()

## 1.3. Loading data ----
#########################

setwd("/Users/dilsaddagtekin/Dropbox/PhD/North_Anatolia/NA_singlesp_2022")

variable_data <- read.csv("variable_data.csv")[,-1]
head(variable_data)

stations_raw <- read.csv("stations.csv")
colnames(stations_raw) <- c('stationNr', 'stationName')
# create stations
stations <- data.frame(stationNr=stations_raw$stationNr[-c(165, 169)],
                       stationName=stations_raw$stationName[-c(165, 169)])
# removed TUL08 and TUL12 stations
length(stations$stationNr) # 171

# add information on 10 sites covering stations
# TAC = Kure Mt. NP
# TAR = Soku WDA
# TBA = Bakacak Mt.
# TDA = Daday
# TEL = Elekdag WDA
# TIL = Ilgaz WDA
# TKA = Kartdag WDA
# TKG = Kurtgirmez
# TTO = Gavurdagi WDA
# TUL = Uluyayla

stations$stationArea <- c(rep("Kure",9), rep("Soku",50), rep("Bakacak",14), rep("Daday",13), rep("Elekdag",13), 
                          rep("Ilgaz",19), rep("Kartdag",20),rep("Kurtgirmez",7), rep("Gavurdagi",12), rep("Uluyayla",14))

stations$stationAreaNr <- c(rep(1,9), rep(2,50), rep(3,14), rep(4,13), rep(5,13), 
                            rep(6,19), rep(7,20),rep(8,7), rep(9,12), rep(10,14))

# write.csv(stations, "stations_with_area.csv")
# sessions: 12/2007 - 11/2018
aprils <- ((2008:2018)*12) + 4
novembers <- ((2008:2018)*12) + 11

sessions <- 24096:24227
sessions <- sessions[-which(sessions %in% c(aprils,novembers))]
length(sessions) # 110

cc_dh <- read.csv("cc_dh_052021.csv")[,-1]
rownames(cc_dh) <- stations$stationNr
colnames(cc_dh) <- sessions
head(cc_dh)


## 1.4. Arranging variable data ----
####################################

# we do not have any records on TUL08 and TUL12 stations, so I decided to discard them here 
variable_data <- variable_data[-c(165,169),]

# arrange everything that needs to be numeric
str(variable_data)

variable_data$Ruggedness <- as.numeric(variable_data$Ruggedness)
variable_data$Elevation <- as.numeric(variable_data$Elevation)

# We standardized the *elevation*, *slope*, *population density*, *distance to populated places* and *distance to nearest road* by making them as close to 0 as possible. 
# we did this by subtracting the mean and dividing by the standard deviation

elevation_mean <- mean(variable_data$Elevation)
elevation_sd <- sd(variable_data$Elevation)
elevation_std <- (variable_data$Elevation - elevation_mean)/elevation_sd

slope_mean <- mean(variable_data$Slope)
slope_sd <- sd(variable_data$Slope)
slope_std <- (variable_data$Slope - slope_mean)/slope_sd

popden_mean <- mean(variable_data$PopDen)
popden_sd <- sd(variable_data$PopDen)
popden_std <- (variable_data$PopDen - popden_mean)/popden_sd

roadden_mean <- mean(variable_data$RoadDen)
roadden_sd <- sd(variable_data$RoadDen)
roadden_std <- (variable_data$RoadDen - roadden_mean)/roadden_sd

distpop_mean <- mean(variable_data$DistPop)
distpop_sd <- sd(variable_data$DistPop)
distpop_std <- (variable_data$DistPop - distpop_mean)/distpop_sd

distroad_mean <- mean(variable_data$DistRoad)
distroad_sd <- sd(variable_data$DistRoad)
distroad_std <- (variable_data$DistRoad - distroad_mean)/distroad_sd

habitat       <- variable_data$Habitat # category
aspect        <- variable_data$Aspect_Cat2 # category
elevation     <- elevation_std
slope         <- slope_std
rugg          <- variable_data$Ruggedness # category
popden        <- popden_std
roadden       <- roadden_std
distpop       <- distpop_std
distroad      <- distroad_std

# prepare variables for season
# we start with Dec 2007 so winter and end with Oct 2018 so summer
# 11 winters and 11 summers: 22 seasons, 11 years

season <- matrix(c("W", "S"), 171,22, byrow=T)     # all seasons


# season as numeric -- 2 level dummy variable: 1/0
season_n <- matrix(NA, nrow(season), ncol(season))
for (i in 1:nrow(season)){
  for (j in 1:ncol(season)) {
    season_n[i,j] <- ifelse(season[i,j]=="W", 1,0)
  }
}
str(season_n)

## 1.5. Prepare data for JAGS ----
##################################

cc_dh <- as.matrix(cc_dh)
# converting detection history matrix into a vector
cc_y <- c(cc_dh)
str(cc_y) # length 18810

# we need to arrange everything now according to this vector and its lengths

# stations names from 1-171 for all secondary occasions we have 110
station <- rep(stations$stationNr, 110)
station <- as.numeric(as.factor(station))
str(station) # length 18810

# area names from 1-10 for all secondary occasions we have 110
area <- rep(stations$stationAreaNr, 110)
area <- as.numeric(as.factor(area))
str(area) # length 18810

# month names from 1-12 for all stations
month <- sessions
month <- rep(c(12,1:3,5:10),11)
month <- rep(month, each=171)
str(month) # length 18810

# secondary occasions from 1-110 for all stations
socc <- rep(1:110, each=171)
str(socc) # length 18810

# year names for all stations
year  <- sessions
year[1] <- 2007
year[2:101] <- rep(2008:2017, each=10)
year[101:110] <- 2018
year <- rep(year, each=171)
str(year) # length 18810

# seasons as dummy variable 1 for winter, 0 for summer
seasons_w <- rep(1, 4)
seasons_s <- rep(0, 6)
seasons <- rep(c(seasons_w, seasons_s),11)
seasons <- rep(seasons, each=171)
str(seasons) # length 18810

# primary occasions from 1-23 for all stations
pocc <- c(rep(1,4), rep(2,6), rep(3,4), rep(4,6), rep(5,4), rep(6,6), 
          rep(7,4), rep(8,6), rep(9,4), rep(10,6), rep(11,4), rep(12,6),
          rep(13,4), rep(14,6), rep(15,4), rep(16,6), rep(17,4), rep(18,6),
          rep(19,4), rep(20,6), rep(21,4), rep(22,6))
pocc <- rep(pocc, each=171)
str(pocc) # length 18810


# Convert categorical variables into numeric for JAGS to recognize them
str(variable_data) # to see which ones are categorical

# habitat -- more than 2 levels
habitat_n <- as.integer(as.factor(habitat)) ; table(habitat_n)
habitat_1 <- NA; habitat_2 <- NA; habitat_3 <- NA; habitat_4 <- NA

for (i in 1:171){
  habitat_1[i] <- ifelse(habitat_n[i]==1,1,0)
  habitat_2[i] <- ifelse(habitat_n[i]==2,1,0)
  habitat_3[i] <- ifelse(habitat_n[i]==3,1,0)
  habitat_4[i] <- ifelse(habitat_n[i]==4,1,0)
}
habitat_n
habitat_4

# aspect -- 2 level dummy variable: 1/0
aspect_n <- NA
for (i in 1:length(aspect)){
  aspect_n[i] <- ifelse(aspect[i]=="N",1,0)
}
aspect_n

# these also need to be compatible with the detection history
# each site covariate for all occasions: 110
aspect.p <- rep(aspect_n, 110)
str(aspect.p) 
habitat.p <- rep(habitat_n, 110)
str(aspect.p) 
elevation.p <- rep(elevation, 110)
str(elevation.p) 
slope.p <- rep(slope, 110)
str(slope.p) 
rugg.p <- rep(rugg, 110)
str(rugg.p) 
popden.p <- rep(popden, 110)
str(popden.p) 
roadden.p <- rep(roadden, 110)
str(roadden.p) 
distpop.p <- rep(distpop, 110)
str(distpop.p) 
distroad.p <- rep(distroad, 110)
str(distroad.p) 


###########################################################################
# 2. JAGS analysis
###########################################################################

## 2.1. Bundle data ----
##########################

# Let's put all these data together and try to compile in order to avoid NAs
# this dat data frame will be used mainly for observational model (p, detection)

dat <- data.frame(y = cc_y)
dat$site        <- station
dat$site        <- as.numeric(as.factor(dat$site)) # we excluded 2 sites, now we need to renumber them
dat$area        <- area
dat$year        <- year
dat$month       <- month
dat$pocc        <- pocc
dat$socc        <- socc
dat$season      <- seasons
dat$aspect      <- aspect.p
dat$habitat     <- habitat.p
dat$elevation   <- elevation.p
dat$slope       <- slope.p
dat$rugg        <- rugg.p
dat$popden      <- popden.p
dat$roadden     <- roadden.p
dat$distpop     <- distpop.p
dat$distroad    <- distroad.p

head(dat)

# taking NAs out
na_out <- which(is.na(dat$y))
table(na_out)
dat <- dat[-na_out,]
str(dat)

####################################################################
# write.csv(dat, "cc_dat_072022.csv")
####################################################################

## 2.2. Analysis ----
##########################

# dh: detection history
# nsite = nrow(cc_dh)        # 171
# nprimary= (ncol(cc_dh)/10)*2    # 22
# nsecondary = ncol(cc_dh)   # 110 OR nsecondary = length(unique(dat$socc))   # 75 

# Bundle data
str(jags.data <- list(y              = dat$y, 
                      nsite          = nrow(cc_dh), 
                      nprimary       = (ncol(cc_dh)/10)*2,
                      nsecondary     = length(unique(dat$socc)),
                      nobs           = nrow(dat),
                      narea          = length(unique(dat$area)),
                      pocc           = dat$pocc,
                      area           = stations$stationAreaNr,
                      site           = dat$site,
                      area.p         = dat$area,
                      season.p       = dat$season,
                      habitat.p      = dat$habitat, 
                      elevation.p    = dat$elevation, 
                      slope.p        = dat$slope, 
                      aspect.p       = dat$aspect, 
                      rugg.p         = dat$rugg,
                      popden.p       = dat$popden, 
                      roadden.p      = dat$roadden,
                      distpop.p      = dat$distpop, 
                      distroad.p     = dat$distroad,
                      season         = season_n,
                      habitat        = habitat_n,
                      elevation      = elevation,
                      slope          = slope,
                      aspect         = aspect_n,
                      rugg           = rugg,
                      popden         = popden,
                      roadden        = roadden,
                      distpop        = distpop,
                      distroad       = distroad))

# Specify model in BUGS language
cat(file = "cc_072022_full.txt"," 
model {

# Model: cc_072022_full
#
# psi1 ~ habitat + popden + habitat:popden + re.psi1(area)
# gamma ~ season + popden + elevation + re.gammma(area)
# eps ~ season + popden + elevation + re.eps(area)
# p ~ season + habitat + season:habitat + re.p(area)
# -----------------------------------------------------------------------------------------
#
# Parameters:
#
# psi1: initial occupancy probability
# gamma: colonization probability
# eps (epsilon) : desertion probability
# p: detection probability
# -----------------------------------------------------------------------------------------
#
# Others:
#
# narea: number of study areas, 10
# nsite: number of camera-trap stations, 171
# nprimary: total number of primary sampling periods, seasons, 22
# nobs: number of observations
# y: detection history
# z: true latent state variable, z=1 means occupied, z=0 not occupied.
# n.occ: number of occupied sites
# n.prop: proportion of occupied sites to all total sites
#
# area: study areas
# season: seasons, 2-level categorical, winter(1) and summer (0)
# popden: rural human population density, continuous, standardized
# elevation: meter, continuous, standardized
# habitat: habitat type, 4-level categorical:
## Broad Leaved Forest (BL) = 1
## Coniferous Forest (CF) = 2 
## Mixed Forest (MF) = 3
## Human Land-Use Areas (O) = 4
# 
# -----------------------------------------------------------------------------------------
#
#### PRIORS ####

# psi1 ~ habitat + popden + habitat:popden + re.psi1(area)

for (h in 1:4) {
  alphapsi[h] ~ dlogis(0, 1)                                    # intercept for each habitat type(4)
  
  # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
  # we have 4 habitat types --> 4 betapsi values for 1 betapsi category
  # we have 1 betapsi category: popden effect
  
  betapsi[h] ~ dnorm(0, 1)                                      # prior for popden coefficient w/ habitat types
}

# gamma ~ season + popden + elevation + re.gamma(area)

alphagamma ~ dlogis(0, 1)                                       # intercept

for(i in 1:3){                                                  # prior for coefficients 
  betagamma[i] ~ dnorm(0, 1)
}

# eps ~ season + popden + elevation + re.eps(area)

alphaeps ~ dlogis(0, 1)                                         # intercept

for(i in 1:3){                                                  # prior for coefficients 
  betaeps[i] ~ dnorm(0, 1)
}

# p ~ season*habitat

for (h in 1:4) {
  alphap[h] ~ dlogis(0, 1)                                      # intercept for each habitat type(4)

  # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
  # we have 4 habitat types --> 4 betap values for 1 betap category
  # we have 1 betap category: betap[p] season effect
  
  betap[h] ~ dnorm(0, 1)                                        # betap season effect
}

# Random effect priors
  
for (t in 1:(narea)){
randompsi[t] ~ dnorm(0, taupsi)
randomgamma[t] ~ dnorm(0, taugamma)
randomeps[t] ~ dnorm(0, taueps)
randomp[t] ~ dnorm(0, taup)

    
}

sigmapsi ~ dunif(0, 10)                     # Priors for standard deviations
sigmagamma ~ dunif(0, 10)
sigmaeps ~ dunif(0, 10)
sigmap ~ dunif(0, 10) 

taupsi <- pow(sigmapsi, -2)
taugamma <- pow(sigmagamma, -2)
taueps <- pow(sigmaeps, -2)
taup <- pow(sigmap, -2)

sigma2psi <- pow(sigmapsi, 2)                  # Temporal variances
sigma2gamma <- pow(sigmagamma, 2)
sigma2eps <- pow(sigmaeps, 2)
sigma2p <- pow(sigmap, 2)

#### MODELS ####

# Ecological submodel: intial occupancy as derived parameter

for (i in 1:nsite){
  z[i,1] ~ dbern(psi1[i])
  logit(psi1[i]) <- alphapsi[habitat[i]] 
  + betapsi[habitat[i]]*popden[i]
  + randompsi[area[i]]
}

# State transitions: colonization and extinction
# for col and ext

for (i in 1:nsite){
  for (s in 2:nprimary){
    logit(gamma[i,s-1]) <- alphagamma 
    + betagamma[1]*season[i,s-1] + betagamma[2]*popden[i] + betagamma[3]*elevation[i]
    + randomgamma[area[i]]
    
    logit(eps[i,s-1]) <- alphaeps 
    + betaeps[1]*season[i,s-1] + betaeps[2]*popden[i] + betaeps[3]*elevation[i]
    + randomeps[area[i]]
  }
}

# for z: true occupancy

for(i in 1:nsite){
  for (t in 2:nprimary){
    
    z[i,t] ~ dbern((z[i,t-1]*(1-eps[i,t-1])) + ((1-z[i,t-1])*gamma[i,t-1]))
  }
}

# Observation model 

for (i in 1:nobs){
  logit(p[i]) <- alphap[habitat.p[i]]
  + betap[habitat.p[i]]*season.p[site[i]]
  + randomp[area.p[i]]
  
    y[i] ~ dbern(z[site[i],pocc[i]]*p[i])  # y: observed occupancy
}

#### DERIVED PARAMETERS ####
# Compute population and sample occupancy

n.occ[1] <- sum(z[1:nsite,1])  # Number of occupied sites in sample
n.prop[1] <- n.occ[1]/nsite

for (t in 2:nprimary){
  n.occ[t] <- sum(z[1:nsite,t])
  n.prop[t] <- n.occ[t]/nsite
}

}")



# Initial values
inits<-function(){list(z=matrix(1,171,22))}

# Parameters monitored
params <- c("alphapsi", "betapsi", 
            "alphagamma", "betagamma", 
            "alphaeps", "betaeps",
            "alphap", "betap",
            "randompsi", "randomgamma", "randomeps", "randomp",
            "n.occ","n.prop")

# MCMC settings 
na <- 100000
nb <- 100000
nt <- 20
ni <- 2000000
nc <- 3


parallel:::setDefaultClusterOptions(setup_strategy = "sequential") # after R update this is needed for parallel run

# call JAGS
start <- Sys.time()

cc_072022_full <- jags(jags.data, inits, params, "cc_072022_full.txt", n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

end <- Sys.time()

end - start

print(cc_072022_full, 3)

str(cc_072022_full$sims.list)

# save(cc_072022_full, file = "cc_072022_full.RData")