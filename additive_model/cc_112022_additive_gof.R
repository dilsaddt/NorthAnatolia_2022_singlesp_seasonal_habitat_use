############################################################################
# This script is for the re-analysis of large-mammals in North Anatolia.
# Capreolus capreolus: single-species dynamic occupancy models.
# Data is from 11/2007 to 11/2018.

# gof code for additive model

# Date: 10/11/2022
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
  library(jagsUI)
  library(MCMCvis)
}

load.libraries()

## 1.3. Loading data ----
#########################

setwd("/path_to_files/")
load("path_to_data/jags_data_cc_112022.RData") # load prepared jags data from cc_112022_data_prep.R

###########################################################################
# 2. JAGS analysis
###########################################################################

# Specify model in BUGS language
cat(file = "cc_112022_additive_gof.txt"," 
model {

# Model: cc_112022_additive_gof
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
  # alphap[h] ~ dlogis(0, 1)                                      # intercept for each habitat type(4)
  
  alphap[h] ~ dlogis(-2, 1)
    
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

sigmapsi ~ dnorm(0, 1)T(0,)                    # Priors for standard deviations
sigmagamma ~ dnorm(0, 1)T(0,)
sigmaeps ~ dnorm(0, 1)T(0,)
sigmap ~ dnorm(0, 1)T(0,)

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
  
   #y[i] ~ dbern(z[site[i],pocc[i]]*p[i])  # y: observed occupancy
  
   pstar[i] <- 1-(1-p[i])^ndays[i] 
   y[i] ~ dbern(z[site[i],pocc[i]]*pstar[i])  # y: observed occupancy
}

#### DERIVED PARAMETERS ####

# derived occupancy state
for(i in 1:nsite){
  
  psi[i, 1] <- psi1[i] # Population occupancy as derived quantity
  for (t in 2:nprimary){
    
    psi[i, t] <- (psi[i,t-1]*(1-eps[i,t-1])) + ((1-psi[i,t-1])*gamma[i,t-1])
  
  }
}

# Compute population and sample occupancy

n.occ[1] <- sum(z[1:nsite,1])  # Number of occupied sites in sample
n.prop[1] <- n.occ[1]/nsite

for (t in 2:nprimary){
  n.occ[t] <- sum(z[1:nsite,t])
  n.prop[t] <- n.occ[t]/nsite
}

#### GOODNESS-OF-FIT GoF ####
# Based on posterior predictive distribution
# -------------------------------------------------------------------------

# Draw a replicate data set under the fitted model

for (i in 1:nobs){
  
    yrep[i] ~ dbern(z[site[i],pocc[i]]*pstar[i])
}

}")



# Initial values
inits<-function(){list(z=matrix(1,171,22))}

# Parameters monitored
params <- c("yrep", "psi", "gamma","eps","z","pstar")

# MCMC settings 
na <- 80000
nb <- 80000
nt <- 20
ni <- 600000
nc <- 3

parallel:::setDefaultClusterOptions(setup_strategy = "sequential") # after R update this is needed for parallel run

# call JAGS
start <- Sys.time()

cc_112022_additive_gof <- jags(jags.data, inits, params, "cc_112022_additive_gof.txt", n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

end <- Sys.time()

end - start


# GOF ---------------------------------------------------------------------

###########################################################################
# FOR ALL MCMC SAMPLES ----------------------------------------------------
###########################################################################

dat_gof <- dat

rownames(dat_gof) <- NULL
str(dat_gof)
dat_gof$site <- as.character(dat_gof$site)
dat_gof$socc <- as.character(dat_gof$socc)
dat_gof$pocc <- as.character(dat_gof$pocc)

dat_gof$soccNr <- ifelse(dat_gof$month == 12, 1,
                         ifelse(dat_gof$month == 1, 2,
                                ifelse(dat_gof$month == 2, 3,
                                       ifelse(dat_gof$month == 3, 4,
                                              ifelse(dat_gof$month == 5, 1,
                                                     ifelse(dat_gof$month == 6, 2,
                                                            ifelse(dat_gof$month == 7, 3,
                                                                   ifelse(dat_gof$month == 8, 4,
                                                                          ifelse(dat_gof$month == 9, 5, 
                                                                                 ifelse(dat_gof$month == 10, 6, NA))))))))))


# create all the output from the model

# taking out the replicated observation data from the model
yrep <- list(NA)
dat_rep <- list(NA)


for (i in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  yrep[[i]] <- cc_112022_additive_gof$sims.list$yrep[i,]
  
  dat_rep[[i]] <- dat_gof[,c(1:9, 19)]
  dat_rep[[i]]$y <- yrep[[i]]
  
}
str(dat_rep[[1]])
str(dat_rep[[2]])

# taking out the estimated detection prob from the model
p_model <- list(NA)
dat_p <- list(NA)

for (i in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  p_model[[i]] <- cc_112022_additive_gof$sims.list$pstar[i,]
  
  dat_p[[i]] <- dat_gof[,c(1:9, 19)]
  dat_p[[i]]$y <- yrep[[i]]
  dat_p[[i]]$pstar <- round(p_model[[i]],2)
  
}

# we need to convert detection histories into 3d array again
# months start with 12/2007

# winter 12/1/2/3
# summer 5/6/7/8/9/10

sites_arr <- 1:171
months_arr <- 1:6
seasons_arr <- 1:11

## observed data: dat_gof$y, only one array
y3dobs <- array(NA, dim = c(171, 6, 22), dimnames = list(1:171, 1:6, 1:22))

for (i in 1:171){
  for (j in 1:6){
    for (t in 1:22){  
      
      y3dobs[i,j,t] <- ifelse(sites_arr[i] %in% dat_gof[which(dat_gof$y==1 & dat_gof$pocc==t & dat_gof$soccNr==j),]$site, 1,
                              ifelse(sites_arr[i] %in% dat_gof[which(dat_gof$y==0 & dat_gof$pocc==t & dat_gof$soccNr==j),]$site, 0, NA))
    }
  }
}


# replicated data: dat_rep, for all MCMC samples
y3drep <- vector(mode = 'list', length = cc_112022_additive_gof$mcmc.info$n.samples)

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  y3drep[[k]] <- array(NaN, dim = c(171, 6, 22), dimnames = list(1:171, 1:6, 1:22))
}

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  for (i in 1:171){
    for (j in 1:6){
      for (t in 1:22){  
        
        y3drep[[k]][i,j,t] <- ifelse(sites_arr[i] %in% dat_rep[[k]][which(dat_rep[[k]]$y==1 & dat_rep[[k]]$pocc==t & dat_rep[[k]]$soccNr==j),]$site, 1,
                                     ifelse(sites_arr[i] %in% dat_rep[[k]][which(dat_rep[[k]]$y==0 & dat_rep[[k]]$pocc==t & dat_rep[[k]]$soccNr==j),]$site, 0, NA))
      }
    }
  }
}

# we also need to convert estimated pstar (detection prob) into 3d array

p_est <- vector(mode = 'list', length = cc_112022_additive_gof$mcmc.info$n.samples)

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  p_est[[k]] <- array(NaN, dim = c(171, 6, 22), dimnames = list(1:171, 1:6, 1:22))
}

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  for (i in 1:171){
    for (j in 1:6){
      for (t in 1:22){  
        
        p_est[[k]][i,j,t] <- ifelse(length(dat_p[[k]][which(dat_p[[k]]$site==i & dat_p[[k]]$pocc ==t & dat_p[[k]]$soccNr==j),]$pstar)!= 0, 
                                    dat_p[[k]][which(dat_p[[k]]$site==i & dat_p[[k]]$pocc ==t & dat_p[[k]]$soccNr  ==j),]$pstar, NA)
      }
    }
  }
}

# (2a) Computations for the GoF of the open part of the model
# (based on number of state transitions)
# ----------------------------------------------------------

nsite = dim(y3drep[[1]])[1] # number of sites
nprimary = dim(y3drep[[1]])[3] # number of seasons (primar occasions)
e = 0.0001 # to avoid dividing by zero

psi <- list(NA)
gamma <- list(NA)
eps <- list(NA)

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  psi[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary)
  gamma[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary-1)
  eps[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary-1)
  
}


for (i in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  psi[[i]] <- cc_112022_additive_gof$sims.list$psi[i,,]
  gamma[[i]] <- cc_112022_additive_gof$sims.list$gamma[i,,]
  eps[[i]] <- cc_112022_additive_gof$sims.list$eps[i,,]
  
}

# Create needed empty matrices

# observed data
zobs <- matrix(NaN, nrow = nsite, ncol=nprimary) # empty matrix to save latent state from observed data

des <- matrix(NA, nrow = nsite, ncol=nprimary-1) 
col <- matrix(NA, nrow = nsite, ncol=nprimary-1)
nondes <- matrix(NA, nrow = nsite, ncol=nprimary-1)
noncol <- matrix(NA, nrow = nsite, ncol=nprimary-1)

tm <- array(NA, dim = c(2,2,nprimary-1))

# replicated data
zobsrep <- list(NA)

desrep <- list(NA)
colrep <- list(NA)
nondesrep <- list(NA)
noncolrep <- list(NA)

tmrep <- list(NA)

noncol.exp <- list(NA)
col.exp <- list(NA)
des.exp <- list(NA)
nondes.exp <- list(NA)

Etm <- list(NA)

x2Open <- list(NA)
x2repOpen <- list(NA)

Chi2Open <- NA
Chi2repOpen <- NA
Chi2ratioOpen <- NA

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  
  zobsrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary)  # empty matrix to save latent state from replicated data
  
  desrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  colrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  nondesrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  noncolrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  
  tmrep[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  
  noncol.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  col.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  des.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  nondes.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  
  Etm[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  
  x2Open[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  x2repOpen[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  
}

# Calculate the Chi-square

# Compute observed z matrix for observed and replicated data
for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  
  for (i in 1:nsite){
    for (t in 1:nprimary){
      zobs[i,t] <- max(y3dobs[i,,t], na.rm=T)  # For observed data
      zobsrep[[k]][i,t] <- max(y3drep[[k]][i,,t], na.rm = T) # For replicated data
    }
    
    # Identify desertions, non-desertions, colonization and non-colonizations
    for (t in 2:nprimary){
      # ... for observed data
      des[i,(t-1)] <- ifelse(zobs[i,t]==0 & zobs[i,t-1]==1,1,0)
      nondes[i,(t-1)] <- ifelse(zobs[i,t]==1 & zobs[i,t-1]==1,1,0)
      col[i,(t-1)] <- ifelse(zobs[i,t]==1 & zobs[i,t-1]==0,1,0)
      noncol[i,(t-1)] <- ifelse(zobs[i,t]==0 & zobs[i,t-1]==0,1,0)
      # ... for replicated data
      desrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==0 & zobsrep[[k]][i,t-1]==1,1,0)
      nondesrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==1 & zobsrep[[k]][i,t-1]==1,1,0)
      colrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==1 & zobsrep[[k]][i,t-1]==0,1,0)
      noncolrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==0 & zobsrep[[k]][i,t-1]==0,1,0)
    }
  }
  
  # Tally up number of transitions and put into a matrix for each year
  for(t in 1:(nprimary-1)){
    # ... for observed data
    tm[1,1,t] <- sum(noncol[,t]) # transition mat for obs. data
    tm[1,2,t] <- sum(col[,t])
    tm[2,1,t] <- sum(des[,t])
    tm[2,2,t] <- sum(nondes[,t])
    # ... for replicated data
    tmrep[[k]][1,1,t] <- sum(noncolrep[[k]][,t]) # transition mat for rep. data
    tmrep[[k]][1,2,t] <- sum(colrep[[k]][,t])
    tmrep[[k]][2,1,t] <- sum(desrep[[k]][,t])
    tmrep[[k]][2,2,t] <- sum(nondesrep[[k]][,t])
  }
  
  # Compute expected numbers of transitions under the model
  # Probability of each individual transition
  for(i in 1:nsite){
    for(t in 1:(nprimary-1)){
      noncol.exp[[k]][i,t] <- (1-psi[[k]][i,t]) * (1-gamma[[k]][i,t])
      col.exp[[k]][i,t] <- (1-psi[[k]][i,t]) * gamma[[k]][i,t]
      des.exp[[k]][i,t] <- psi[[k]][i,t] * eps[[k]][i,t]
      nondes.exp[[k]][i,t] <- psi[[k]][i,t] * (1-eps[[k]][i,t])
    }
  }
  
  
  
  # Sum up over sites to obtain the expected number of those transitions
  for(t in 1:(nprimary-1)){
    Etm[[k]][1,1,t] <- sum(noncol.exp[[k]][,t])
    Etm[[k]][1,2,t] <- sum(col.exp[[k]][,t])
    Etm[[k]][2,1,t] <- sum(des.exp[[k]][,t])
    Etm[[k]][2,2,t] <- sum(nondes.exp[[k]][,t])
  }
  
  
  # Compute Chi-square discrepancy
  for(t in 1:(nprimary-1)){
    # ... for observed data
    x2Open[[k]][1,1,t] <- ((tm[1,1,t] - Etm[[k]][1,1,t]) ^ 2) / (Etm[[k]][1,1,t]+e)
    x2Open[[k]][1,2,t] <- ((tm[1,2,t] - Etm[[k]][1,2,t]) ^ 2) / (Etm[[k]][1,2,t]+e)
    x2Open[[k]][2,1,t] <- ((tm[2,1,t] - Etm[[k]][2,1,t]) ^ 2) / (Etm[[k]][2,1,t]+e)
    x2Open[[k]][2,2,t] <- ((tm[2,2,t] - Etm[[k]][2,2,t]) ^ 2) / (Etm[[k]][2,2,t]+e)
    # ... for replicated data
    x2repOpen[[k]][1,1,t] <- ((tmrep[[k]][1,1,t]-Etm[[k]][1,1,t]) ^ 2)/(Etm[[k]][1,1,t]+e)
    x2repOpen[[k]][1,2,t] <- ((tmrep[[k]][1,2,t]-Etm[[k]][1,2,t]) ^ 2)/(Etm[[k]][1,2,t]+e)
    x2repOpen[[k]][2,1,t] <- ((tmrep[[k]][2,1,t]-Etm[[k]][2,1,t]) ^ 2)/(Etm[[k]][2,1,t]+e)
    x2repOpen[[k]][2,2,t] <- ((tmrep[[k]][2,2,t]-Etm[[k]][2,2,t]) ^ 2)/(Etm[[k]][2,2,t]+e)
  }
  
  # Add up overall test statistic and compute fit stat ratio (open part)
  Chi2Open[[k]] <- sum(x2Open[[k]][,,])       # Chisq. statistic for observed data
  Chi2repOpen[[k]] <- sum(x2repOpen[[k]][,,]) # Chisq. statistic for replicated data
  Chi2ratioOpen[[k]] <- Chi2Open[[k]] / Chi2repOpen[[k]]
}

# plot Open part
pl <- range(c(Chi2Open, Chi2repOpen))
plot(Chi2Open, Chi2repOpen,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(2685, 2755, paste('Bpv = ', round(mean(Chi2repOpen >
                                              Chi2Open), 2)), cex = 2)


# (2b) Computations for the GoF of the closed part of the model
# (based on the number of times detected, i.e., detection freqiencies)
# --------------------------------------------------------------------
nsite = dim(y3drep[[1]])[1] # number of sites
nsecondary = dim(y3drep[[1]])[2] # number of mpnths (secondary occasions)
nprimary = dim(y3drep[[1]])[3] # number of seasons (primary occasions)
e = 0.0001 # to avoid dividing by zero

z <- list(NA)
pstar <- list(NA)

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  z[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary)
  pstar[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary)
  
}


for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  z[[k]] <- cc_112022_additive_gof$sims.list$z[k,,]
  pstar[[k]] <- p_est[[k]]
  
}


# Create needed empty matrices

# observed data
detfreq <- matrix(NA, nrow = nsite, ncol=nprimary) # empty matrix to save latent state from observed data

# replicated data
detfreqrep <- list(NA)

tmp <- list(NA)

E <- list(NA)

x2Closed <- list(NA)
x2repClosed <- list(NA)
ftClosed <- list(NA)
ftrepClosed <- list(NA)

Chi2Closed <- NA
FTClosed <- NA
Chi2repClosed <- NA
FTrepClosed <- NA
Chi2ratioClosed <- NA
FTratioClosed <- NA

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  
  detfreqrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary)  # empty matrix to save latent state from replicated data
  tmp[[k]] <- array(NaN, dim=c(nsite, nsecondary, nprimary))
  E[[k]] <- matrix(NA, nrow = nsite, ncol=nprimary)
  
  x2Closed[[k]] <- matrix(NA, nrow = nsite, ncol=nprimary)
  x2repClosed[[k]] <- matrix(NA, nrow = nsite, ncol=nprimary)
  ftClosed[[k]] <- matrix(NA, nrow = nsite, ncol=nprimary)
  ftrepClosed[[k]] <- matrix(NA, nrow = nsite, ncol=nprimary)
}

# Calculate Chi-square and Freeman-Tuey

for (k in 1:cc_112022_additive_gof$mcmc.info$n.samples){
  
  # Compute detection frequencies for observed and replicated data
  for (i in 1:nsite){
    for (t in 1:nprimary){
      # Det. frequencies for observed and replicated data
      detfreq[i,t] <- sum(y3dobs[i,,t], na.rm = T)
      detfreqrep[[k]][i,t] <- sum(y3drep[[k]][i,,t], na.rm = T)
      # Expected detection frequencies under the model
      for (j in 1:nsecondary){
        tmp[[k]][i,j,t] <- z[[k]][i, t] * pstar[[k]][i,j,t]
      }
      E[[k]][i,t] <- sum(tmp[[k]][i,,t], na.rm = T)     # Expected number of detections
      # Chi-square and Freeman-Tukey discrepancy measures
      # ..... for actual data set
      x2Closed[[k]][i,t] <- ((detfreq[i,t] - E[[k]][i,t]) ^ 2) / (E[[k]][i,t]+e)
      ftClosed[[k]][i,t] <- (sqrt(detfreq[i,t]) - sqrt(E[[k]][i,t])) ^ 2
      # ..... for replicated data set
      x2repClosed[[k]][i,t] <- ((detfreqrep[[k]][i,t] - E[[k]][i,t]) ^ 2) / (E[[k]][i,t]+e)
      ftrepClosed[[k]][i,t] <- (sqrt(detfreqrep[[k]][i,t]) - sqrt(E[[k]][i,t])) ^ 2
    }
  }
  
  # Add up Chi-square and FT discrepancies and compute fit stat ratio (closed part)
  Chi2Closed[[k]] <- sum(x2Closed[[k]][,])
  FTClosed[[k]] <- sum(ftClosed[[k]][,])
  Chi2repClosed[[k]] <- sum(x2repClosed[[k]][,])
  FTrepClosed[[k]] <- sum(ftrepClosed[[k]][,])
  Chi2ratioClosed[[k]] <- Chi2Closed[[k]] / Chi2repClosed[[k]]
  FTratioClosed[[k]] <- FTClosed[[k]] / FTrepClosed[[k]]
}

# Closed part of model: Chi-squared
pl <- range(c(Chi2Closed, Chi2repClosed))
plot(Chi2Closed, Chi2repClosed,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Closed part of model (Chi-squared)", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(230, 390, paste('Bpv = ', round(mean(Chi2repClosed > Chi2Closed), 2)), cex = 2)

# Closed part of model: Freeman-Tukey
pl <- range(c(FTClosed, FTrepClosed))
plot(FTClosed, FTrepClosed,
     xlab = "FT observed data", ylab = "FT expected data",
     main = "Closed part of model (Freeman-Tukey)", xlim = pl, ylim = pl,
     frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(95, 170, paste('Bpv = ', round(mean(FTrepClosed > FTClosed), 2)), cex = 2)
