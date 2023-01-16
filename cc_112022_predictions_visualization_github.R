###########################################################################
# Predictions and visualization
###########################################################################

# load(path/to/model/output/model_output.RData)
# Here our model output name is cc_112022_additive

###########################################################################

## Loading libraries ----
##############################
library(ggplot)

setwd("path_to_data_files/")

## 1.3. Loading data ----
#########################

load("path_to_model_output")

## 1.3.1. Variable data ----
###########################

variable_data <- read.csv("variable_data.csv")[,-1]
head(variable_data)

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

popden_mean <- mean(variable_data$PopDen)
popden_sd <- sd(variable_data$PopDen)
popden_std <- (variable_data$PopDen - popden_mean)/popden_sd

###########################################################################
# Predictions
###########################################################################

## 1. Colonization probability ----
###########################################################################

## 1.1. Colonization ~ season + popden + elevation + re(area) ----
###########################################################################

popden_d <- seq(min(popden), max(popden), length.out = 10)
elevation_d <- seq(min(elevation), max(elevation), length.out = 10)

# checking popden effect
# to check the popden effect we need to fix elevation to its mean value, which is 0 here because we standardized it.

#length.out popden X  seasons X mcmc list
# we don't need to plot different areas, it was just for accounting in the model
col_popden <- array(NA, dim=c(10, 2, cc_112022_additive$mcmc.info$n.samples))

for(p in 1:10){ # for popden
  for(i in 1:10){ # for area
    
    col_popden[p,1,] <- plogis(cc_112022_additive$sims.list$alphagamma
                               + cc_112022_additive$sims.list$betagamma[,2]*popden_d[p]
                               + cc_112022_additive$sims.list$betagamma[,3]*0 # mean elev value = 0
                               + cc_112022_additive$sims.list$randomgamma[,i]) # summer = 0 so betagamma1 will be gone
    
    col_popden[p,2,] <- plogis(cc_112022_additive$sims.list$alphagamma
                               + cc_112022_additive$sims.list$betagamma[,1]*1 # winter = 1 
                               + cc_112022_additive$sims.list$betagamma[,2]*popden_d[p]
                               + cc_112022_additive$sims.list$betagamma[,3]*0 # mean elev value = 0
                               + cc_112022_additive$sims.list$randomgamma[,i]) 
  }
}

str(col_popden)

# then we take the mean of the mcmc list
pm.col_popden <- apply(col_popden, c(1,2), mean)
str(pm.col_popden)

# then calculate the credible intervals
CRI.col_popden <- apply(col_popden, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.col_popden)

col_popden.prob <- expand.grid(popden = popden_d, season=c("Summer to Winter", "Winter to Summer"))

col_popden.prob$pred[1:10] <- pm.col_popden[,1]
col_popden.prob$pred[11:20] <- pm.col_popden[,2]

col_popden.prob$pred <- round(col_popden.prob$pred,3)

col_popden.prob$lower[1:10] <- CRI.col_popden[1,1:10,1]
col_popden.prob$lower[11:20] <- CRI.col_popden[1,1:10,2]

col_popden.prob$upper[1:10] <- CRI.col_popden[2,1:10,1]
col_popden.prob$upper[11:20] <- CRI.col_popden[2,1:10,2]

col_popden.prob$species <- "Roe Deer"

head(col_popden.prob)
str(col_popden.prob)

ggplot(col_popden.prob)+
  geom_ribbon(aes(x= popden, ymin = lower, ymax = upper, fill=season, group=season),alpha=0.3)+
  geom_point(aes(popden, pred, col=season, group=season ), size=5) +
  geom_line(aes(popden, pred, col=season, group=season), lwd=2) +
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Colonization Probabilitiy", " ", (gamma)))) +
  xlab("Rural human population density") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "top")


# checking elevation effect
# to check the elevation effect we need to fix popden to its mean value, which is 0 here because we standardized it.

#length.out elevation X  seasons X mcmc list
# we don't need to plot different areas, it was just for accounting in the model
col_elev <- array(NA, dim=c(10, 2, cc_112022_additive$mcmc.info$n.samples))

for(e in 1:10){ # for elevation
  for(i in 1:10){ # for area
    
    col_elev[e,1,] <- plogis(cc_112022_additive$sims.list$alphagamma
                             + cc_112022_additive$sims.list$betagamma[,2]*0 # mean popden value = 0
                             + cc_112022_additive$sims.list$betagamma[,3]*elevation_d[e]
                             + cc_112022_additive$sims.list$randomgamma[,i]) # summer = 0 so betagamma1 will be gone
    
    col_elev[e,2,] <- plogis(cc_112022_additive$sims.list$alphagamma
                             + cc_112022_additive$sims.list$betagamma[,1]*1 # winter = 1 
                             + cc_112022_additive$sims.list$betagamma[,2]*0 # mean popden value = 0
                             + cc_112022_additive$sims.list$betagamma[,3]*elevation_d[e]
                             + cc_112022_additive$sims.list$randomgamma[,i]) 
  }
}

str(col_elev)

# then we take the mean of the mcmc list
pm.col_elev <- apply(col_elev, c(1,2), mean)
str(pm.col_elev)

# then calculate the credible intervals
CRI.col_elev <- apply(col_elev, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.col_elev)

col_elev.prob <- expand.grid(elevation = elevation_d, season=c("Summer to Winter", "Winter to Summer"))

col_elev.prob$pred[1:10] <- pm.col_elev[,1]
col_elev.prob$pred[11:20] <- pm.col_elev[,2]

col_elev.prob$pred <- round(col_elev.prob$pred,3)

col_elev.prob$lower[1:10] <- CRI.col_elev[1,1:10,1]
col_elev.prob$lower[11:20] <- CRI.col_elev[1,1:10,2]

col_elev.prob$upper[1:10] <- CRI.col_elev[2,1:10,1]
col_elev.prob$upper[11:20] <- CRI.col_elev[2,1:10,2]

col_elev.prob$species <- "Roe Deer"

head(col_elev.prob)
str(col_elev.prob)

ggplot(col_elev.prob)+
  geom_ribbon(aes(x= elevation, ymin = lower, ymax = upper, fill=season, group=season),alpha=0.3)+
  geom_point(aes(elevation, pred, col=season, group=season ), size=5) +
  geom_line(aes(elevation, pred, col=season, group=season), lwd=2) +
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Colonization Probabilitiy", " ", (gamma)))) +
  xlab("Elevation (m)") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "top")


## 1.2. Desertion season + popden + elevation + re(area) ----
###########################################################################

popden_d <- seq(min(popden), max(popden), length.out = 10)
elevation_d <- seq(min(elevation), max(elevation), length.out = 10)

# checking popden effect
# to check the popden effect we need to fix elevation to its mean value, which is 0 here because we standardized it.

#length.out popden X  seasons X mcmc list
# we don't need to plot different areas, it was just for accounting in the model
ext_popden <- array(NA, dim=c(10, 2, cc_112022_additive$mcmc.info$n.samples))

for(p in 1:10){ # for popden
  for(i in 1:10){ # for area
    
    ext_popden[p,1,] <- plogis(cc_112022_additive$sims.list$alphaeps
                               + cc_112022_additive$sims.list$betaeps[,2]*popden_d[p]
                               + cc_112022_additive$sims.list$betaeps[,3]*0 # mean elev value = 0
                               + cc_112022_additive$sims.list$randomeps[,i]) # summer = 0 so betagamma1 will be gone
    
    ext_popden[p,2,] <- plogis(cc_112022_additive$sims.list$alphaeps
                               + cc_112022_additive$sims.list$betaeps[,1]*1 # winter = 1 
                               + cc_112022_additive$sims.list$betaeps[,2]*popden_d[p]
                               + cc_112022_additive$sims.list$betaeps[,3]*0 # mean elev value = 0
                               + cc_112022_additive$sims.list$randomeps[,i]) 
  }
}

str(ext_popden)

# then we take the mean of the mcmc list
pm.ext_popden <- apply(ext_popden, c(1,2), mean)
str(pm.ext_popden)

# then calculate the credible intervals
CRI.ext_popden <- apply(ext_popden, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.ext_popden)

ext_popden.prob <- expand.grid(popden = popden_d, season=c("Summer to Winter", "Winter to Summer"))

ext_popden.prob$pred[1:10] <- pm.ext_popden[,1]
ext_popden.prob$pred[11:20] <- pm.ext_popden[,2]

ext_popden.prob$pred <- round(ext_popden.prob$pred,3)

ext_popden.prob$lower[1:10] <- CRI.ext_popden[1,1:10,1]
ext_popden.prob$lower[11:20] <- CRI.ext_popden[1,1:10,2]

ext_popden.prob$upper[1:10] <- CRI.ext_popden[2,1:10,1]
ext_popden.prob$upper[11:20] <- CRI.ext_popden[2,1:10,2]

ext_popden.prob$species <- "Roe Deer"

head(ext_popden.prob)
str(ext_popden.prob)

ggplot(ext_popden.prob)+
  #geom_errorbar(aes(areaNr,ymin = lower, ymax = upper, col=season, group=season), size = 1.5, width = 0.2) +
  geom_ribbon(aes(x= popden, ymin = lower, ymax = upper, fill=season, group=season),alpha=0.3)+
  geom_point(aes(popden, pred, col=season, group=season ), size=5) +
  geom_line(aes(popden, pred, col=season, group=season), lwd=2) +
  #scale_color_hue(direction = -1) +
  # geom_text(data=col.prob[1:2,],aes(1:2,pred,label=pred),hjust=0.5, vjust=-3, size=6, fontface="bold")+
  #add_phylopic(roedeer, x=1.5, y=0.50, ysize=1.5,alpha = 0.3, color = "black") +
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Desertion Probabilitiy", " ", (epsilon)))) +
  #xlab(expression("Rural human population density "(per~'4'~km^2))) +
  xlab("Rural human population density") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "top")


# checking elevation effect
# to check the elevation effect we need to fix popden to its mean value, which is 0 here because we standardized it.

#length.out elevation X  seasons X mcmc list
# we don't need to plot different areas, it was just for accounting in the model
ext_elev <- array(NA, dim=c(10, 2, cc_112022_additive$mcmc.info$n.samples))

for(e in 1:10){ # for elevation
  for(i in 1:10){ # for area
    
    ext_elev[e,1,] <- plogis(cc_112022_additive$sims.list$alphaeps
                             + cc_112022_additive$sims.list$betaeps[,2]*0 # mean popden value = 0
                             + cc_112022_additive$sims.list$betaeps[,3]*elevation_d[e]
                             + cc_112022_additive$sims.list$randomeps[,i]) # summer = 0 so betagamma1 will be gone
    
    ext_elev[e,2,] <- plogis(cc_112022_additive$sims.list$alphaeps
                             + cc_112022_additive$sims.list$betaeps[,1]*1 # winter = 1 
                             + cc_112022_additive$sims.list$betaeps[,2]*0 # mean popden value = 0
                             + cc_112022_additive$sims.list$betaeps[,3]*elevation_d[e]
                             + cc_112022_additive$sims.list$randomeps[,i]) 
  }
}

str(ext_elev)

# then we take the mean of the mcmc list
pm.ext_elev <- apply(ext_elev, c(1,2), mean)
str(pm.ext_elev)

# then calculate the credible intervals
CRI.ext_elev <- apply(ext_elev, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.ext_elev)

ext_elev.prob <- expand.grid(elevation = elevation_d, season=c("Summer to Winter", "Winter to Summer"))

ext_elev.prob$pred[1:10] <- pm.ext_elev[,1]
ext_elev.prob$pred[11:20] <- pm.ext_elev[,2]

ext_elev.prob$pred <- round(ext_elev.prob$pred,3)

ext_elev.prob$lower[1:10] <- CRI.ext_elev[1,1:10,1]
ext_elev.prob$lower[11:20] <- CRI.ext_elev[1,1:10,2]

ext_elev.prob$upper[1:10] <- CRI.ext_elev[2,1:10,1]
ext_elev.prob$upper[11:20] <- CRI.ext_elev[2,1:10,2]

ext_elev.prob$species <- "Roe Deer"

head(ext_elev.prob)
str(ext_elev.prob)


ggplot(ext_elev.prob)+
  geom_ribbon(aes(x= elevation, ymin = lower, ymax = upper, fill=season, group=season),alpha=0.3)+
  geom_point(aes(elevation, pred, col=season, group=season ), size=5) +
  geom_line(aes(elevation, pred, col=season, group=season), lwd=2) +
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Desertion Probabilitiy", " ", (epsilon)))) +
  xlab("Elevation (m)") +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "top")

## 1.3. Detection ~ season + habitat + season:habitat ----
###########################################################################

# Broad Leaved Forest (BL) = 1
# Coniferous Forest (CF) = 2 
# Mixed Forest (MF) = 3
# Human Land-Use Areas (O) = 4

# Winter = 1
# Summer = 0

# season X habitat X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

det <- array(NA, dim = c(2,4,cc_112022_additive$mcmc.info$n.samples))

for (i in 1:10){ # for area
  
  # summer, summer = 0 so no betap
  det[1,1,] <- plogis(cc_112022_additive$sims.list$alphap[,1] # BL
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  det[1,2,] <- plogis(cc_112022_additive$sims.list$alphap[,2] # CF
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  det[1,3,] <- plogis(cc_112022_additive$sims.list$alphap[,3] # MF
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  det[1,4,] <- plogis(cc_112022_additive$sims.list$alphap[,4] # O
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  # winter, winter = 1
  det[2,1,] <- plogis(cc_112022_additive$sims.list$alphap[,1] # BL
                      + cc_112022_additive$sims.list$betap[,1]*1
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  det[2,2,] <- plogis(cc_112022_additive$sims.list$alphap[,2] # CF
                      + cc_112022_additive$sims.list$betap[,2]*1
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  det[2,3,] <- plogis(cc_112022_additive$sims.list$alphap[,3] # MF
                      + cc_112022_additive$sims.list$betap[,3]*1
                      +cc_112022_additive$sims.list$randomp[,i]) 
  
  det[2,4,] <- plogis(cc_112022_additive$sims.list$alphap[,4] # O
                      + cc_112022_additive$sims.list$betap[,4]*1
                      +cc_112022_additive$sims.list$randomp[,i]) 
}

str(det)  

# then we take the mean of the mcmc list
pm.det <- apply(det, c(1,2), mean)
str(pm.det)

# then calculate the credible intervals
CRI.det <- apply(det, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.det)

det.prob <- expand.grid(habitat=c("Broad Leaved Forest", "Coniferous Forest", "Mixed Forest", "Human Land-Use Areas"), season=c("Summer to Winter", "Winter to Summer"))

det.prob$pred[1:4] <- pm.det[1,]
det.prob$pred[5:8] <- pm.det[2,]

det.prob$lower[1:4] <- CRI.det[1,1,1:4]
det.prob$lower[5:8] <- CRI.det[1,2,1:4]

det.prob$upper[1:4] <- CRI.det[2,1,1:4]
det.prob$upper[5:8] <- CRI.det[2,2,1:4]

head(det.prob)

ggplot(det.prob)+
  geom_errorbar(aes(season,ymin = lower, ymax = upper, col=habitat, group=habitat), size = 1, width = 0.5, position = position_dodge(width = 0.8))  +
  scale_fill_continuous(type = "viridis") +
  geom_point(aes(season, pred, col=habitat, group=habitat ), size=5, position = position_dodge(width = 0.8))  +
  ylim(c(0,1)) + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Detection Probabilitiy", " ", (p)))) +
  xlab('Season Transition') +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "top")

## 1.4. Proportion of occupied sites ----
#########################################

plot(cc_112022_additive$mean$n.prop)

prop_occu <- data.frame(season=1:22)
prop_occu$pred <- cc_112022_additive$mean$n.prop
prop_occu$lower <- cc_112022_additive$q2.5$n.prop
prop_occu$upper <- cc_112022_additive$q97.5$n.prop
prop_occu$species <- "Roe Deer"
head(prop_occu)

ggplot(prop_occu) + 
  geom_ribbon(aes(x=1:22, ymin=lower, ymax=upper), fill="#cc4c02") +
  geom_point(aes(1:22, pred)) + 
  geom_line(aes(1:22, pred))+ 
  scale_x_continuous(breaks = 1:22) + theme_minimal() + 
  ylim(0,1) +
  theme(panel.grid.minor = element_blank())+
  ylab("Proportion of occupied sites") +
  xlab('Seasons') +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "none")
