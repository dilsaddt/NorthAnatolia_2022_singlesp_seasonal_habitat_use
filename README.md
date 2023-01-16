# NorthAnatolia_2022_singlesp_seasonal_habitat_use

The code for analysing seasonal habitat use of large mammals in northwestern Anatolia by using camera-trapping data 
and single-species dynamic occupancy models.

In total eight species were analyzed: 
brown bear (Ursus arctos), Eurasian lynx (Lynx lynx), gray wolf (Canis lupus), red fox (Vulpes vulpes), 
wild boar (Sus scrofa), roe deer (Capreolus capreolus), European hare (Lepus europaeus), and red deer (Cervus elaphus).

Here, as an example only the analysis of roe deer (Capreolus capreolus) is given. 
For more information and any questions please do not hesitate to get in touch: dilsad.dagtekin@ieu.uzh.ch

This is an updated analysis based on reviewer comments on November 2022.


## Files

### Data folder

cc_dh_112022.csv = detection history of the Capreolus capreolus for JAGS analysis, rows are sites, columns are months (occasions)

cc_dat_112022.csv = vectorized detection history of the Capreolus capreolus with site- and time covariates for JAGS analysis

stations.csv = list of station names to number them

stations_with_area.csv = same with stations.csv including the area information, used as random effect in the models

camtrap_effort_Ndays.csv = number of days a camera trap was working at a given site/month

variable_data.csv = data frame site-specific environmental and anthropogenic covariates

### Additive folder

cc_11022_additive.txt = JAGS model text file for the additive model

cc_11022_additive_gof.txt = JAGS model text file for additive model's goodness-of-fit test
