# NorthAnatolia_2022_singlesp_seasonal_habitat_use

The code for analysing seasonal habitat use of large mammals in northwestern Anatolia by using camera-trapping data 
and single-species dynamic occupancy models.

In total eight species were analyzed: 
brown bear (Ursus arctos), Eurasian lynx (Lynx lynx), gray wolf (Canis lupus), red fox (Vulpes vulpes), 
wild boar (Sus scrofa), roe deer (Capreolus capreolus), European hare (Lepus europaeus), and red deer (Cervus elaphus).

Here, as an example only the analysis of roe deer (Capreolus capreolus) is given. 
For more information and any questions please do not hesitate to get in touch: dilsad.dagtekin@ieu.uzh.ch

This is an updated analysis based on reviewer comments on July 2022.


## Files

cc_dh_052021.csv = detection history of the Capreolus capreolus for JAGS analysis, rows are sites, columns are months (occasions)

cc_dat_072022.csv = vectorized detection history of the Capreolus capreolus with site- and time covariates for JAGS analysis

stations.csv = list of station names to number them

stations_with_area.csv = same with stations.csv including the area information, used as random effect in the models

variable_data.csv = data frame site-specific environmental and anthropogenic covariates

cc_072022_full.txt = JAGS model text file

cc_072022_full_gof.txt = JAGS model text file for goodness-of-fit test
