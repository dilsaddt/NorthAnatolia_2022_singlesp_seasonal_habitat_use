# NorthAnatolia_2022_singlesp_seasonal_habitat_use

The code for analysing seasonal habitat use of large mammals in northwestern Anatolia by using camera-trapping data 
and single-species dynamic occupancy models.

In total eight species were analyzed: 
brown bear (_Ursus arctos_), Eurasian lynx (_Lynx lynx_), gray wolf (_Canis lupus_), red fox (_Vulpes vulpes_), 
wild boar (_Sus scrofa_), roe deer (_Capreolus capreolus_), European hare (_Lepus europaeus_), and red deer (_Cervus elaphus_).

Here, as an example only the analysis of roe deer (_Capreolus capreolus_) is given. 
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

### Model folders

cc_11022_model_type.txt = JAGS model text file for the named model type

cc_11022_model_type_gof.txt = JAGS model text file for the named model type's goodness-of-fit test
