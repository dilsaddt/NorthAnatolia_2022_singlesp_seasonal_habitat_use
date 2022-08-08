###########################################################################
# 3. Model diagnostics
###########################################################################

# load(path/to/model/output/model_output.RData)
# Here our model output name is cc_072022_full

## Convergence and distribution check ----
###############################################

# checking Rhat values, if there ara any over 1.1
hist(cc_072022_full$summary[,8])
length(which(cc_072022_full$summary[,8]>1.1))

# visual traceplot checking
par(mfrow=c(3,3))
traceplot(cc_072022_full, parameters = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"))

par(mfrow = c(1,1))
whiskerplot(cc_072022_full, parameters = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"))
whiskerplot(cc_072022_full, parameters = c("n.occ"))

MCMCplot(cc_072022_full$samples, params = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"))
MCMCplot(cc_072022_full$samples, params = c("randompsi", "randomgamma", "randomeps", "randomp"))

MCMCtrace(cc_072022_full$samples, params = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"), pdf=T, ind = T, Rhat = T, n.eff = T)
MCMCtrace(cc_072022_full$samples, params = c("randompsi", "randomgamma", "randomeps", "randomp"), pdf=F, ind = T, Rhat = T, n.eff = T)


# Parameter distributions - caterpillar plots in ms -----------------------

MCMCplot(cc_072022_full$samples, 
         params = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"), horiz=F,
         labels = c(expression(paste(psi,"1_", alpha, "_","BL")),
                    expression(paste(psi,"1_", alpha, "_","CF")),
                    expression(paste(psi,"1_", alpha, "_","MF")),
                    expression(paste(psi,"1_", alpha, "_","O")),
                    expression(paste(psi,"1_", beta, "_","BL")),
                    expression(paste(psi,"1_", beta, "_","CF")),
                    expression(paste(psi,"1_", beta, "_","MF")),
                    expression(paste(psi,"1_", beta, "_","O")),
                    expression(paste(gamma, "_", alpha)),
                    expression(paste(gamma, "_", beta, "_", "Winter")),
                    expression(paste(gamma, "_", beta, "_", "PopDen")),
                    expression(paste(gamma, "_", beta, "_", "Elevation")),
                    expression(paste(epsilon, "_", alpha)),
                    expression(paste(epsilon, "_", beta, "_", "Winter")),
                    expression(paste(epsilon, "_", beta, "_", "PopDen")),
                    expression(paste(epsilon, "_", beta, "_", "Elevation")),
                    expression(paste(p,"_", alpha, "_","BL")),
                    expression(paste(p,"_", alpha, "_","CF")),
                    expression(paste(p,"_", alpha, "_","MF")),
                    expression(paste(p,"_", alpha, "_","O")),
                    expression(paste(p,"_", beta, "_","Winter_BL")),
                    expression(paste(p,"_", beta, "_","Winter_CF")),
                    expression(paste(p,"_", beta, "_","Winter_MF")),
                    expression(paste(p,"_", beta, "_","Winter_O"))),
         guide_lines = T,
         offset = 0.1,
         sz_labels = 4, 
         sz_ax_txt = 5,
         sz_tick_txt = 4,
         sz_med =6, 
         sz_thick = 20,
         sz_thin = 8,
         mar = c(6, 8, 2, 2))


# random effects

MCMCplot(cc_072022_full$samples, 
         params = c("randompsi"), horiz=F,
         labels = c(expression(paste("RE_", psi, "[1]")),
                    expression(paste("RE_", psi, "[2]")),
                    expression(paste("RE_", psi, "[3]")),
                    expression(paste("RE_", psi, "[4]")),
                    expression(paste("RE_", psi, "[5]")),
                    expression(paste("RE_", psi, "[6]")),
                    expression(paste("RE_", psi, "[7]")),
                    expression(paste("RE_", psi, "[8]")),
                    expression(paste("RE_", psi, "[9]")),
                    expression(paste("RE_", psi, "[10]"))),
         guide_lines = T,
         offset = 0.1,
         sz_labels = 4, 
         sz_ax_txt = 5,
         sz_tick_txt = 4,
         sz_med =6, 
         sz_thick = 20,
         sz_thin = 8,
         mar = c(6, 8, 2, 2))

MCMCplot(cc_072022_full$samples, 
         params = c( "randomgamma"), horiz=F,
         labels = c(expression(paste("RE_", gamma, "[1]")),
                    expression(paste("RE_", gamma, "[2]")),
                    expression(paste("RE_", gamma, "[3]")),
                    expression(paste("RE_", gamma, "[4]")),
                    expression(paste("RE_", gamma, "[5]")),
                    expression(paste("RE_", gamma, "[6]")),
                    expression(paste("RE_", gamma, "[7]")),
                    expression(paste("RE_", gamma, "[8]")),
                    expression(paste("RE_", gamma, "[9]")),
                    expression(paste("RE_", gamma, "[10]"))),
         guide_lines = T,
         offset = 0.1,
         sz_labels = 4, 
         sz_ax_txt = 5,
         sz_tick_txt = 4,
         sz_med =6, 
         sz_thick = 20,
         sz_thin = 8,
         mar = c(6, 8, 2, 2))

MCMCplot(cc_072022_full$samples, 
         params = c("randomeps"), horiz=F,
         labels = c(expression(paste("RE_", epsilon, "[1]")),
                    expression(paste("RE_", epsilon, "[2]")),
                    expression(paste("RE_", epsilon, "[3]")),
                    expression(paste("RE_", epsilon, "[4]")),
                    expression(paste("RE_", epsilon, "[5]")),
                    expression(paste("RE_", epsilon, "[6]")),
                    expression(paste("RE_", epsilon, "[7]")),
                    expression(paste("RE_", epsilon, "[8]")),
                    expression(paste("RE_", epsilon, "[9]")),
                    expression(paste("RE_", epsilon, "[10]"))),
         guide_lines = T,
         offset = 0.1,
         sz_labels = 4, 
         sz_ax_txt = 5,
         sz_tick_txt = 4,
         sz_med =6, 
         sz_thick = 20,
         sz_thin = 8,
         mar = c(6, 8, 2, 2))

MCMCplot(cc_072022_full$samples, 
         params = c("randomp"), horiz=F,
         labels = c(expression(paste("RE_", p, "[1]")),
                    expression(paste("RE_", p, "[2]")),
                    expression(paste("RE_", p, "[3]")),
                    expression(paste("RE_", p, "[4]")),
                    expression(paste("RE_", p, "[5]")),
                    expression(paste("RE_", p, "[6]")),
                    expression(paste("RE_", p, "[7]")),
                    expression(paste("RE_", p, "[8]")),
                    expression(paste("RE_", p, "[9]")),
                    expression(paste("RE_", p, "[10]"))),
         guide_lines = T,
         offset = 0.1,
         sz_labels = 4, 
         sz_ax_txt = 5,
         sz_tick_txt = 4,
         sz_med =6, 
         sz_thick = 20,
         sz_thin = 8,
         mar = c(6, 8, 2, 2))

