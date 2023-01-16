###########################################################################
# 3. Model diagnostics
###########################################################################

# load(path/to/model/output/model_output.RData)
# Here our model output name is cc_112022_additive
# This is an example script, shown with the output of additive model
# Same code can be run with interaction models, but, one should be careful about the parameter number and the names then

## Convergence and distribution check ----
###############################################

# checking Rhat values, if there ara any over 1.1
hist(cc_112022_additive$summary[,8])
length(which(cc_112022_additive$summary[,8]>1.1))

# visual traceplot checking
ni=1000000
MCMCtrace(cc_112022_additive$samples, params = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"), pdf=F, ind = T, Rhat = T, n.eff = T, iter = ni)
MCMCtrace(cc_112022_additive$samples, params = c("randompsi", "randomgamma", "randomeps", "randomp"), pdf=F, ind = T, Rhat = T, n.eff = T, iter = ni)


# Parameter distributions -----------------------

MCMCplot(cc_112022_additive$samples, 
         params = c("alphapsi", "betapsi", "alphagamma", "betagamma", "alphaeps", "betaeps", "alphap", "betap"), 
         labels = c(expression(paste(psi,"1_", alpha, "_","BL")),
                    expression(paste(psi,"1_", alpha, "_","CF")),
                    expression(paste(psi,"1_", alpha, "_","MF")),
                    expression(paste(psi,"1_", alpha, "_","HLU")),
                    expression(paste(psi,"1_", beta, "_","BL")),
                    expression(paste(psi,"1_", beta, "_","CF")),
                    expression(paste(psi,"1_", beta, "_","MF")),
                    expression(paste(psi,"1_", beta, "_","HLU")),
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
                    expression(paste(p,"_", alpha, "_","HLU")),
                    expression(paste(p,"_", beta, "_","Winter_BL")),
                    expression(paste(p,"_", beta, "_","Winter_CF")),
                    expression(paste(p,"_", beta, "_","Winter_MF")),
                    expression(paste(p,"_", beta, "_","Winter_HLU"))),
         horiz=T,
         ref_ovl = T,
         guide_lines = T, 
         xlim = c(-6,6),
         xlab = NULL,
         offset = 0.1,
         sz_labels = 4, 
         #sz_ax_txt = 5,
         sz_tick_txt = 3,
         sz_med = 6, 
         sz_thick = 20,
         sz_thin = 8,
         mar = c(6, 8, 2, 2))



# random effects

MCMCplot(cc_112022_additive$samples, 
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

MCMCplot(cc_112022_additive$samples, 
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

MCMCplot(cc_112022_additive$samples, 
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

MCMCplot(cc_112022_additive$samples, 
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
