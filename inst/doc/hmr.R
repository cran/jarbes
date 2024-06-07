## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = FALSE-------------------------------------------------------------
library(kableExtra)

text.matrix = c("Yes", "No", 
                  "AD-Data",  "Empty Cell",
                  "IPD-RWD",  "IPD-RWD"
                )

dim(text.matrix) = c(2, 3)

text.matrix %>% 
  kbl() %>%
  kable_paper() %>%
  kable_classic_2(full_width = F) %>% 
    add_header_above(c("Meet inclusion criteria", "RCTs", "OS")) %>% 
  add_header_above(c("", "Study Type" = 2))  


## ----echo=FALSE, message=FALSE, warning=FALSE,results = "hide"----------------
# Packages
library(ggplot2)
library(gridExtra)
library(mcmcplots)
library(ggplot2)
library(grid)
library(gridExtra)
library(R2jags)
library(GGally)

## ----echo=FALSE, message=FALSE, warning=FALSE,results = "hide"----------------
library(jarbes)
data("healing")
AD = healing[, c("y_c", "n_c", "y_t", "n_t", "Study")]
names(AD) = c("yc", "nc", "yt", "nt", "Study")
m.0 = metarisk(AD)
#summary(m.0)

OR = mean(m.0$BUGSoutput$sims.matrix[,"Odds.pool"])
CI95 = quantile(m.0$BUGSoutput$sims.matrix[,"Odds.pool"], probs = c(0.025, 0.975))


## ----echo=FALSE, message=FALSE, warning=FALSE,results = "hide"----------------

data("healingipd")
IPD = healingipd[, c("healing.without.amp", "PAD", "neuropathy", 
                      "first.ever.lesion", "no.continuous.care", 
                      "male", "diab.typ2", "insulin", "HOCHD", 
                      "HOS", "CRF", "dialysis", "DNOAP", "smoking.ever", 
                      "diabdur", "wagner.class")]


N.IPD = dim(IPD)[1]

 wagner.groups = table(IPD$wagner.class)
IPD.and.wagner = table(IPD$PAD, IPD$wagner.class)


## ----hmr_fit, message=FALSE, warning=FALSE, results = "hide"------------------
#' hmr
#'
set.seed(2022)

hmr.model.1 = hmr(AD,                      # Data frame for Aggregated 
           two.by.two = FALSE,      # If TRUE indicates that the trial results are with names: yc, nc, yt, nt
           dataIPD = IPD,           # Data frame of the IPD 
           re = "sm",               # Random effects model: "sm" scale mixtures 
           link = "logit",          # Link function of the random effects
           sd.mu.1 = 1,             # Scale parameter for the prior of mu.1
           sd.mu.2 = 1,             # Scale parameter for the prior of mu.2 
           sd.mu.phi = 1,           # Scale parameter for the prior of mu.phi 
           sigma.1.upper = 5,       # Upper bound of the prior of sigma.1  
           sigma.2.upper = 5,       # Upper bound of the prior of sigma.2
           sigma.beta.upper = 5,    # Upper bound of the prior of sigma.beta
           sd.Fisher.rho = 1.25,    # Scale parameter for the prior of rho under Fisher's transformation 
           df.estimate = TRUE,      # If TRUE the degrees of freedom are estimated
           df.lower = 3,            # Lower bound of the df's prior
           df.upper = 10,           # Upper bound of the df's prior
           nr.chains = 2,           # Number of MCMC chains
           nr.iterations = 10000,   # Total number of iterations
           nr.adapt = 1000,        # Number of iteration discarded for burnin period
           nr.thin = 1)             # Thinning rate

## -----------------------------------------------------------------------------
summary(hmr.model.1, digits= 2)

## ----hmr1, fig.cap = "Conflict of evidence analysis. Left panel: Prior to posterior sensitivity analysis of bias mean between the RCTs and the IPD-RWD. Right panel: posterior distribution of the outliers detection weights." , echo=TRUE, warning=FALSE, results='hide',message=FALSE, fig.width = 7, fig.height = 5----
# Diagnostic plot
# Analysis of conflict of evidence ................
diagnostic(hmr.model.1, study.names = AD$Study, 
    title.plot.mu.phi = "Prior to Posterior\nSensitivity Analysis",
    title.plot.weights = "Outlier Detection",
    post.p.value.cut = 0.1,
           lwd.forest = 1, shape.forest = 1,
    size.forest = 0.4)


## ----hmr2, fig.cap = "Forest plot of posteriors of regression coefficients of the IPD-RWD. The most relevant risk factors identified in this analysis were: the classification of Wagner (1-2 vs. 3-4-5.)and PAD (no vs. yes)." , echo=TRUE, warning=FALSE, results='hide', message=FALSE, fig.width = 6, fig.height = 7----

attach.jags(hmr.model.1, overwrite = TRUE)

# Figure: Posterior distribution of the regression coefficients IPD
# Forest plot for the 95% posterior intervals of the regression coefficients

# Variable names
var.names = names(IPD[-1])
var.names = c("PAD", "neuropathy", "first.ever.lesion", "no.continuous.care", 
  "male", "diab.typ2", "insulin", "HOCHD", "HOS", "CRF", "dialysis", 
  "DNOAP", "smoking.ever", "diabdur", "Wagner.4")

# Coefficient names
v = paste("beta", names(IPD[-1]), sep = ".")

mcmc.x.2 = as.mcmc.rjags(hmr.model.1)

greek.names = paste(paste("beta[",1:15, sep=""),"]", sep="")
par.names = paste(paste("beta.IPD[",1:15, sep=""),"]", sep="")

caterplot(mcmc.x.2,  
          parms = par.names,
          col = "black", lty = 1, 
          labels = greek.names,
          greek = T,
          labels.loc="axis", cex =0.7, 
          style = "plain",reorder = F, denstrip = F)

caterplot(mcmc.x.2,  
          parms = par.names,
          col = "black", lty = 1, 
          labels = var.names,
          labels.loc="above",
          greek = F,
          cex =0.7, 
          style = "plain",reorder = F, denstrip = F,
          add = T)

abline(v=0, lty = 2, lwd = 2)


## ----hmr3, fig.cap = "Summary results of generalizing relative treatment effects: The RCTs' results are displayed as a forest plot. The fitted hierarchical meta-regression model is summarized with the solid  lines representing the posterior median and 95% intervals. The vertical dashed lines represent the location of the posterior means of the baseline risk for different data and study types." , echo=TRUE, warning=FALSE, results='hide',  results='hide', message=FALSE, fig.width = 7, fig.height = 5----

# Generalization of treatment effects 
plot(hmr.model.1,
     x.lim = c(-5, 3),
     y.lim = c(-2, 6),
     x.lab = "Event rate of The Control Group (logit scale)",
     y.lab = "No improvement = Effectiveness -> Improvement",
     title.plot = "HMR: Effectiveness Against Baseline Risk",
     Study.Types = c("AD-RCTs", "IPD-RWD")
     )

## ----hmr4, fig.cap = "Posterior contourns (50%, 75% and 95%) for the effectivenes for subgroups identified in the Hierarchical Meta-Regression analysis. Left panel: Subgroup of patients with PDA. Right panel: Subgroup of patients with Wagner score > 2." , echo=TRUE, warning=FALSE, results='hide', message=FALSE, fig.width = 7, fig.height = 5----
 
# PDA
p.PDA = effect(hmr.model.1, 
       title.plot = "Subgroup with PDA",
       k = 1,            # Regression coefficient
       x.lim = c(-7, 2.5), 
       y.lim = c(-.5, 2.5), 
       y.lab = "Effectiveness (logit scale)",
       x.lab = "Baseline risk (logit scale)",
       kde2d.n= 150, S = 15000,
       color.line = "blue",
       color.hist = "lightblue",
       font.size.title = 8)

# Wanger
p.Wagner = effect(hmr.model.1, 
  title.plot = "Subgroup with Wagner Score > 2",
  k = 15,            # Regression coefficient
  x.lim = c(-7, 2.5), 
  y.lim = c(-.5, 2.5), 
  y.lab = "Effectiveness (logit scale)",
  x.lab = "Baseline risk (logit scale)",
  kde2d.n= 150, S = 15000, 
  color.line = "red", 
  color.hist = "#DE1A1A",    #medium red
  display.probability = FALSE, 
  line.no.effect = 0,
  font.size.title = 8)

grid.arrange(p.PDA, p.Wagner, ncol = 2, nrow = 1)


## ----hmr5, fig.cap = "Posterior contourns (50%, 75% and 95%) for the effectivenes for subgroups identified in the Hierarchical Meta-Regression analysis. Left panel: Subgroup of patients with PDA. Right panel: Subgroup of patients with Wagner score > 2." , echo=TRUE, warning=FALSE, results='hide', message=FALSE, fig.width = 7, fig.height = 5----
 
# PDA
p.PDA = effect(hmr.model.1, 
       title.plot = "Subgroup with PDA",
       k = 1,            # Regression coefficient
  x.lim = c(0, 1), 
  y.lim = c(0, 5), 
  y.lab = "Effectiveness (Odds Ratio)",
  x.lab = "Baseline risk (probability)",
  kde2d.n= 150, S = 15000, 
  color.line = "red", 
  color.hist = "#DE1A1A",    #medium red
  display.probability = TRUE, 
  line.no.effect = 1,
  font.size.title = 8)

# Wanger
p.Wagner = effect(hmr.model.1, 
  title.plot = "Subgroup with Wagner Score > 2",
  k = 15,            # Regression coefficient
  x.lim = c(0, 1), 
  y.lim = c(0, 5), 
  y.lab = "Effectiveness (Odds Ratio)",
  x.lab = "Baseline risk (probability)",
  kde2d.n= 150, S = 15000, 
  color.line = "blue",
  color.hist = "lightblue",
  display.probability = TRUE, 
  line.no.effect = 1,
  font.size.title = 8)

grid.arrange(p.PDA, p.Wagner, ncol = 2, nrow = 1)


