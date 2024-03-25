## PAE data anlaysis ##
## Alex Stringer & Tugba Akkaya Hocagil
## 2023

## Packages ##
pkgs <- c(
	"tidyverse",
	"mgcv"
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE, quietly = TRUE)
  }
}

## Load the semibmd package ----
if (!require("semibmd",character.only = TRUE,quietly = TRUE)) {
	stop("semibmd package not found. Install with the command:\n\n remotes::install_github('awstringer1/semibmd')\n")
}


## Set paths ##
# Set this to a base directory of your choice.
basepath <- "~/work/projects/benchmark-dose"
# Data should be called "pae-data.csv" and stored in basepath/data
datapath <- file.path(basepath,"data")
resultspath <- file.path(basepath,"figures")
stopifnot(dir.exists(datapath))
if (!dir.exists(resultspath)) dir.create(resultspath)

## Shouldn't need to change anything beyond this point ----

## Data ##

dat <- read_csv(file.path(datapath,"pae-data.csv"))
dat$CogScoreScaled <- scale(dat$CogScore)
dat$Cohort <- factor(dat$CohortIndicator)

## Model ##
mod <- benchmark_dose_tmb(
  monosmooths = list(s(Alch3, bs="bs",k=50)),
  smooths = list(s(PS3, by=Cohort)),
  linearterms = NULL,
  data=dat,
  exposure = 'Alch3',
  response = 'CogScoreScaled',
  x0=0,
  p0=.025,
  BMR=.01,
  maxitr = 10,
  bayes_boot = 1e05,
  verbose = FALSE,
  scale_data = FALSE # Already pre-scaled
) # Takes about 5 minutes on 2021 M1 Max Macbook Pro

summary(mod)
get_all_bmdl(mod)
# relative computation times
comptimes <- get_computation_times(mod)
boottime <- comptimes['bmd_samples']
reltimes <- 100*round(comptimes/boottime,4)
total_noboot <- comptimes['model'] + comptimes['posterior_samples'] + comptimes['bmdl_score']

# plots for paper
plotdat <- plot(mod,plot=FALSE)
plotdatmono <- plotdat$mono

pdf(file.path(resultspath,"dose-response-plot.pdf"),width=5,height=5)
with(dat,plot(CogScoreScaled~Alch3,
              main = "Dose-response curve, cog. score vs drinks/day",
              xlab = "log(1+drinks per day)",
              ylab = "Standardized cognitive function score",
              pch = "."))

with(plotdatmono,lines(x,estimate))
with(plotdatmono,lines(x,lower,lty='dashed'))
with(plotdatmono,lines(x,upper,lty='dashed'))

abline(v=get_bmd(mod)[1],lty='solid')
abline(v=get_all_bmdl(mod)['score'],lty='dashed')
abline(v=get_all_bmdl(mod)['delta'],lty='dotted')
abline(v=get_all_bmdl(mod)['bmdl_bayes'],lty='dotdash')
dev.off()

# posterior samples of the BMD
bmd_samps <- mod$info$bmd_samps
modapprox <- get_approximations(mod)
bmdmean <- get_bmd(mod)[1]
bmdsd <- sqrt(modapprox[1]) / abs(modapprox[2])

pdf(file.path(resultspath,"bmdhistogram.pdf"),width=5,height=5)
hist(bmd_samps,
     freq = FALSE,
     main = "Posterior dist. of the BMD, 100,000 samples",
     xlab = "Sampled benchmark dose",
     ylab = "Density",
     breaks = 100
     )
curve(dnorm(x,bmdmean,bmdsd),add=TRUE)
abline(v=get_all_bmdl(mod)['delta'],lty='dotted')
abline(v=get_all_bmdl(mod)['score'],lty='dashed')
abline(v=get_all_bmdl(mod)['bmdl_bayes'],lty='dotdash')
dev.off()

# residual plot
pdf(file.path(resultspath, "residualfittedplot.pdf"), width = 5, height = 5)
plot(residuals(mod) ~ fitted(mod),
  main = "Residuals vs fitted values, PAE dose-response model",
  xlab = "Fitted Values",
  ylab = "Residuals"
)
abline(h = 0, lty = "dashed")
dev.off()

pdf(file.path(resultspath, "residualqqplot.pdf"), width = 5, height = 5)
qqnorm(residuals(mod),
  main = "QQ-plot of residuals, PAE dose-response model",
  xlab = "N(0,1) Quantiles",
  ylab = "Empirical Quantiles of Residuals"
)
qqline(residuals(mod))

dev.off()
