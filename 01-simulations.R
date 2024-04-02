## Simulations ###
# For paper: Semi parametric benchmark dose analysis with monotone additive models
# Section 4
# Alex Stringer
# 2023

## Load Libraries ----

pkgs <- c(
	"tidyverse",
	"scales",
	"mgcv",
	"scam",
	"Matrix",
	"parallel",
	"Rcpp",
	"RcppEigen",
	"nloptr"
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE, quietly = TRUE)
  }
}
# Set number of cores
options(mc.cores = parallel::detectCores())

## Load the semibmd package ----
if (!require("semibmd",character.only = TRUE,quietly = TRUE)) {
	stop("semibmd package not found. Install with the command:\n\n remotes::install_github('awstringer1/semibmd')\n")
}

# Set sizes for saving figures
GGPLOTTEXTSIZE <- 16
GGPLOTWIDTH <- 7
GGPLOTHEIGHT <- 10

## Set Paths ----

# Set these to what you want to call the results
simname <- "202403-28-v1.RData"
tablename <- "202403-28-v1.csv"

# Change this to where you want the results saved
globalpath <- "~/work/projects/benchmark-dose"
stopifnot(dir.exists(globalpath))
resultspath <- file.path(globalpath,"results")
stopifnot(dir.exists(resultspath))
plotpath <- file.path(globalpath,"figures")
stopifnot(dir.exists(plotpath))

## Below this point, no changes are necessary ----

## Set parameters ##
n <- c(200,500,1000)
sigma <- c(.1,.2,.5)

BMR <- .01 # The actual BMR, not a multiple of p0
scale <- c(.1,.5,1,2,5)
numsims <- 100000
bayes_boot <- 1000

GLOBAL_P0 <- .01 # Set this for the duration
set_parameters <- function(scale=1,sigma=1,p0=GLOBAL_P0,BMR=.025,plot_it=FALSE,knots=10) {
  xmin <- 0
  xmax <- 1
  # True function f
  stopifnot(scale>0)
  f <- function(x) exp(-x*scale)
  finv <- function(x) -log(x)/scale
  
  x0 <- 0
  stopifnot(x0>=xmin)
  tau0 <- f(x0) + sigma * qnorm(p0) # Response level
  stopifnot(abs(pnorm(tau0,f(x0),sigma) - p0) < sqrt(.Machine$double.eps)) # Definition of tau
  
  A <- qnorm(p0+BMR) - qnorm(p0)
  # Benchmark dose
  xb <- finv(f(x0) - sigma*A)
  stopifnot(abs(pnorm(tau0,f(xb),sigma) - (p0+BMR)) < sqrt(.Machine$double.eps)) # Definition of benchmark dose
  
  if (plot_it) {
    curve(f,xmin,xmax,ylim=c(0,1))
    abline(v=xb,lty='dashed')
    cat("xb = ",xb,"; f'(xb)/sigma = ",numDeriv::grad(f,xb)/sigma,".\n",sep="")
  } else {
    list(
      xmin = xmin,xmax = xmax,
      sigma = sigma,
      p0 = p0,BMR = BMR,
      f = f,
      x0 = x0,
      tau0 = tau0,
      A = A,
      xb = xb,
      knots=knots
    )
  }
}

## Simulate Data ##
simulate_data <- function(n,params) {
  xcov <- with(params,seq(xmin,xmax,length.out=n)) # Deterministic covariate
  # Response
  y <- with(params,rnorm(n,f(xcov),sigma))
  # y <- y - mean(y)
 data.frame(y=y,x=xcov)
}


dosim <- function(lst) {
  n <- lst[1]
  scale <- lst[4]
  sigma <- lst[2]
  BMR <- lst[3]
  knots <- 10
  numsamps <- lst[6]
  sim <- lst[7]
  cat("Simulation:",sim,"of",length(nsimlist),"| n =",n,"| scale =",scale,"| sigma =",sigma,"| BMR =",BMR,"| \n")

  params <- set_parameters(scale=scale,sigma = sigma,BMR=BMR,knots = knots)
  simdata <- simulate_data(n,params)

  mod <- tryCatch(benchmark_dose_tmb(monosmooths = list(s(x,bs='bs',k=knots)),
			    smooths = NULL,
			    linearterms = NULL,
			   data = simdata,
			   exposure = 'x',
			   response = 'y',
			   x0 = 0,
			   p0 = GLOBAL_P0,
			   BMR = BMR,
			   verbose = FALSE,
			   eps = 1e-06,
			   maxitr = 100,
			   bayes_boot = numsamps,
			   scale_data = FALSE
	),error = function(e) e)
  if (inherits(mod,'condition')) return(NULL)

  bmdest <- tryCatch(get_bmd(mod)[1],error = function(e) e)
  if (inherits(bmdest,'condition')) return(NULL)

  if (get_errors(mod)) return(NULL)

	# Bias #
	estimate <- get_bmd(mod)[1]
	true <- params$xb
	bias <- estimate - true
	# Coverage #
	all_bmdl <- get_all_bmdl(mod)
	covr_score <- as.logical(all_bmdl['score'] <= true)
	covr_delta <- as.logical(all_bmdl['delta'] <= true)
	covr_bayes <- as.logical(all_bmdl['bmdl_bayes'] <= true)

	# Computation times
	comptimes <- get_computation_times(mod)
	time_model <- unname(comptimes['model'])
	time_bmd <- unname(comptimes['bmd_estimate'])
	time_delta <- unname(comptimes['bmdl_delta'])
	time_score <- unname(comptimes['bmdl_score'])
	time_bayes <- unname(comptimes['bmd_samples'])

	c('n' = n,'scale' = scale,'sigma' = sigma,'BMR' = BMR,
			'bias' = bias,'estimate' = estimate,'true' = true,
			'bmdl_score' = unname(all_bmdl['score']),
			'bmdl_delta' = unname(all_bmdl['delta']),
			'bmdl_bayes' = unname(all_bmdl['bmdl_bayes']),
			'covr_score' = covr_score,'covr_delta' = covr_delta,
			'covr_bayes' = covr_bayes,
			'time_model' = time_model,'time_bmd' = time_bmd,'time_delta' = time_delta,'time_score' = time_score,
			'time_bayes' = time_bayes,
			'rel_time_score' = time_score/time_delta,
			'rel_time_bayes' = time_bayes/time_delta
	)
}

# Do the simulations

simstodoframe <- expand.grid(n=n,sigma=sigma,BMR=BMR,scale=scale,sim=1:numsims,bayes_boot=bayes_boot)
nsimlist <- vector(mode='list',length=nrow(simstodoframe))
for (i in 1:length(nsimlist)) nsimlist[[i]] <- c(as.numeric(simstodoframe[i, ]),i)

cat("Doing",length(nsimlist),"simulations...\n")
tm <- Sys.time()
sims <- mclapply(nsimlist,dosim)
simtime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
cat("Done sims, they took",simtime,"seconds.\n")


cat("Processing simulations...\n")
processed_sims <- dplyr::bind_rows(sims) %>% as_tibble()
cat("Finished processing simulations.\n")

cat("Saving sims...\n")
save(processed_sims,file=file.path(resultspath,simname))
load(file.path(resultspath,simname))
cat("Done saving sims. File:",file.path(resultspath,simname),"\n")
cat("Saving output table...\n")

outputtable <- processed_sims %>%
	filter_all(all_vars(is.finite(.))) %>%
	group_by(n,scale,sigma) %>%
	summarize(
			bias_mean = 100*mean(bias,na.rm=TRUE),
		  bias_sd = 100*(sd(bias,na.rm=TRUE)/sqrt(n())),
		  covr_delta_mean = mean(covr_delta,na.rm=TRUE),
		  covr_delta_sd = sqrt(covr_delta_mean*(1-covr_delta_mean)/n()),
		  covr_score_mean = mean(covr_score,na.rm=TRUE),
		  covr_score_sd = sqrt(covr_score_mean*(1-covr_score_mean)/n()),
		  covr_bayes_mean = mean(covr_bayes,na.rm=TRUE),
		  covr_bayes_sd = sqrt(covr_bayes_mean*(1-covr_bayes_mean)/n()),
		  rel_time_score_mean = mean(rel_time_score,na.rm=TRUE),
		  rel_time_score_sd = sd(rel_time_score,na.rm=TRUE)/sqrt(n()),
		  rel_time_bayes_mean = mean(rel_time_bayes,na.rm=TRUE),
		  rel_time_bayes_sd = sd(rel_time_bayes,na.rm=TRUE)/sqrt(n()),
		  percent_failed_mean = (numsims-n())/numsims,
		  percent_failed_sd = sqrt(percent_failed_mean*(1-percent_failed_mean) / n()),   
		  percent_useless_delta_mean = mean(bmdl_delta <= 0,na.rm=TRUE),
		  percent_useless_delta_sd = sqrt(percent_useless_delta_mean*(1-percent_useless_delta_mean)/n())  
	) %>%
	mutate(covr_delta_mean = 100*covr_delta_mean,
	       covr_delta_sd = 100*covr_delta_sd,
	       covr_score_mean = 100*covr_score_mean,
	       covr_score_sd = 100*covr_score_sd,
	       covr_bayes_mean = 100*covr_bayes_mean,
	       covr_bayes_sd = 100*covr_bayes_sd,
	       percent_failed_mean = 100*percent_failed_mean,
	       percent_failed_sd = 100*percent_failed_sd,
	       percent_useless_delta_mean = 100*percent_useless_delta_mean,
	       percent_useless_delta_sd = 100*percent_useless_delta_sd
	)


readr::write_csv(outputtable,file = file.path(resultspath,tablename))

cat("Finished saving output table, file: ",file.path(resultspath,tablename),"\n")

## Plots ##
# Make plots of the bias/CP, easier for the reader to visualize.
scale_lab <- function(s) paste0("f(x) = exp(-",s,"x)")
sigma_lab <- function(s) paste0(expression(sigma~"="), s)

coverageplot <- outputtable %>%
	select(n, scale, sigma, contains("covr")) %>%
	pivot_longer(contains("covr"), names_to = "temp", values_to = "coverage") %>%
	mutate(type = str_match(temp, "delta|score|bayes")[ ,1],
				 whichthing = str_match(temp, "mean|sd")[ ,1]) %>%
	pivot_wider(id_cols = all_of(c("n", "scale", "sigma", "type")), names_from = whichthing, values_from = coverage) %>%
	ggplot(aes(x = factor(n), y = mean, group = factor(type), colour = factor(type))) + 
	theme_bw() + 
	facet_grid(scale ~ sigma, labeller = label_bquote(cols = sigma~"="~.(sigma), rows = "f(x) = exp(-"*.(scale)*"x)")) + 
	geom_line() +
	geom_line(aes(y = mean + 2 * sd), linetype = "dashed") +
	geom_line(aes(y = mean - 2 * sd), linetype = "dashed") +
	labs(title = expression("Empirical Coverage: proportion of"~widehat(x)[l]*""<=x[b]),
			 x = "Sample size, n", y = "Empirical Coverage, %",
			 colour = expression(widehat(x)[l])) +
	coord_cartesian(ylim = 100*c(0.8,1)) + 
	scale_y_continuous(breaks = seq(80, 100, by = 5)) +
	# scale_colour_brewer(type = "qual", palette = 2, labels = c("bayes" = "Bootstrap", "score" = "Pivot", "delta" = "Delta")) + 
	# scale_colour_brewer(type = "seq", palette = "Greys", labels = c("bayes" = "Bootstrap", "score" = "Pivot", "delta" = "Delta")) + 
	scale_colour_manual(
		labels = c("score" = "Pivot", "bayes" = "Bootstrap", "delta" = "Delta"),
		breaks = c("score", "bayes", "delta"),
		values = c("#000000", "#636363", "#bdbdbd")
	) +
	geom_hline(yintercept = 97.5, linetype = "dotted") + 
	theme(legend.position = "bottom", text = element_text(size = GGPLOTTEXTSIZE))

ggsave(file.path(resultspath, "coverageplot.pdf"), plot = coverageplot, width = GGPLOTWIDTH, height = GGPLOTHEIGHT)

# For the supplement: boxplots of the lengths of the CIs (BMD - BMDL) and the bias
lengthtable <- processed_sims %>%
	select(n, scale, sigma, estimate, bmdl_score, bmdl_delta, bmdl_bayes) %>%
	pivot_longer(contains("bmdl"), names_to = "method", values_to = "bmdl") %>%
	mutate(method = str_replace(method, "bmdl_", ""), length = estimate - bmdl)

lengthplotgood <- lengthtable %>%
	filter(length < .1) %>%
	ggplot(aes(x = method, y = length)) + 
	theme_bw() + 
	facet_grid(scale ~ n, labeller = label_bquote(cols = n~"="~.(n), rows = "f(x) = exp(-"*.(scale)*"x)")) + 
	geom_boxplot() +
	scale_x_discrete(labels = c("bayes" = "Boot", "delta" = "Delta", "score" = "Pivot")) + 
	labs(title = "Lengths of one-sided BMD confidence intervals",
			x = "Method", y = expression("Length,"~widehat(x)[b] - widehat(x)[l])) +
	theme(text = element_text(size = GGPLOTTEXTSIZE))


lengthplotbad <- lengthtable %>%
	filter(length < 1e06) %>%
	ggplot(aes(x = method, y = length)) + 
	theme_bw() + 
	facet_grid(scale ~ n, labeller = label_bquote(cols = n~"="~.(n), rows = "f(x) = exp(-"*.(scale)*"x)")) + 
	geom_boxplot() +
	scale_x_discrete(labels = c("bayes" = "Boot", "delta" = "Delta", "score" = "Pivot")) + 
	scale_y_continuous(labels = scales::comma_format()) +
	labs(title = "Lengths of one-sided BMD confidence intervals",
			x = "Method", y = expression("Length,"~widehat(x)[b] - widehat(x)[l])) +
	theme(text = element_text(size = GGPLOTTEXTSIZE))

ggsave(file.path(resultspath, "lengthplotgood.pdf"), plot = lengthplotgood, width = GGPLOTWIDTH, height = GGPLOTHEIGHT)
ggsave(file.path(resultspath, "lengthplotbad.pdf"), plot = lengthplotbad, width = GGPLOTWIDTH, height = GGPLOTHEIGHT)

biasplot <- processed_sims %>%
	ggplot(aes(x = factor(n), y = bias)) + 
	theme_bw() + 
	facet_grid(scale ~ sigma, labeller = label_bquote(cols = sigma~"="~.(sigma), rows = "f(x) = exp(-"*.(scale)*"x)")) + 
	geom_boxplot() +
	labs(title = expression("Empirical Bias,"~widehat(x)[b]~-~x[b]),
			 x = "Sample size, n", y = "Empirical Bias") +
	theme(text = element_text(size = GGPLOTTEXTSIZE))

biasplotzoom <- processed_sims %>%
filter(abs(bias) < .05) %>%
	ggplot(aes(x = factor(n), y = bias)) + 
	theme_bw() + 
	facet_grid(scale ~ sigma, labeller = label_bquote(cols = sigma~"="~.(sigma), rows = "f(x) = exp(-"*.(scale)*"x)")) + 
	geom_boxplot() +
	labs(title = expression("Empirical Bias,"~widehat(x)[b]~-~x[b]),
			 x = "Sample size, n", y = "Empirical Bias") +
	theme(text = element_text(size = GGPLOTTEXTSIZE))

	

ggsave(file.path(resultspath, "biasplot.pdf"), plot = biasplot, width = GGPLOTWIDTH, height = GGPLOTHEIGHT)
ggsave(file.path(resultspath, "biasplotzoom.pdf"), plot = biasplotzoom, width = GGPLOTWIDTH, height = GGPLOTHEIGHT)
