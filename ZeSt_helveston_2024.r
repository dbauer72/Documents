#####################################################
#### Helveston et al. example                  ######
#####################################################

## load data directly 
load(file = "helveston_data_wide.RData")
sink("helveston.txt")

## data structure for mlogit 
library(mlogit)
cars_data <- mlogit.data(data_wide, varying = 4:51, sep = "_", shape = "wide", choice = "choice", id.var ="id")

head(cars_data[1:9,1:9])

# MNL estimation for panel data set 

(form <- as.formula(paste("choice ~",  paste(vars, collapse = " + ") ,  "| 0")))
P <- length(vars)

st <- Sys.time()
pr <- mlogit(form,cars_data, probit=FALSE, panel = TRUE, rpar = c(opCost = 'n'), start = c(true_b,0.05), R=1000)
fin <- Sys.time()
summary(pr)
(time_elaps <- fin - st)

# comparison: probit without panel 
st <- Sys.time()
pr2 <- mlogit(form,cars_data, probit=TRUE, panel = FALSE, start = true_b, R=200)
fin <- Sys.time()
summary(pr2)
(time_elaps <- fin - st)



#######################################
# set up MNP model without mixing ----------------------------------------------------------
#######################################

make_Hb <- function(P) {
  diag(P)
}

make_fb <- function(P) {
  matrix(0, P, 1)
}

make_HO <- function(P_re) {
  if (P_re > 0) {
    ltho <- P_re * (P_re + 1) / 2
    if (P_re>1){
      ind <- fdiag(P_re)+1
    } else {
      ind <- 1
    }
    return(diag(ltho)[, ind, drop = FALSE])
  } else {
    return(matrix(0, 0, 0))
  }
}

make_fO <- function(P_re) {
  if (P_re > 0) {
    ltho <- P_re * (P_re + 1) / 2
    matrix(0, ltho, 1)
  } else {
    return(matrix(0, 0, 0))
  }
}


# re = NULL 
re = NULL
P_re <- length(re)

# load package 
library(Rprobit)

mod_hel <- mod_cl$new(
  alt = 3,
  Hb = make_Hb(P),
  fb = make_fb(P),
  HO = make_HO(P_re),
  fO = make_fO(P_re),
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE
)

mod_hel$fL[1] = 1.264911
mod_hel$fL[4] = 1.264911
mod_hel$fL[6] = 1.264911

mod_hel

# cross sectional data set 
data_wide_CS <- data_wide
data_wide_CS[,"idc"] <-1 
data_wide_CS[,"id"] <-1:5760 

hel_obj <- setup_Rprobit(
  form = form,
  data_raw = data_wide_CS,
  re = NULL,
  id = "id",
  mod = mod_hel
)

hel_obj$control$probit <- TRUE

theta_0 <- c(true_b)
hel_obj$theta_0 <- theta_0
hel_obj$theta <- theta_0

hel_fit <- fit_Rprobit(hel_obj, init_method = "theta")
summary(hel_fit)

hel_fit$info$estimation_time

### now make opCost random 
re = "opCost"
P_re <- length(re)

mod_hel <- mod_cl$new(
  alt = 3,
  Hb = make_Hb(P),
  fb = make_fb(P),
  HO = make_HO(P_re),
  fO = make_fO(P_re),
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE
)

mod_hel$fL[1] = 1.264911
mod_hel$fL[4] = 1.264911
mod_hel$fL[6] = 1.264911

# cross sectional estimation
hel_obj_CS <- setup_Rprobit(
  form = form,
  data_raw = data_wide_CS,
  re = re,
  id = "id",
  mod = mod_hel
)

hel_obj_CS$control$probit <- FALSE
theta_0 <- c(true_b, rep(0.05, P_re))
hel_obj_CS$theta_0 <- theta_0

hel_fit_CS <- fit_Rprobit(hel_obj_CS, init_method = "theta", cml_pair_type = 0)
summary(hel_fit_CS)

hel_fit_CS$info$estimation_time


## now use panel data. 
hel_obj <- setup_Rprobit(
  form = form,
  data_raw = data_wide,
  re = re,
  id = "id",
  mod = mod_hel
)

hel_obj$control$probit <- FALSE
theta_0 <- c(true_b, rep(0.05, P_re))
hel_obj$theta_0 <- theta_0

hel_fit0 <- fit_Rprobit(hel_obj, init_method = "theta", cml_pair_type = 0)
summary(hel_fit0)

hel_fit0$info$estimation_time

hel_fit1 <- fit_Rprobit(hel_obj, init_method = "theta", cml_pair_type = 1)
summary(hel_fit1)

hel_fit1$info$estimation_time

hel_obj$theta_0 <- hel_fit1$theta 
hel_obj$theta <- hel_fit1$theta 

hel_fit2 <- fit_Rprobit(hel_obj, init_method = "theta", cml_pair_type = 0)
summary(hel_fit2)

hel_fit2$info$estimation_time


plot(hel_fit0, margin = 3)

# compare estimates 
beta <- coef(pr)

cbind(beta[1:16]/beta[2],hel_fit$theta[1:16]/hel_fit$theta[2],hel_fit0$theta[1:16]/hel_fit0$theta[2],hel_fit1$theta[1:16]/hel_fit1$theta[2])




##########################
## Latent class model 
##########################

mode_LC <- mod_latclass_cl$new(
  Hb   = matrix(0,16,1),
  fb   = matrix(0,16,1),
  HO   = matrix(1,1,0),
  fO   =  matrix(1,1,1)*0.01,
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE
)

mode_LC$fL[1] = 1.264911
mode_LC$fL[4] = 1.264911
mode_LC$fL[6] = 1.264911

mode_LC$Hb[1,1] <- 1
mode_LC$fb <- as.matrix(hel_fit0$theta[1:16],ncol=1)
mode_LC$fb[1] <- 0 

mode_LC$num_class <- 3

mode_LC$Hb
mode_LC$fb

theta_0LC <- c(-0.5,-0.3,-0.1,0,0)

hel_obj$mod <- mode_LC

hel_obj$control$probit <- FALSE
hel_obj$control$el = 1
hel_obj$control$control_weights$cml_pair_type <- 0
hel_obj$theta_0 <- theta_0LC
hel_obj$theta <- theta_0LC

hel_LC <- fit_LC_Rprobit(hel_obj, init_method = "theta",cml_pair_type = 1)
summary(hel_LC)

hel_LC$info$estimation_time
plot(hel_LC)

pr_LC <- predict(hel_LC)

#########################################
## non-parametric estimation 
#########################################
params = matrix(seq(from=-0.9,to=.5,by= 0.02),nrow=1)

mode_np <- mod_nonpara_cl$new(alt  = 3,
                              Hb   = matrix(0,16,1),
                              fb   = matrix(0,16,1),
                              HO   = matrix(0,1,0),
                              fO   =  matrix(1,1,1)*0.01,
                              HL   = matrix(0,6,0),
                              fL   = matrix(0,6,1),
                              ordered = FALSE
)

mode_np$fL[1] = 1.264911
mode_np$fL[4] = 1.264911
mode_np$fL[6] = 1.264911

mode_np$Hb[1,1] <- 1
mode_np$fb <- as.matrix(hel_fit0$theta[1:16],ncol=1)
mode_np$fb[1] <- 0 

mode_np$set_grid_points(params)

hel_obj_np <- setup_Rprobit(
  form = form,
  data_raw = data_wide,
  re = re,
  id = "id",
  mod = mode_np
)

hel_obj_np$control$probit <- FALSE
hel_obj_np$control$control_weights$cml_pair_type <- 1 # full pairwise CML function to be used.
hel_obj_np$mod <- mode_np

Rpro_np <- fit_nonpara_Rprobit(hel_obj_np, init_method = "random", control_nlm = NULL, cml_pair_type = 1)

probs_np <- Rpro_np$theta
plot(params,probs_np,col="red")

# add density for LC model. 
plot(hel_LC)
points(params,probs_np*111,col="red")

sink()
