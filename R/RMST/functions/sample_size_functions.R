#compares sample sizes for superiority trials and noninferiority trials, RMST-difference vs. LRT. For different scales but same shapes.

compare_sample_size <- function(scale_t = 0.01,
                                scale_c = 0.02,
                                shape_t = 1,
                                shape_c = 1,
                                accrual_time = 0.00000001,
                                follow_up_time = 12,
                                horizon,
                                sides = 1,
                                alpha = 0.025,
                                power = 0.8,
                                margin = 0,
                                npsurvSS_RMSTD = T,
                                npsurvSS_LRT_unrestricted = T,
                                npsurvSS_LRT_restricted = T,
                                RMSTdesign_closed_form = T,
                                RMSTdesign_simulation = T,
                                powerRMST = T,
                                ssrmst = T,
                                survmixer = T
){
  total_time <- accrual_time + follow_up_time
  # initialize df for storing sample sizes
  sample_sizes_df <- data.frame(matrix(ncol = 8))
  rownames(sample_sizes_df) <- "n or power"
  colnames(sample_sizes_df) <- c("npsurvSS RMSTD", "npsurvSS LRT unrestricted", "npsurvSS LRT restricted", "RMSTdesign closed form", "RMSTdesign simulation", "powerRMST power", "ssrmst power", "survmixer")
  
  if(npsurvSS_RMSTD | npsurvSS_LRT_unrestricted | npsurvSS_LRT_restricted){
    npsurvSS_arm_t <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_t, surv_shape = shape_t,
                                 loss_scale = 0.0000000000000000000000001, loss_shape = 1, total_time = total_time)
    npsurvSS_arm_c <- create_arm(size = 1, accr_time = accrual_time, surv_scale = scale_c, surv_shape = shape_c,
                                 loss_scale = 0.0000000000000000000000001, loss_shape = 1, follow_time = follow_up_time)
    if(npsurvSS_RMSTD){
      sample_sizes_df["npsurvSS RMSTD"] <- round(size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "rmst difference", milestone = horizon), power = 0.8, alpha = alpha, sides = sides)['n'])
    }
    if(npsurvSS_LRT_unrestricted){
      sample_sizes_df["npsurvSS LRT unrestricted"] <- round(size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "weighted logrank"), power = 0.8, alpha = alpha, sides = sides)['n'])
    }
    if(npsurvSS_LRT_restricted){
      npsurvSS_arm_t$follow_time <- npsurvSS_arm_c$follow_time <- horizon
      npsurvSS_arm_t$total_time <- npsurvSS_arm_c$total_time <- horizon + accrual_time
      sample_sizes_df["npsurvSS LRT restricted"] <- round(size_two_arm(npsurvSS_arm_t, npsurvSS_arm_c, test = list(test = "weighted logrank"), power = 0.8, alpha = alpha, sides = sides)['n'])
    }
    print(sample_sizes_df)
  }
  if(RMSTdesign_closed_form & margin == 0){
    sample_sizes_df["RMSTdesign closed form"] <- round(RMSTpow(survdefT = survdef(haz = scale_t), survdefC = survdef(haz = scale_c), k1 = accrual_time, k2 = follow_up_time, power = power, tau = horizon, sim = F, method = 'tau_star', 
                                                               two.sided = sides == 2, alpha = alpha)$n)
    print(sample_sizes_df)
  }
  if(RMSTdesign_simulation & margin == 0){
    sample_sizes_df["RMSTdesign simulation"]  <- round(RMSTpow(survdefT = survdef(haz = scale_t), survdefC = survdef(haz = scale_c), k1 = accrual_time, k2 = follow_up_time, power = power, tau = horizon, sim = T, method = 'tau_star', 
                                                               two.sided = sides == 2, alpha = alpha, M = 10)$n)
    print(sample_sizes_df)
  }
  if(powerRMST & npsurvSS_RMSTD){
    sample_sizes_df["powerRMST power"] <-  powerRMST(n = sample_sizes_df[["npsurvSS RMSTD"]], ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=margin, one_sided_alpha=alpha/sides, seed=NULL, M = 10)$power
    print(sample_sizes_df)
  }
  if(ssrmst & npsurvSS_RMSTD){
    sample_sizes_df["ssrmst power"] <- ssrmst(ac_number = sample_sizes_df[["npsurvSS RMSTD"]], ac_period = accrual_time, tot_time = total_time, tau = horizon, scale0 = 1/scale_c, scale1 = 1/scale_t, margin=0, 
                                              one_sided_alpha=alpha/sides, ntest=10)$power1
    print(sample_sizes_df)
  }
  if(survmixer & margin == 0){
    sample_sizes_df["survmixer"] <-  round(survm_samplesize(
      ascale0_r = 1/scale_c,
      ascale0_nr = 1/scale_c,
      ascale1_r = 1/scale_t,
      ascale1_nr = 1/scale_t,
      delta_p = 0,
      p0 = 1,
      m0_r = 1,
      m0_nr = 1,
      diffm_r = 1,
      diffm_nr = 1,
      S0_r = 1,
      S0_nr = 1,
      diffS_r = 1,
      diffS_nr = 1,
      Delta_r = 1,
      Delta_nr = 1,
      ascale_cens = 1000,
      tau = horizon,
      bshape0 = 1,
      bshape1 = 1,
      alpha = alpha/sides,
      beta = 0.2,
      set_param = 0
    )[1, 2])
  }
  print(sample_sizes_df)
  return(sample_sizes_df)
}

# Calculate design parameters delta, sigma2, and tsigma2
# Calculate design parameters delta, sigma2, and tsigma2
calc_design <- function(arm0, arm1, test) {
  
  p1  <- arm1$size / (arm0$size + arm1$size)
  p0  <- 1 - p1
  tsigma2 <- NULL
  
  # survival difference and ratio
  if (grepl("survival", test$test)) {
    
    if (! "milestone" %in% names(test)) {
      stop(paste("Please provide milestone for ",
                 test$test,
                 ".",
                 sep=""),
           call.=F)
    }
    milestone <- test$milestone
    surv0 <- psurv(milestone, arm0, lower.tail=F)
    surv1 <- psurv(milestone, arm1, lower.tail=F)
    
    if (test$test == "survival difference") { # survival difference
      delta   <- surv0 - surv1
      sigma2  <- sigma2j_surv(arm0, milestone) / p0 + sigma2j_surv(arm1, milestone) / p1
    } else if (test$test == "survival ratio") { # survival ratio
      delta   <- log(surv0) - log(surv1)
      sigma2  <- sigma2j_cumh(arm0, milestone) / p0 + sigma2j_cumh(arm1, milestone) / p1
    } else {
      stop("Please specify valid survival contrast.", call.=F)
    }
    
  } else if (grepl("rmst", test$test)) { # rmst difference and ratio
    
    if (! "milestone" %in% names(test)) {
      stop(paste("Please provide milestone for ",
                 test$test,
                 ".",
                 sep=""),
           call.=F)
    }
    milestone <- test$milestone
    rmst0 <- deltaj_rmst(milestone, arm0)
    rmst1 <- deltaj_rmst(milestone, arm1)
    
    if (test$test == "rmst difference") { # rmst difference
      delta   <- rmst0 - rmst1
      sigma2  <- sigma2j_rmst(arm0, milestone) / p0 + sigma2j_rmst(arm1, milestone) / p1
    } else if (test$test == "rmst ratio") { # rmst ratio
      delta   <- log(rmst0) - log(rmst1)
      sigma2  <- sigma2j_rmst(arm0, milestone) / p0 / rmst0^2 +
        sigma2j_rmst(arm1, milestone) / p1 / rmst1^2
    } else {
      stop("Please specify a valid rmst contrast.", call.=F)
    }
    
  } else if (grepl("percentile", test$test)) { # percentile difference or ratio
    
    if (! "percentile" %in% names(test)) {
      stop(paste("Please provide percentile for ",
                 test$test,
                 ".",
                 sep=""),
           call.=F)
    }
    percentile <- test$percentile
    perc0 <- qsurv(percentile, arm0)
    perc1 <- qsurv(percentile, arm1)
    
    if (test$test == "percentile difference") { # percentile difference
      delta   <- perc0 - perc1
      sigma2  <- sigma2j_perc(arm0, percentile) / p0 + sigma2j_perc(arm1, percentile) / p1
    } else if (test$test == "percentile ratio") { # percentile ratio
      delta   <- log(perc0) - log(perc1)
      sigma2  <- sigma2j_perc(arm0, percentile) / p0 / perc0^2 +
        sigma2j_perc(arm1, percentile) / p1 / perc1^2
    } else {
      stop("Please specify a valid percentile contrast.", call.=F)
    }
  } else if (test$test == "hazard ratio") {
    delta     <- c(arm0$surv_shape * log( arm1$surv_scale / arm0$surv_scale ))[1] # log hazard-ratio
    sigma2    <- sigma2_clhr(arm0, arm1)
  } else if (test$test == "weighted logrank") {
    # weight
    if (! "weight" %in% names(test)) {
      test$weight <- "1"
    }
    # delta and sigma2
    if (! "mean.approx" %in% names(test)) {
      test$mean.approx <- "asymptotic"
    }
    delta     <- delta_wlr(arm0, arm1, test$weight, test$mean.approx)
    sigma2    <- sigma2_wlr(arm0, arm1, test$weight, test$mean.approx)
    # tsigma2
    if (! "var.approx" %in% names(test)) {
      test$var.approx <- "1"
    }
    if (test$var.approx == "1") {
      tsigma2 <- sigma2
    } else {
      tsigma2 <- tsigma2_wlr(arm0, arm1, test$weight, test$var.approx)
    }
  } else {
    stop("Please specify a valid test.", call.=F)
  }
  if(is.null(tsigma2)) {
    tsigma2 <- sigma2
  }
  
  return(list(delta=delta, sigma2=sigma2, tsigma2=tsigma2))
  
}

power_two_arm <- function(arm0,
                          arm1,
                          test=list(test="weighted logrank"),
                          alpha=0.025,
                          sides=1) {
  
  if (! inherits(test[[1]], "list")) { # one test to perform
    
    n       <- arm0$size + arm1$size
    design  <- calc_design(arm0, arm1, test)
    out     <- stats::pnorm((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(n) * design$delta) /
                              sqrt(design$tsigma2),
                            lower.tail=F) +
      (sides==2) * stats::pnorm((sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) - sqrt(n) * design$delta) /
                                  sqrt(design$tsigma2),
                                lower.tail=F)
    return(out)
    
  } else { # multiple tests to perform
    
    out <- c()
    for (i in 1:length(test)) {
      label = ifelse("label" %in% names(test[[i]]), test[[i]]$label, i)
      out <- rbind(out, c(label, power_two_arm(arm0, arm1, test[[i]], alpha, sides)))
    }
    out <- data.frame(out)
    names(out) <- c("test", "power")
    return(out)
    
  }
  
}

#' Sample size
#'
#' Calculate required sample size and expected number of events for a
#' two-arm survival study.
#' @param arm0  object of class 'arm'.
#' @param arm1  object of class 'arm'.
#' @param test  list or list of lists. Each list must contain at minimum
#'   the key 'test' describing the type of statistical test. Default test
#'   is the "weighted logrank". Kaplan-Meier based tests ("survival difference",
#'   "survival ratio", "rmst difference", "rmst ratio", "percentile difference",
#'   and "percentile ratio") require the user to define an additional key,
#'   either the desired 'milestone' or 'percentile'. The weighted log-rank test
#'   does not require additional keys. However, user may choose which weight function
#'   ("1"=unweighted, "n"=Gehan-Breslow, "sqrtN"=Tarone-Ware, "FH_[a]_[b]"=
#'   Fleming-Harrington with p=a and q=b) and which approximation for the
#'   large-sample mean ("asymptotic", "generalized schoenfeld", "event driven",
#'   "freedman", "rubinstein") and variance ("1", "block[ randomization]", "simple[ randomization]") 
#'   they wish to use. Default choice is 'weight'="1", 'mean.approx'="asymptotic", and 'var.approx'="1".
#'   For more details regarding the different mean and variance approximations
#'   for the weight log-rank test, please see Yung and Liu (2020). If there are multiple lists, 
#'   then users may provide a 'label' for each list to be displayed in the output.
#' @param power 1 - type 2 error rate
#' @param alpha type 1 error rate
#' @param sides 1=1-sided test, 2=2-sided test
#' @return
#'   \item{n0}{sample size for \code{arm0}}
#'   \item{n1}{sample size for \code{arm1}}
#'   \item{n}{total sample size}
#'   \item{d0}{expected number of events for \code{arm0}}
#'   \item{d1}{expected number of events for \code{arm1}}
#'   \item{d}{total expected number of events; can be used to convert a time-driven
#'   trial to an event-driven trial.}
#' @seealso \code{\link{create_arm}} for creating an object of class 'arm'.
#' @references Yung, G and Liu, Y. (2020). Sample size and power for the weighted
#' log-rank test and Kaplan-Meier based tests with allowance for non-proportional
#' hazards. \emph{Biometrics} 76(3):939-950.
#' @examples
#' arm0 <- create_arm(size=120, accr_time=6, surv_scale=0.05, loss_scale=0.005, follow_time=12)
#' arm1 <- create_arm(size=120, accr_time=6, surv_scale=0.03, loss_scale=0.005, follow_time=12)
#' size_two_arm(arm0, arm1)
#' size_two_arm(arm0, arm1, list(test="weighted logrank",
#'   weight="n",
#'   mean.approx="generalized schoenfeld",
#'   var.approx="block"))
#' size_two_arm(arm0, arm1, list(test="survival difference", milestone=12))
#' size_two_arm(arm0, arm1, list(test="rmst ratio", milestone=12))
#' size_two_arm(arm0, arm1, list(test="percentile difference", percentile=0.25))
#' size_two_arm(arm0, arm1, list(
#'   list(test="weighted logrank", label="Logrank"),
#'   list(test="survival difference", milestone=12, label="12-month survival difference")))

#' @export
size_two_arm <- function(arm0,
                         arm1,
                         test=list(test="weighted logrank"),
                         power=0.8,
                         alpha=0.025,
                         sides=1) {
  
  if (! inherits(test[[1]], "list")) { # one test to perform
    
    # sample size for 1-sided test
    p1      <- arm1$size / (arm0$size + arm1$size)
    p0      <- 1 - p1
    design  <- calc_design(arm0, arm1, test)
    out     <- ( sqrt(design$sigma2) * stats::qnorm(1 - alpha / sides) + sqrt(design$tsigma2) * stats::qnorm(power) )^2 /
      design$delta^2 *
      c(p0, p1)
    # out     <- ceiling(out) # rounding
    
    # refine sample size for 2-sided test
    if (sides==2) {
      i         <- 1
      arm0$size <- out[1]
      arm1$size <- out[2]
      cont      <- T
      while (cont) {
        temp <- power_two_arm(arm0, arm1, test, alpha, sides)
        if (temp > power) {
          i         <- i + 1
          arm0$size <- out[1] - i * p0
          arm1$size <- out[2] - i * p1
        } else {
          out   <- out - (i-1) * c(p0, p1)
          cont  <- F
        }
      }
    }
    
    out <- c(out, # n0, n1
             sum(out), # n
             out[1] * prob_event(arm0), # d0
             out[2]  *prob_event(arm1), # d1
             out[1] * prob_event(arm0) + out[2] * prob_event(arm1)) # d
    # if (out[4] %% 1 < out[5] %% 1) { # rounding
    #   out[4] <- floor(out[4])
    #   out[5] <- ceiling(out[5])
    # } else {
    #   out[4] <- ceiling(out[4])
    #   out[5] <- floor(out[5])
    # }
    names(out) <- c("n0", "n1", "n", "d0", "d1", "d")
    return(out)
    
  } else { # multiple tests to perform
    
    out <- c()
    for (i in 1:length(test)) {
      label = ifelse("label" %in% names(test[[i]]), test[[i]]$label, i)
      out <- rbind(out, c(label, size_two_arm(arm0, arm1, test[[i]], power, alpha, sides)))
    }
    out <- data.frame(out)
    names(out)[1] <- "test"
    return(out)
    
  }
  
}


powerRMST <-
  function(n, ac_period, tot_time, tau, scale0, scale1, shape=1, margin=0, allocation1=0.5, one_sided_alpha=0.025, seed=NULL, M=20000){
    ac_rate = n / ac_period
    
    if (ac_period > tot_time){
      n0 = round(ac_rate*tot_time*(1-allocation1))
      n1 = round(ac_rate*tot_time*allocation1)
    }
    if (ac_period <= tot_time){
      n0 = round(ac_rate*ac_period*(1-allocation1))
      n1 = round(ac_rate*ac_period*allocation1)
    }
    
    ###--- test (main part) ---
    answer     = NULL
    check      = NULL
    event_arm0 = NULL
    event_arm1 = NULL
    
    if (!is.null(seed)){
      set.seed(seed)
    }
    
    for (w in 1:M){
      
      ##-- data frame --
      E             = rweibull(n0, shape, scale0)
      C             = tot_time - runif(n0, 0, ac_period)
      time          = pmin(E,C)
      status        = as.numeric(E<=C)
      arm           = rep(0,n0)
      data0         = data.frame(time, status, arm)
      ind0          = data0$status==1
      
      
      E             = rweibull(n1, shape, scale1)
      C             = tot_time - runif(n1, 0, ac_period)
      time          = pmin(E,C)
      status        = as.numeric(E<=C)
      arm           = rep(1,n1)
      data1         = data.frame(time, status, arm)
      ind1          = data1$status==1
      
      data   = rbind(data0, data1)
      data   = data[data$time>0, ]
      
      data$status[which(data$time>tau)] <- 0
      data$time[which(data$time>tau)] <- as.numeric(tau)
      
      res.km <- summary(survfit(Surv(time, status)~arm, data=data), rmean=tau, print.rmean=T)
      
      rmstC <- res.km$table[1,5]
      rmstC_SE <- res.km$table[1,6]
      rmstE <- res.km$table[2,5]
      rmstE_SE <- res.km$table[2,6]
      
      z <- qnorm(1-one_sided_alpha)
      lower <- rmstE-rmstC - z*sqrt(rmstC_SE^2+rmstE_SE^2)
      
      if (lower >= margin){
        answer[w] = 1
      } else {
        answer[w] = 0
      }
    }
    
    ###--- power ---
    power = sum(answer)/M
    
    
    ###--- output ---
    out = matrix(0,1,3)
    
    out[1,1] = n0+n1
    out[1,2] = n0
    out[1,3] = n1
    
    
    rownames(out) = c("Sample size")	
    colnames(out) = c("Total", "arm0", "arm1")
    
    Z = list()
    
    Z$result      = out
    Z$power       = power
    Z$ac_rate     = ac_rate
    Z$ac_period   = ac_period
    Z$tot_time    = tot_time
    Z$margin      = margin
    Z$tau         = tau
    
    Z
  }
