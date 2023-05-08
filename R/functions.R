## FUNCTIONS #########

# sensRICLPM input needs to be a data frame with all observations of one variable X in the first few columns and the other variable Y in the second half of columns
# data: (example:) if there are three waves of data, column 1-3 should contain the observed scores on X, and column 4-6 should contain the observed scores on Y)
# ME_prop_x: the steps in the proportions of measurement error in the observed scores of X that will be included in the sensitivity analysis
# ME_prop_y: the steps in the proportions of measurement error in the observed scores of Y that will be included in the sensitivity analysis
# constraints: can be set to "none", "stationarity", "lagged", "within", "residual", and "equal ME variances"


#' Title
#'
#' @param data
#' @param ME_prop_x
#' @param ME_prop_y
#' @param combinations
#' @param constraints
#'
#' @return
#' @export
#'
#' @examples data(amotivationexample)
sensRICLPM <- function(data,
                       ME_prop_x = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9),
                       ME_prop_y = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9),
                       # var_name_x = "X",
                       # var_name_y = "Y",
                       combinations = "all",
                       constraints = "none"){

  library(magrittr)
  library(dplyr)
  library(purrr)
  library(lavaan)
  library(ggplot2)
  # add input checks

  if(combinations == "all"){
    conditions <- expand.grid(MEpropx = ME_prop_x , MEpropy = ME_prop_y) # all possible combinations of ME proportions
  }else if(combinations == "simple"){
    conditions <- data.frame(MEpropx = ME_prop_x,
                             MEpropy = ME_prop_y) # only equal ME proportions
  }
  # add number of waves and sample size to conditions
  conditions$waves <- timepoints(data)
  conditions$samplesize <- nrow(data)

  # consistent names for observed variables
  var_names <- sapply(c("X", "Y"), paste0, 1:conditions[["waves"]][[1]])
  names(data) <- var_names

  conditions$constraints <- constraints

  # Compute ME variances from proportions
  conditions$MEx <- purrr::map(conditions$MEpropx, .f = compute_ME_var_x, data, constraints)
  conditions$MEy <- purrr::map(conditions$MEpropy, .f = compute_ME_var_y, data, constraints)


  # generate lavaan syntax for all combinations of ME proportions
  syntaxes <- lapply(X = asplit(conditions, MARGIN = 1),
                     FUN = create_lavaan,
                     constraints = constraints
  )
  # extract all model syntaxes
  conditions$syntax <- purrr::map(syntaxes, function(i){
    pluck(i, "est_synt")
  })


  message(rlang::format_error_bullets(c(
    i = "Fitting lavaan models..."
  )))

  # fit all models
  fits <- suppressWarnings(purrr::map(conditions$syntax, lavaan::lavaan, data))


  message(rlang::format_error_bullets(c(
    i ="Saving estimates...")))

  # check convergence
  conditions$converged <- suppressWarnings(
    purrr::map_lgl(fits, function(i){
      lavInspect(i, what ="converged")}))

  # check admissibility
  conditions$admissible <- suppressWarnings(
    purrr::map_lgl(fits, function(i){
      lavInspect(i, what ="post.check")}))


  # save parameter estimates
  conditions$results <- purrr::map(fits, function(i) parameterEstimates(i, remove.nonfree = TRUE))

  # return(conditions)

  ## PLOTS ##
  # create dataframe for estimates of lagged parameters


  parameterlist <- list(

    alphas <- purrr::map(conditions$results, function(i) # extract all alphas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~", startsWith(rhs, prefix = "wX"))
    ),
    betas <- purrr::map(conditions$results, function(i)  # extract all betas (cross-lagged)
      dplyr::filter(i, startsWith(lhs, prefix = "wY"), op == "~", startsWith(rhs, prefix = "wX"))
    ),
    gammas <- purrr::map(conditions$results, function(i)  # extract all gammas (cross-lagged)
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~", startsWith(rhs, prefix = "wY"))
    ),
    deltas <- purrr::map(conditions$results, function(i)  # extract all deltas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "wY"), op == "~", startsWith(rhs, prefix = "wY"))
    ),
    RIX <- purrr::map(conditions$results, function(i) # extract all random intercept variance estimates for X
      dplyr::filter(i, startsWith(lhs, prefix = "RI_X"), op == "~~", startsWith(rhs, prefix = "RI_X"))
    ),
    RIY <- purrr::map(conditions$results, function(i) # extract all random intercept variance estimates for Y
      dplyr::filter(i, startsWith(lhs, prefix = "RI_Y"), op == "~~", startsWith(rhs, prefix = "RI_Y"))
    ),
    wX <- purrr::map(conditions$results, function(i) # extract all within variance estimates for X
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~~", startsWith(rhs, prefix = "wX"))
    ),
    wY <- purrr::map(conditions$results, function(i) # extract all within variance estimates for Y
      dplyr::filter(i, startsWith(lhs, prefix = "wY"), op == "~~", startsWith(rhs, prefix = "wY"))
    ),
    RIcov <- purrr::map(conditions$results, function(i) # extract all alphas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "RI_X"), op == "~~", startsWith(rhs, prefix = "RI_Y"))
    ),
    wcov <- purrr::map(conditions$results, function(i) # extract all alphas (autoregressive)
      dplyr::filter(i, startsWith(lhs, prefix = "wX"), op == "~~", startsWith(rhs, prefix = "wY"))
    ))

  # create list containing estimates for all lagged effects, SEs & CIs (grouped by wave)

  plotframe <- create_outputtable(parameterlist, conditions)
  # plotframe_s <- plotframe %>% filter(parameter %in% c("alpha", "beta", "gamma","delta"))

  # ggplots

  if(constraints == "stationarity"){
    p1 <- plotframe %>%
      filter(converged== "TRUE" & parameter %in% c("alpha", "beta", "gamma", "delta"))%>%
      ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+
      geom_line(linewidth = .5)+
      geom_point()+
      geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
      xlab("Proportion of measurement error in x")+
      ylab("Estimate")+
      labs(color = "Proportion of measurement error in y")+
      theme_bw()+
      theme(legend.position="bottom")+
      ggtitle("Sensitivity of lagged effect estimates to measurement error")+
      facet_wrap(~parameter)

    p2 <-  plotframe %>%
      filter(admissible== "TRUE"& parameter %in% c("alpha", "beta", "gamma", "delta"))%>%
      ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+ # make them shades of grey
      geom_line(linewidth = .5)+
      geom_point()+
      geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
      xlab("Proportion of measurement error in X")+
      ylab("Estimate")+
      labs(color = "Proportion of measurement error in Y")+
      theme_bw()+
      theme(legend.position="bottom")+
      ggtitle("Sensitivity of lagged effect estimates to measurement error")+
      facet_wrap(~parameter)

  }else{
    p1 <- plotframe %>%
      filter(converged== "TRUE" & parameter %in% c("alpha", "beta", "gamma", "delta"))%>%
      ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+
      geom_line(linewidth = .5)+
      geom_point()+
      geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
      xlab("Proportion of measurement error in X")+
      ylab("Estimate")+
      labs(color = "Proportion of measurement error in Y")+
      theme_bw()+
      theme(legend.position="bottom")+
      ggtitle("Sensitivity of lagged effect estimates to measurement error")+
      facet_grid(wave ~ parameter)

    p2 <-  plotframe %>%

      filter(admissible== "TRUE" & parameter %in% c("alpha", "beta", "gamma", "delta"))%>%
      ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+ # make them shades of grey
      geom_line(linewidth = .5)+
      geom_point()+
      scale_fill_brewer(palette="Set1")+
      #  geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
      xlab("Proportion of measurement error in X")+
      ylab("Estimate")+
      labs(color = "Proportion of measurement error in Y")+
      theme_bw()+
      theme(legend.position="bottom")+
      ggtitle("Sensitivity of lagged effect estimates to measurement error")+
      facet_grid(wave ~ parameter)

  }
  result <- list(outputtable = plotframe,
                 p1 = p1,
                 p2 = p2)


}



# create list containing estimates for all lagged effects, SEs & CIs (grouped by wave)

# parframe <- as.data.frame(purrr::map(parameterlist, create_plotdata))

# some reshaping
#  parframe <- dplyr::select(parframe, -starts_with("wave."), -starts_with("id.")) # remove duplicate columns
# plotframe <- cbind(conditions, parframe) # join condition table and data frame containing lagged estimates

#  plotframe <- dplyr::rename(plotframe, # consistent variable names for following pivot_longer
#                            est_1.0 = est_1,
#                           est_2.0 = est_2,
#                          se_1.0 = se_1,
#                         se_2.0 = se_2,
#                        pvalue_1.0 = pvalue_1,
#                       pvalue_2.0 = pvalue_2,
#                      cilower_1.0 = cilower_1,
#                     cilower_2.0 = cilower_2,
#                    ciupper_1.0 = ciupper_1,
#                   ciupper_2.0 = ciupper_2)

# plotframe<- pivot_longer(plotframe, !c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged,
#                                       admissible, results, id, syntax), names_to = c(".value", "wave", "parameter"),
#                        names_pattern = "(.+)_(\\d).(.+)")  # make longer format


# cosmetics
# plotframe$parameter <- dplyr::recode(plotframe$parameter, "0" = "alpha", "1" = "beta", "2" = "gamma", "3" = "delta")# works
#  plotframe$wave <- as.numeric(plotframe$wave)

#return(plotframe)

# ggplots

#  if(constraints == "stationarity"){
#    p1 <- plotframe %>% dplyr::filter(converged== "TRUE")%>%
#      ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+
#      geom_line(linewidth = .5)+
#      geom_point()+
#      geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
#      xlab("Proportion of measurement error in x")+
#      ylab("Estimate")+
#      labs(color = "Proportion of measurement error in y")+
#      theme_bw()+
#      theme(legend.position="bottom")+
#      ggtitle("Sensitivity of lagged effect estimates to measurement error")+
#      facet_wrap(~parameter)

#   p2 <-  plotframe %>% dplyr::filter(admissible== "TRUE")%>%
##    ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+ # make them shades of grey
#    geom_line(linewidth = .5)+
#    geom_point()+
#    geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
#    xlab("Proportion of measurement error in x")+
##    ylab("Estimate")+
#    labs(color = "Proportion of measurement error in y")+
#   theme_bw()+
#   theme(legend.position="bottom")+
#   ggtitle("Sensitivity of lagged effect estimates to measurement error (admissible results)")+
#   facet_wrap(~parameter)

#  }else{
#   p1 <- plotframe %>% dplyr::filter(converged== "TRUE")%>%
#      ggplot(aes(x = MEpropx, y = value, color = factor(MEpropy)))+
#     geom_line(linewidth = .5)+
#    geom_point()+
#   geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
#  xlab("Proportion of measurement error in x")+
#      ylab("Estimate")+
#     labs(color = "Proportion of measurement error in y")+
#    theme_bw()+
#   theme(legend.position="bottom")+
#  ggtitle("Sensitivity of lagged effect estimates to measurement error")+
# facet_grid(wave ~ parameter)

#   p2 <-  plotframe %>% dplyr::filter(admissible== "TRUE")%>%
#    ggplot(aes(x = MEpropx, y = est, color = factor(MEpropy)))+ # make them shades of grey
#   geom_line(linewidth = .5)+
#  geom_point()+
# geom_pointrange(aes(ymin = cilower, ymax = ciupper))+
#      xlab("Proportion of measurement error in x")+
#     ylab("Estimate")+
#    labs(color = "Proportion of measurement error in y")+
#   theme_bw()+
#  theme(legend.position="bottom")+
# ggtitle("Sensitivity of lagged effect estimates to measurement error (admissible results)")+
# facet_grid(wave ~ parameter)
#}


# return(p1)

#  result <- list(outputtable = plotframe,
#                p1 = p1,
#               p2 = p2)


#}





# timepoints() saves the number of measurement waves
timepoints<- function(d){
  nwaves <- length(diag(cov(d)))/2
  return(nwaves)
}




# compute_ME_var_x() computes values that ME variances of x need to be fixed to from different proportions of ME & total #variance in the data
compute_ME_var_x <- function(ME_prop_x, data, constraints){
  waves <- timepoints(data)

  if (constraints == "equal_ME_variances"){ # ME variances constrained to be equal over waves
    ME_var_x_mat <-matrix(
      (sum(diag(cov(data[1:waves]))) / waves) %*% ME_prop_x, # computes ME variance from average total variance in x
      nrow = length(diag(cov(data[1:waves]))), byrow=T)
  }else{
    ME_var_x_mat <-matrix(
      diag(cov(data[1:waves])) %x% ME_prop_x, # computes ME variance for x at each wave
      nrow=length(diag(cov(data[1:waves]))), byrow=T)
  }
}


# compute_ME_var_y() computes values that ME variances of y need to be fixed to from different proportions of ME & total #variance in the data
compute_ME_var_y <- function(ME_prop_y, data, constraints){
  waves <- timepoints(data)
  if (constraints == "equal_ME_variances"){ # ME variances constrained to be equal over waves
    matrix(
      (sum(diag(cov(data[-(1:waves)]))) / waves) %*% ME_prop_y, # computes ME variance from average total variance in y
      nrow = length(diag(cov(data[-(1:waves)]))), byrow=T)
  }else{
    matrix(
      diag(cov(data[-(1:waves)])) %x% ME_prop_y,  # computes ME variance for y at each wave
      nrow=length(diag(cov(data[-(1:waves)]))), byrow=T)
  }
}



# ME_var_pop() computes the value that population ME variances need to be so that there will be certain proportions of ME
ME_var_pop <- function(ME.p, RI.var, W.var){
  (ME.p%*%RI.var + ME.p%*%W.var)/(1 - ME.p)
}



# create_plotdata() creates a long format table from a lavaan.data.frame object, including estimates for lagged effects, SEs, pvalues, and confidence intervals grouped by wave,
#create_plotdata <- function(parameter){

#  plotdata <- cbind(
#    est <- purrr::map(parameter, function(i){
#      pluck(i, "est") # extract estimates
#    }),
#    se <- purrr::map(parameter, function(i){
#      pluck(i, "se")  # extract SEs
#    }),
#    pvalue <- purrr::map(parameter, function(i){
#      pluck(i, "pvalue")  # extract pvalues
#    }),
#    ci_lower <- purrr::map(parameter, function(i){
#      pluck(i, "ci.lower") # extract lower bounds of confidence intervals
##    }),
#    ci_higher <- purrr::map(parameter, function(i){
#      pluck(i, "ci.upper") # extract upper bounds of confidence intervals
#    })
#  )
#  plotdata <- as.data.frame(plotdata)
#  names(plotdata) <- c("est", "se", "pvalue", "cilower", "ciupper") # add column names
#  plotdata$id <- 1:nrow(plotdata) # add id
#  plotdata <- tidyr:: unnest_wider(plotdata, col = !id, names_sep = "_") # unnest wave estimates
# plotdata <- pivot_longer(plotdata, !id, names_to = c(".value", "wave"), names_sep = "_") # make longer format
# return(plotdata)
#}


# create_outputtable() creates the outputtable from a lavaan.data.frame (list) object
#create_outputtable <- function(parameter, conditions){
#
#  parframe <- as.data.frame(purrr::map(parameter, create_plotdata)) # grab relevant estimates

# consitent column names:
# parnames <- c("alpha", "beta", "gamma", "delta", "RIX", "RIY", "wX", "wY", "RIcov", "wcov")
#  estnames <- c("est", "se", "pvalue", "cilower", "ciupper")
#  colnames(parframe) <-  as.vector(outer(estnames, parnames, paste, sep="_")) # add column names

#  plotframe <- cbind(conditions, parframe) # join condition table and data frame containing estimates

#  plotframe <- tidyr:: unnest_wider(plotframe,!c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged,
##                                                 admissible, results, syntax), names_sep = ".") # unnest wave estimates

#
# outputtable<- pivot_longer(plotframe, names_to = c(".value", "parameter", "wave"), !c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged,
#                                                                                      admissible, results, syntax),
#                          names_pattern = "(.+)_(.+).(\\d)")  # make longer format

#}


# create_outputtable() creates the outputtable from a lavaan.data.frame (list) object
create_outputtable <- function(parameter, conditions){

  parframe <- as.data.frame(purrr::map(parameter, create_plotdata)) # grab relevant estimates

  # consitent column names:
  parnames <- c("alpha", "beta", "gamma", "delta", "RIX", "RIY", "wX", "wY", "RIcov", "wcov")
  estnames <- c("est", "se", "pvalue", "cilower", "ciupper")
  colnames(parframe) <-  as.vector(outer(estnames, parnames, paste, sep="_")) # add column names

  plotframe <- cbind(conditions, parframe) # join condition table and data frame containing estimates

  plotframe <- tidyr:: unnest_wider(plotframe,!c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged,
                                                 admissible, results, syntax), names_sep = ".") # unnest wave estimates


  outputtable<- pivot_longer(plotframe, names_to = c(".value", "parameter", "wave"), !c(MEpropx, MEpropy, waves, samplesize, constraints, MEx, MEy, converged,
                                                                                        admissible, results, syntax),
                             names_pattern = "(.+)_(.+).(\\d)")  # make longer format

}

# create_plotdata() creates a long format table from a lavaan.data.frame object, including estimates for lagged effects, SEs, pvalues, and confidence intervals grouped by wave,
create_plotdata <- function(parameter){

  plotdata <- cbind(
    est <- purrr::map(parameter, function(i){
      pluck(i, "est") # extract estimates
    }),
    se <- purrr::map(parameter, function(i){
      pluck(i, "se")  # extract SEs
    }),
    pvalue <- purrr::map(parameter, function(i){
      pluck(i, "pvalue")  # extract pvalues
    }),
    ci_lower <- purrr::map(parameter, function(i){
      pluck(i, "ci.lower") # extract lower bounds of confidence intervals
    }),
    ci_higher <- purrr::map(parameter, function(i){
      pluck(i, "ci.upper") # extract upper bounds of confidence intervals
    })
  )
  plotdata <- as.data.frame(plotdata)

}


# extra function
ME_var_pop2 <- function(ME.p, RI.var, W.var, waves){
  matrix(
    (ME.p%*%RI.var + ME.p%*%W.var)/(1 - ME.p), nrow=waves, byrow = T)
}


# create_plotdata() creates a long format table from a lavaan.data.frame object, including estimates for lagged effects, SEs, pvalues, and confidence intervals grouped by wave,
create_plotdata <- function(parameter){

  plotdata <- cbind(
    est <- purrr::map(parameter, function(i){
      pluck(i, "est") # extract estimates
    }),
    se <- purrr::map(parameter, function(i){
      pluck(i, "se")  # extract SEs
    }),
    pvalue <- purrr::map(parameter, function(i){
      pluck(i, "pvalue")  # extract pvalues
    }),
    ci_lower <- purrr::map(parameter, function(i){
      pluck(i, "ci.lower") # extract lower bounds of confidence intervals
    }),
    ci_higher <- purrr::map(parameter, function(i){
      pluck(i, "ci.upper") # extract upper bounds of confidence intervals
    })
  )
  plotdata <- as.data.frame(plotdata)

}


# extra function
ME_var_pop2 <- function(ME.p, RI.var, W.var, waves){
  matrix(
    (ME.p%*%RI.var + ME.p%*%W.var)/(1 - ME.p), nrow=waves, byrow = T)
}

# lavaan.R

# create_lavaan() produces model syntax for RI-CLPMs with ME variances fixed to different values
create_lavaan <- function(conditions,
                          constraints
)
{

  # Generate default variable names
  name_var <- LETTERS[24:25] # X and Y


  # Create matrix of names for observed variable, within, measurement error, and between components
  name_obs <- sapply(name_var, paste0, 1:conditions[["waves"]])
  name_within <- sapply(name_var, function(x) {
    paste0("w", x, 1:conditions[["waves"]])
  })


  name_RI <- paste0("RI_", name_var)


  # Create estimation parameter table
  est_tab <- rbind(
    lav_RI(conditions = conditions, name_RI, name_obs),
    est_RI_var(name_RI),
    est_RI_cov(name_RI),
    est_within(name_within, name_obs, constraints),
    est_lagged(conditions = conditions, name_within, constraints),
    est_within_var1(name_within, constraints),
    est_within_cov1(name_within, constraints),
    est_within_var2(conditions = conditions, name_within, constraints),
    est_within_cov2(conditions = conditions, name_within, constraints),
    fix_MEx(conditions = conditions, name_obs),
    fix_MEy(conditions = conditions, name_obs),
    if (constraints == "stationarity"){
      rbind(
        lav_within_cov(conditions = conditions),
        lav_stationarity(conditions = conditions)
      )
    }
  )
  rownames(est_tab) <- NULL

  # Create lavaan syntax

  est_synt <- paste0( # Paste over parameters
    paste0( # Paste over columns
      est_tab[, 1],
      est_tab[, 2],
      est_tab[, 3],
      est_tab[, 4],
      est_tab[, 5]
    ),
    collapse = "\n"
  )


  # Create condition list with extra element space
  list(
    waves = conditions[["waves"]],
    ME_proportion_in_x = conditions[["MEpropx"]],
    ME_proportion_in_y = conditions[["MEpropy"]],
    est_synt = est_synt,
    est_tab = est_tab,
    estimates = NA,
    uncertainty = NA,
    errors = NA,
    not_converged = NA,
    inadmissible = NA
  )
}


lav_RI <- function(conditions, name_RI, name_obs) {
  lhs <- rep(name_RI, each = conditions[["waves"]])
  op <- rep("=~", times = 2 * conditions[["waves"]])
  pv <- rep("1", times = 2 * conditions[["waves"]])
  con <- rep("*", times = 2 * conditions[["waves"]])
  rhs <- c(unlist(name_obs))
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_RI_var <- function(name_RI) {
  lhs <- rhs <- name_RI
  op <- "~~"
  con <- "*"
  pv <- "NA"
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}

est_RI_cov <- function(name_RI) {
  lhs <- name_RI[[1]]
  rhs <- name_RI[[2]]
  op <- "~~"
  con <- "*"
  pv <- "NA"
  free <- TRUE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_within <- function(name_within, name_obs, constraints) {
  lhs <- c(name_within)
  op <- "=~"
  if (constraints == "stationarity") {
    pv <- "NA"
    free <- TRUE
  } else {
    pv <- "1"
    free <- FALSE
  }
  con <- "*"
  rhs <- c(name_obs)
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_lagged <- function(conditions, name_within, constraints) {
  lhs <- rep(c(t(name_within))[-(1:2)], each = 2)
  op <- "~"
  con <- "*"
  if (constraints == "lagged" || constraints == "within" ||
      constraints == "stationarity"){
    pv <- c("a", "b", "c", "d") # Labels for constraints
    free <- TRUE
    rhs <- c(apply(name_within[-conditions[["waves"]], ], 1, rep, times = 2))
  } else {
    pv <- "NA"
    free <- TRUE
    rhs <- c(apply(name_within[-conditions[["waves"]], ], 1, rep, times = 2))
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}



est_within_var1 <- function(name_within, constraints) {
  lhs <- rhs <- c(t(name_within[1, ]))
  op <- "~~"
  con <- "*"
  if (constraints == "stationarity") {
    pv <- "1"
    # pv <-     pv <- c(
    #  paste0("varx"),
    # paste0("vary"))
    free <- TRUE

  } else {
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}


est_within_cov1 <- function(name_within, constraints) {
  lhs <- name_within[1, 1]
  rhs <- name_within[1, 2]
  op <- "~~"
  con <- "*"
  if (constraints == "stationarity") { # Label
    pv <- "cov1"
    free <- TRUE
  } else { # Freely estimate
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}


est_within_var2 <- function(conditions, name_within, constraints) {
  lhs <- rhs <- c(name_within[-1, ])
  op <- "~~"
  con <- "*"
  if (constraints == "residuals" || constraints == "within") { # Constrain over time
    pv <- rep(c("rvarx", "rvary"), each = (conditions[["waves"]] - 1))
    free <- TRUE
  } else if (constraints == "stationarity") {
    pv <- c(
      paste0("rvarx", 2:conditions[["waves"]]),
      paste0("rvary", 2:conditions[["waves"]]))
    free <- TRUE
  }
  else { # Freely estimate
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}


est_within_cov2 <- function(conditions, name_within, constraints){
  lhs <- name_within[-1, 1]
  rhs <- name_within[-1, 2]
  op <- "~~"
  con <- "*"
  if (constraints == "residuals" || constraints == "within") { # Constrain over time
    pv <- "rcov"
    free <- TRUE
  } else if (constraints == "stationarity") { # Label
    pv <- paste0("rcov", 2:conditions[["waves"]])
    free <- TRUE
  }
  else { # Freely estimate
    pv <- "NA"
    free <- TRUE
  }
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE))
}


fix_MEx <- function(conditions, name_obs){
  lhs <- rhs <- name_obs[, "X"]
  op <- "~~"
  con <- "*"
  pv <- unlist(conditions[["MEx"]])
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE))
}


fix_MEy <- function(conditions, name_obs){
  lhs <- rhs <- name_obs[, "Y"]
  op <- "~~"
  con <- "*"
  pv <- unlist(conditions[["MEy"]])
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free, stringsAsFactors = FALSE))
}



# computes covariance of within components at each wave
lav_within_cov <- function(conditions) {
  lhs <- paste0("cov", 2:(conditions[["waves"]] - 1))
  op <- ":="
  pv <- con <- ""
  rhs1 <- paste0("a*c + b*d + a*d*cov", 1:(conditions[["waves"]] - 2))
  rhs2 <- paste0(" + b*c*cov", 1:(conditions[["waves"]]- 2))
  rhs3 <- paste0(" + rcov", 2:(conditions[["waves"]] - 1))
  rhs <- paste0(rhs1, rhs2, rhs3)
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}

# constrain residual variances

lav_stationarity <- function(conditions){   # double check this
  lhs <- c(
    paste0("rvarx", 2:conditions[["waves"]]),
    paste0("rvary", 2:conditions[["waves"]])
  )
  op <- "=="
  rhs <- c(
    paste0("1 - (a^2 + b^2 + 2*a*b*cov", 1:(conditions[["waves"]] - 1), ")"),
    paste0("1 - (c^2 + d^2 + 2*c*d*cov", 1:(conditions[["waves"]] - 1), ")")
  )
  pv <- con <- ""
  free <- FALSE
  return(cbind.data.frame(lhs, op, pv, con, rhs, free,
                          stringsAsFactors = FALSE
  ))
}




