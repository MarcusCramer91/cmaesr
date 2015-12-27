#' @title Stopping condition: maximal iterations.
#'
#' @description Stop on maximal number of iterations.
#'
#' @param max.iter [integer(1)]\cr
#'   Maximal number of iterations.
#'   Default is \code{100}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnMaxIters = function(max.iter = 100L) {
  assertInt(max.iter, na.ok = FALSE)
  force(max.iter)
  return(makeStoppingCondition(
    name = "maxIter",
    message = sprintf("MaxIter: reached maximal number of iterations/generations %i.", max.iter),
    stop.fun = function(envir = parent.frame()) {
      return(envir$iter > max.iter)
    }
  ))
}

#' @title Stopping condition: indefinite covariance matrix.
#'
#' @description Stop if covariance matrix is not positive definite anymore.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnIndefCovMat = function() {
  return(makeStoppingCondition(
    name = "indefCovMat",
    message = "Covariance matrix is not numerically positive definite.",
    stop.fun = function(envir = parent.frame()) {
      e.values = envir$e$values
      return(any(e.values <= sqrt(.Machine$double.eps) * abs(max(e.values))))
    }
  ))
}

#' @title Stopping condition: optimal params.
#'
#' @description Stop if euclidean distance of parameter is below
#' some tolerance value.
#'
#' @param opt.param [\code{numeric}]\cr
#'   Known optimal parameter settings.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptParam = function(opt.param, tol = 1e-8) {
  assertNumeric(opt.param, any.missing = FALSE, all.missing = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.param)
  force(tol)
  return(makeStoppingCondition(
    name = "optParamTol",
    message = sprintf("Optimal parameters approximated nicely (gap < %.2f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(sqrt(sum(envir$best.param - opt.param)^2) < tol)
    }
  ))
}

#' @title Stopping condition: optimal objective value.
#'
#' @description Stop if best solution is close to optimal objective value.
#'
#' @param opt.value [\code{numeric(1)}]\cr
#'   Known optimal objective function value.
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \eqn{1e^{-8}}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOptValue = function(opt.value, tol = 1e-8) {
  assertNumber(opt.value, na.ok = FALSE)
  assertNumber(tol, lower = 0, na.ok = FALSE, finite = TRUE)
  force(opt.value)
  force(tol)
  return(makeStoppingCondition(
    name = "optValTol",
    message = sprintf("Optimal function value approximated nicely (gap < %.10f).", tol),
    stop.fun = function(envir = parent.frame()) {
      return(abs(envir$best.fitness - opt.value) < tol)
    }
  ))
}

#' @title Stopping condition: maximal time.
#'
#' @description Stop if maximal running time budget is reached.
#'
#' @param budget [\code{integer(1)}]\cr
#'   Time budget in seconds.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnTimeBudget = function(budget) {
  assertInt(budget, na.ok = FALSE, lower = 1L)
  force(budget)
  return(makeStoppingCondition(
    name = "timeBudget",
    message = sprintf("Time budget of %i [secs] reached.", budget),
    stop.fun = function(envir = parent.frame()) {
      return(difftime(Sys.time(), envir$start.time, units = "secs") > budget)
    }
  ))
}

#' @title Stopping condition: maximal funtion evaluations.
#'
#' @description Stop if maximal number of function evaluations is reached.
#'
#' @param max.evals [\code{integer(1)}]\cr
#'   Maximal number of allowed function evaluations.
#' @return [\code{cma_stopping_condition}]
#' @export
stopOnMaxEvals = function(max.evals) {
  assertInt(max.evals, na.ok = FALSE, lower = 1L)
  force(max.evals)
  return(makeStoppingCondition(
    name = "maxEvals",
    message = sprintf("Maximal number of %i function evaluations reached.", max.evals),
    stop.fun = function(envir = parent.frame()) {
      return(envir$n.evals >= max.evals)
    }
  ))
}

#' @title Stopping condition: low standard deviation.
#'
#' @description Stop if the standard deviation falls below a tolerance value
#' in all coordinates?
#'
#' @param tol [\code{integer(1)}]\cr
#'   Tolerance value.
#' @return [\code{cma_stopping_condition}]
#' @export
#FIXME: default value is 10^(-12) * sigma. Here we have no access to the sigma value.
stopOnTolX = function(tol = 10^(-12)) {
  assertInt(tol, na.ok = FALSE)
  force(tol)
  return(makeStoppingCondition(
    name = "tolX",
    message = sprintf("Standard deviation below tolerance in all coordinates."),
    stop.fun = function(envir = parent.frame()) {
      return(all(envir$D < tol) && all((envir$sigma * envir$p.c) < tol))
    }
  ))
}

#' @title Stopping condition: principal axis.
#'
#' @description Stop if addition of 0.1 * sigma in a principal axis
#' direction does not change mean value.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnNoEffectAxis = function() {
  return(makeStoppingCondition(
    name = "noEffectAxis",
    message = "Addition of 0.1 times sigma does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      ii = (envir$iter %% envir$n) + 1L
      ui = envir$e$vectors[, ii]
      if (any(envir$e$values[ii] < 0)) return(TRUE)
      lambdai = sqrt(envir$e$values[ii])
      m.old = envir$m.old
      return(sum((m.old - (m.old + 0.1 * envir$sigma * lambdai * ui))^2) < .Machine$double.eps)
    }
  ))
}

#' @title Stopping condition: standard deviation in coordinates.
#'
#' @description Stop if addition of 0.2 * standard deviations in any
#' coordinate does not change mean value.
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnNoEffectCoord = function() {
  return(makeStoppingCondition(
    name = "noEffectCoord",
    message = "Addition of 0.2 times sigma in any coordinate does not change mean value.",
    stop.fun = function(envir = parent.frame()) {
      m.old = envir$m.old
      return(sum((m.old - (m.old + 0.2 * envir$sigma))^2) < .Machine$double.eps)
    }
  ))
}

#' @title Stopping condition: high condition number.
#'
#' @description Stop if condition number of covariance matrix exceeds
#' tolerance value.
#'
#' @param tol [\code{numeric(1)}]\cr
#'   Tolerance value.
#'   Default is \code{1e14}.
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnCondCov = function(tol = 1e14) {
  assertNumber(tol, na.ok = FALSE, lower = 0, finite = TRUE)
  force(tol)
  return(makeStoppingCondition(
    name = "conditionCov",
    message = sprintf("Condition number of covariance matrix exceeds %f", tol),
    stop.fun = function(envir = parent.frame()) {
    #C = covmat
    #catch invalid values for covmat
    if (any(is.na(envir$C) | is.nan(envir$C) | is.infinite(envir$C))) return(TRUE)
    else return(kappa(envir$C) > tol)

    }
  ))
}


#===============================================================================
#=============================Online Convergence Detection======================
#===============================================================================
#' @title Stopping condition: Online Convergence Detection.
#'
#' @description Stop if OCD....
#'
#' @return [\code{cma_stopping_condition}]
#' @family stopping conditions
#' @export
stopOnOCD = function(varLimit, nPreGen,maxGen = NULL)
{
  # Check if varLimit is a single numeric
  assertNumber(varLimit, na.ok = FALSE)
  # Check if nPreGen is a single integerish value
  assertInt(nPreGen, na.ok = FALSE)
  # initialize significane level alpha with default value 0.05
  alpha = 0.05
  # Check if maxGen is a single integerish value
  if(!is.null(maxGen)) {
    assertInt(maxGen, na.ok = FALSE)
  }else{
    maxGen = Inf
  }
  # initialize p-values of Chi-squared variance test
  pvalue_current_gen_chi = numeric()
  pvalue_preceding_gen_chi = numeric()
  # initialize p-values of the t-test on the regression coefficient
  pvalue_current_gen_t = numeric()
  pvalue_preceding_gen_t = numeric()
  # return stopping condition being compatible with cma-es implementation by Jakob Bossek
  return(makeStoppingCondition(
    name = "OCD",
    message = sprintf("OCD successfully: Variance limit %f", varLimit),
    param.set = list(varLimit, nPreGen),
    stop.fun = function(envir = parent.frame()) {
      # Check if the number of iterations exceeds the user-defined number of maxGen. If TRUE, stop cma-es
      
      if(envir$iter >= maxGen){
        return(envir$iter >= maxGen)
      }
      
      # Check if number of iterations is greater than user-defined nPreGen
      if(envir$iter > nPreGen){
        # Here: In single objective optimization, the indicator of interest is the best fitness value of each generation.
        # PF_i is the best fitness value of the i-th generation which is used as a reference value
        # for calculating the indicator values of the last nPreGen generations
        PF_i = envir$best.fitness
        # PI_all is a vector with one entry for each generation, except the first generation.
        # PI_all stores the difference between the performance indicator values of the last nPreGen generations and the current generation i.
        PI_all = sapply((envir$generation.bestfitness)[-length(envir$generation.bestfitness)], function(x) x-PF_i, simplify = TRUE)
        # normalize PI_all in range upper.bound - lower.bound, i.e. the range of the objective values after nPreGen generations.
        # This value is fixed for all upcomming generations
        PI_all = PI_all/(envir$upper.bound-envir$lower.bound)
        # PI_current_gen is a subset of PI_all which stores the last nPreGen indicator values with respect to the current generation i.
        PI_current_gen = PI_all[(envir$iter-nPreGen):(envir$iter -1)]
        if((envir$iter - nPreGen) <= 1){
          # PI_preceding_gen is a subset of PI_all which stores the last nPreGen indicator values with respect to the last generation i-1.
          PI_preceding_gen = PI_current_gen
        }else{
          PI_preceding_gen =  PI_all[(envir$iter - (nPreGen+1)):(envir$iter - 2)]
        }
        # perform chi2 variance tests and return corresponding p-values
        pvalue_current_gen_chi = pChi2(varLimit, PI_current_gen)
        pvalue_preceding_gen_chi = pChi2(varLimit, PI_preceding_gen)
        # perform two-sided t-test and return corresponding p-values
        pvalue_current_gen_t = pReg(PI_current_gen)
        pvalue_preceding_gen_t = pReg(PI_preceding_gen)
        # log termination condition in cma_es
        if (pvalue_current_gen_chi <= alpha && pvalue_preceding_gen_chi <= alpha) envir$stopped.on.chi = envir$stopped.on.chi + 1
        if (pvalue_current_gen_t > alpha && pvalue_preceding_gen_t > alpha) envir$stopped.on.t = envir$stopped.on.t + 1
        # return TRUE, i.e. stop cmaes exectuion, if p-value is below specified significance level alpha
        return (pvalue_current_gen_chi <= alpha && pvalue_preceding_gen_chi <= alpha || pvalue_current_gen_t > alpha && pvalue_preceding_gen_t > alpha)
      }
      else{
        return(FALSE)
      }
    }
  ))
}

#' @export
pChi2 <- function (varLimit, PI) {
  # Determine degrees of freedom
  N = length(PI)-1
  # calculate test statistic
  Chi = (var(PI)*N)/varLimit
  # get p-value of corresponding chi2 test
  p = pchisq(Chi, N, lower.tail = TRUE)
  return (p)
}

pReg <- function (PI) {
  # Determin degrees of freedom
  N = length(PI)-1
  # standardize PI
  if(sd(PI)==0) {
    return (1)
  }else{
    PI = (PI-mean(PI))/sd(PI)
  }
  # initialize X, i.e. a vector of the generations of PI
  X = seq(1,length(PI),1)
  # standardize X
  X = (X-mean(X))/sd(X)
  # linear regression without intercept
  beta = (solve(X%*%X))%*%(X%*%PI)
  # residuals
  residuals = PI - X*beta
  # mse
  mse = (residuals%*%residuals) / N
  # compute test statistic
  t = beta/sqrt(mse*solve(X%*%X))
  # look up t distribution for N degrees of freedom
  p_value = 2*pt(-abs(t), df=N)
  return (p_value)
}

