stopOnOCD = function(varLimit, nPreGen, maxGen)
  {
  # Check if varLimit is a single numeric
  assertNumber(varLimit, na.ok = FALSE)
  # Check if nPreGen is a single integerish value
  assertInt(nPreGen, na.ok = FALSE)
  # initialize significane level alpha with default value 0.05
  alpha = 0.05
  # Check if maxGen is a single integerish value
  assertInt(maxGen, na.ok = FALSE)
  # return stopping condition being compatible with cma-es implementation by Jakob Bossek
  return(makeStoppingCondition(
    name = "Online Convergence Detection",
    message = sprintf("OCD successfully: Variance limit %f", varLimit),
    stop.fun = function(envir = parent.frame()) {
      # Check if the number of iterations exceeds the user-defined number of maxGen. If TRUE, stop cma-es
      if(envir$iter >= maxGen){
        return(envir$iter >= maxGen)
      }
      # Check if number of iterations is greater than user-defined nPreGen
        if(envir$iter > nPreGen){
          # PI_all is a vector with one entry for each generation. 
          # PI_all stores the difference between the performance indicator value of the preceding and the current generation.
          # Here: In single objective optimization, the performance indicator of interest is the best fitness value of each generation.
          PI_all = rbind(envir$generation.fitness[[1]], cbind(unlist(envir$generation.fitness))) - rbind(cbind(unlist(envir$generation.fitness)), envir$best.fitness)
          # PI_current_gen is a subset of PI_all which stores the last nPreGen indicator values with respect to the current generation i.
          PI_current_gen = PI_all[(envir$iter-nPreGen):envir$iter]
          if((envir$iter - nPreGen) > 1){
            # PI_preceding_gen is a subset of PI_all which stores the last nPreGen indicator values with respect to the last generation i-1.
            PI_preceding_gen =  PI_all[(envir$iter - (nPreGen+1)):(envir$iter - 1)]
          }else{
            PI_preceding_gen = PI_current_gen
          }
          # perform chi2 variance tests and return corresponding p-values
          pvalue_current_gen = pChi2(varLimit, PI_current_gen)
          pvalue_preceding_gen = pChi2(varLimit, PI_preceding_gen)
          # return TRUE, i.e. stop cmaes exectuion, if p-value is below specified significance level alpha
          return (pvalue_current_gen <= alpha && pvalue_preceding_gen <= alpha)
        }
        else{
          return(FALSE)
          }
       }
  ))
}

pChi2 <- function (varLimit, PI) {
  # Determine degrees of freedom
  N = length(PI)-1
  # calculate test statistic
  Chi = (var(PI)*N)/varLimit
  # get p-value of corresponding chi2 test
  p = pchisq(Chi, N, lower.tail = TRUE)
  return (p)
}

