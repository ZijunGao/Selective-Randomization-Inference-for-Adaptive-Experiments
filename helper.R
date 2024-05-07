# Helper functions for conditional randomization test for the enrichment analysis

# Data generator.
# Input:
# N1: the number of units enrolled in stage 1.
# n0: the number of units with covariate X = 0 in stage 1.
# n1: the number of units with covariate X = 1 in stage 1.
# N2: the number of units enrolled in stage 2.
# response.distribution: a function specifying the distribution of the control potential outcomes. Default is standard normal. Example: function(n){return(rnorm(n))}.
# tau: the constant treatment effect. Default is 0.
# treatment.design: "CRD" or "Bernoulli". For "CRD", a random half of the units are assigned to the treatment group. For "Bernoulli", each subject's treatment assignment indicator is i.i.d. Bernoulli(0.5). Default is "CRD".
# inclusion.function: a function specifying the rule of the second-stage subject inclusion based on the first-stage data. Example: function(...){return(3)}.
# inclusion.set: a function matching the output of the function "inclusion.function" to the subgroups to be included in stage 2. Example: function(inclusion){if(inclusion == 0){return(numeric(0))}; if(inclusion == 1){return(c(0,1))}}
# test.function: a function specifying the test statistic. Default is the mean function.
# Example:
# enrichment.design.data.gnr(N1 = 100, N2 = 50, n0 = 50, n1 = 50,
#                            inclusion.function = function(...){return(3)})
enrichment.design.data.gnr = function(N1, N2, n0 = N1/2, n1 = N1/2,
                                      response.distribution = function(n){return(rnorm(n))},
                                      tau = 0,
                                      treatment.design = "CRD",
                                      inclusion.function,
                                      inclusion.set = default.inclusion.set, 
                                      test.function = SATE){
  # stage 1
  X = c(rep(0, n0), rep(1, n1)) # two sub-groups
  Y.0 = response.distribution(N1) # control potential outcome
  A = treatment(N = N1, N1 = N1 / 2, treatment.design = treatment.design) # treatment assignment indicator
  while(min(c(sum(A * X), sum((1 - A) * X), sum(A * (1 - X)), sum((1 - A) * (1 - X)))) == 0){
    A = treatment(N = N1, N1 = N1 / 2, treatment.design = treatment.design) # repeat until each sub-group contains at least one treated unit and one control unit
  }
  data.stage.1 = data.gnr(n = N1, Y.0 = Y.0, A = A, tau = tau, X = X)
  # inclusion principle
  inclusion = inclusion.function(Y = data.stage.1$Y,
                                 A = data.stage.1$A,
                                 X = data.stage.1$X)
  hypothesis = inclusion.set(inclusion = inclusion)
  
  # stage 2
  Y.0 = response.distribution(N2)
  if(identical(hypothesis, numeric(0))){
    return(list("inclusion" = 0))
  }else if(identical(hypothesis, 1)){
    X = rep(1, N2)
  }else if(identical(hypothesis, 0)){
    X = rep(0, N2)
  }else if(identical(hypothesis, c(0, 1))){
    X = rep(0, N2)
    X[sample(N2, N2 * n1 / N1)] = 1
  }
  A = treatment(N = N2, N1 = N2 / 2, treatment.design = treatment.design)
  data.stage.2 = data.gnr(n = N2, Y.0 = Y.0, A = A, tau = tau, X = X)
  
  Y.final = c(data.stage.1$Y[which(data.stage.1$X %in% hypothesis)], data.stage.2$Y)
  A.final =  c(data.stage.1$A[which(data.stage.1$X %in% hypothesis)], data.stage.2$A)
  result = list("test.statistic" = test.function(Y = Y.final, A = A.final),
                "data.stage.1" = data.stage.1,
                "data.stage.2" = data.stage.2,
                "inclusion" = inclusion,
                "hypothesis" = hypothesis)
  return(result)
}


# Randomization tests for the adaptive two-stage enrichment analysis.
# Wrapper function.

# Input:
# result: a list of data produced by the function "enrichment.design.data.gnr".
# tau.seq: a sequence of constant treatment effects to test. Default is 0.
# method: type of the randomization test. Default is "CRT". Currently available options:
# "RT": randomization test.
# "RT.second.stage": randomization test only using the second-stage data.
# "CRT": conditional randomization test.
# sampler: Monte Carlo sampling method for "CRT". Available options:
# "rejections.sampler" (default): rejection sampling.
# "MCMC": Markov Chain Monte Carlo method, specifically the random scan, random group, random walk Metropolis-Hastings within Gibbs sampling algorithm.
# m.permutation: the number of permutations to be performed.
# return.all: if True, intermediate results are returned. See "Output" for more details.
# window.size: used by "Gibbs", "RWM", "RWMGibbs". "Gibbs", "RWMGibbs" divide all units into sub-groups of window.size. "RWM" randomly chooses a set of window.size treatments and permutes them.
# Output:
# p.val: p-value of the randomization test.
# return.all: here!!!
# Example:
# start.time = proc.time()
# set.seed(318)
# hist(replicate(1000,
#                inference.enrichment.analysis(result = enrichment.design.data.gnr(N1 = 40, N2 = 25, tau = 0, inclusion.function = inclusion.function),
#                                                           m.permutation = 100,
#                                                           tau.seq = 0,
#                                                           method = "CRT",
#                                                           sampler = "MH",
#                                                           count.max = 100 * m.permutation)), main = "", xlab = "p-value", breaks = 20, xlim = c(0,1), freq = F) # method = "CRT", "RT.second.stage", "RT"; "rejection.sampler", "MCMC"
# end.time = proc.time(); print(end.time[3] - start.time[3])
inference.enrichment.analysis = function(result,
                                         tau = 0,
                                         method = "CRT",
                                         exact = F,
                                         sampler = "rejection.sampler",
                                         design = "CRD",
                                         m.permutation = 100,
                                         count.max = 100 * m.permutation,
                                         return.all = F,
                                         alternative = "greater",
                                         window.size = 5,
                                         random.initialize = T,
                                         burn.in.size = 100){
  if(!exact){
    result = inference.enrichment.analysis.approx(result = result,
                                                  tau.seq = tau,
                                                  method = method,
                                                  sampler = sampler,
                                                  m.permutation = m.permutation,
                                                  design = design,
                                                  count.max = count.max,
                                                  return.all = return.all,
                                                  alternative = alternative,
                                                  window.size = window.size,
                                                  random.initialize = random.initialize,
                                                  burn.in.size = burn.in.size)
    
    return(result)
  }
  result = inference.enrichment.analysis.exact(result = result,
                                               tau.seq = tau,
                                               method = method,
                                               design = design,
                                               return.all = return.all,
                                               alternative = alternative)
  return(result)
}


# Randomization test applied to the adaptive two-stage enrichment analysis.
# Input:
# result: a list of data produced by the function "enrichment.design.data.gnr".
# tau.seq: a sequence of constant treatment effects. Default is 0.
# method: type of the randomization test. Default is "CRT". Available options:
# "RT": randomization test.
# "RT.second.stage": randomization test only using the second-stage data.
# "CRT" (default): conditional randomization test.
# sampler: sampling method for "CRT". Available options:
# "rejections.sampler" (default): rejection sampling.
# "MH" or "MH": Markov Chain Monte Carlo, Metropolis-Hastings.
# design: treatment assignment. Available options:
# "CRD" (default): completely randomized design.
# m.permutation: the number of permutations to be performed.
# return.all: if True, intermediate results are returned. See "Output" for more details.
# window.size: used by Gibbs sampling. Randomly choose a set of window.size treatments and randomly permute them.
# Output:
# p.val: a vector of p-values for each constant treatment effect specified in tau.seq.
# test.statistic: a matrix of test statistics: each row stands for a trial, and each column stands for a treatment effect specified in tau.seq. Only returned if return.all is set to be True.
# iteration: a vector of the number of completed permutations for each constant treatment effect specified in tau.seq. Only returned if return.all is set to be True.
# Example:
# See examples of inference.enrichment.analysis for details.

inference.enrichment.analysis.approx = function(result,
                                                tau.seq = 0,
                                                method = "CRT",
                                                design = "CRD",
                                                sampler = "rejection.sampler",
                                                m.permutation = 100,
                                                count.max = 100 * m.permutation,
                                                return.all = F,
                                                alternative = "greater", 
                                                window.size = 5,
                                                random.initialize = T,
                                                burn.in.size = 100){
  # preprocess
  N1 = length(result$data.stage.1$Y); N2 = length(result$data.stage.2$Y)
  inclusion = result$inclusion; test.statistic = result$test.statistic
  selected.stage.1 = which(result$data.stage.1$X %in% inclusion.set(inclusion)) # units in the first stage tested by the adaptive hypothesis
  selected.stage.2 = which(result$data.stage.2$X %in% inclusion.set(inclusion)) # units in the second stage tested by the adaptive hypothesis
  
  # impute the potential outcomes of the units tested by the adaptive hypothesis
  result$data.stage.1$Y.0.impute = result$data.stage.1$Y - matrix(tau.seq, byrow = T, nrow = N1, ncol = length(tau.seq)) * result$data.stage.1$A # each row stands for a unit; each column stands for a treatment effect level
  result$data.stage.1$Y.1.impute = result$data.stage.1$Y + matrix(tau.seq, byrow = T, nrow = N1, ncol = length(tau.seq)) * (1 - result$data.stage.1$A)
  result$data.stage.2$Y.0.impute = result$data.stage.2$Y - matrix(tau.seq, byrow = T, nrow = N2, ncol = length(tau.seq)) * result$data.stage.2$A
  result$data.stage.2$Y.1.impute = result$data.stage.2$Y + matrix(tau.seq, byrow = T, nrow = N2, ncol = length(tau.seq)) * (1 - result$data.stage.2$A)
  
  # randomization test
  if(method == "CRT"){
    if(sampler == "rejection.sampler"){
      test.statistic.ref = rejection.sampler(stage.1.A = result$data.stage.1$A,
                                             stage.2.A = result$data.stage.2$A,
                                             stage.1.X = result$data.stage.1$X,
                                             stage.2.X = result$data.stage.2$X,
                                             stage.1.Y.0.impute = result$data.stage.1$Y.0.impute,
                                             stage.1.Y.1.impute = result$data.stage.1$Y.1.impute,
                                             stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                                             stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                                             inclusion = inclusion,
                                             selected.stage.1 = selected.stage.1,
                                             selected.stage.2 = selected.stage.2,
                                             tau.seq = tau.seq,
                                             m.permutation = m.permutation,
                                             count.max = count.max)
    }else if(sampler %in% c("MCMC", "MH", "RWM")){
      result.full = RWM(stage.1.A = result$data.stage.1$A,
                        stage.2.A = result$data.stage.2$A,
                        stage.1.X = result$data.stage.1$X,
                        stage.2.X = result$data.stage.2$X,
                        stage.1.Y.0.impute = result$data.stage.1$Y.0.impute,
                        stage.1.Y.1.impute = result$data.stage.1$Y.1.impute,
                        stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                        stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                        inclusion = inclusion,
                        selected.stage.1 = selected.stage.1,
                        selected.stage.2 = selected.stage.2,
                        tau.seq = tau.seq,
                        m.permutation = m.permutation,
                        count.max = count.max,
                        window.size = window.size,
                        random.initialize = random.initialize,
                        burn.in.size = burn.in.size)
      test.statistic.ref = result.full$test.statistic.ref
    }
    if(sum(is.na(test.statistic.ref)) > 0){
      # use RT.second.stage as the backup when the CRT computation budget is exhausted
      method = "CRT to RT.second.stage"
      p.val.CRT = p.val.function(test.statistic = test.statistic,
                                 test.statistic.ref = test.statistic.ref,
                                 m.permutation = m.permutation,
                                 alternative = alternative)
    }
  }
  if(method %in% c("RT.second.stage", "CRT to RT.second.stage")){
    test.statistic = SATE(Y = result$data.stage.2$Y[selected.stage.2],
                          A = result$data.stage.2$A[selected.stage.2])
    if(is.na(test.statistic)){print("Please input a valid test statistic based on the original stage 2 data")}
    test.statistic.ref = RT.second.stage(stage.2.A = result$data.stage.2$A,
                                         stage.2.X = result$data.stage.2$X,
                                         stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                                         stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                                         selected.stage.2 = selected.stage.2,
                                         tau.seq = tau.seq,
                                         m.permutation = m.permutation,
                                         count.max = count.max)
  }
  if(method == "RT"){
    test.statistic.ref = RT(stage.1.A = result$data.stage.1$A,
                            stage.2.A = result$data.stage.2$A,
                            stage.1.X = result$data.stage.1$X,
                            stage.2.X = result$data.stage.2$X,
                            stage.1.Y.0.impute = result$data.stage.1$Y.0.impute,
                            stage.1.Y.1.impute = result$data.stage.1$Y.1.impute,
                            stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                            stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                            selected.stage.1 = selected.stage.1,
                            selected.stage.2 = selected.stage.2,
                            tau.seq = tau.seq,
                            m.permutation = m.permutation,
                            count.max = count.max)
  }
  p.val = p.val.function(test.statistic = test.statistic,
                         test.statistic.ref = test.statistic.ref,
                         m.permutation = m.permutation,
                         alternative = alternative)
  if(method == "CRT to RT.second.stage"){
    p.val[!is.na(p.val.CRT)] = p.val.CRT[!is.na(p.val.CRT)] # only use RT.second.stage for the hypotheses where the CRT reaches the maximum allowed number of randomizations
  }
  if(return.all){
    if(method != "CRT"){
      return(list(p.val = p.val, test.statistic = test.statistic.ref))
    }else{
      return(list(p.val = p.val, result.full = result.full))
    } 
  }
  return(p.val)
}



# Exact randomization test applied to the adaptive two-stage enrichment analysis.
# Input:
# result: a list of data produced by the function "enrichment.design.data.gnr".
# method: type of the randomization test. Default is "CRT". Available options:
# "RT": randomization test.
# "RT.second.stage": randomization test only using the second-stage data.
# "CRT" (default): conditional randomization test.
# exact: exact permutation test is adopted if True.
# design: treatment assignment. Available options:
# "CRD" (default): completely randomized design.
# return.all: if True, intermediate results are returned. See "Output" for more details.
# Output:
# p.val: a vector of p-values for each constant treatment effect specified in tau.seq.
# test.statistic: a matrix of test statistics: each row stands for a trial, and each column stands for a treatment effect specified in tau.seq. Only returned if return.all is set to be True.
# iteration: a vector of the number of acceptable treatment assignments for each constant treatment effect specified in tau.seq. Only returned if return.all is set to be True.
# Example:c
# See examples of inference.enrichment.analysis for details.
inference.enrichment.analysis.exact = function(result,
                                               tau.seq = 0,
                                               method = "CRT",
                                               design = "CRD",
                                               return.all = F,
                                               alternative = "greater"){
  # preprocess
  N1 = length(result$data.stage.1$Y); N2 = length(result$data.stage.2$Y)
  if(N1 > 10 | N2 > 10){stop("Please input N1, N2 <= 10")}
  inclusion = result$inclusion; test.statistic = result$test.statistic
  selected.stage.1 = which(result$data.stage.1$X %in% inclusion.set(inclusion)) # units in the first stage tested by the adaptive hypothesis
  selected.stage.2 = which(result$data.stage.2$X %in% inclusion.set(inclusion)) # units in the second stage tested by the adaptive hypothesis
  
  # impute the potential outcomes of the units tested by the adaptive hypothesis
  result$data.stage.1$Y.0.impute = result$data.stage.1$Y - matrix(tau.seq, byrow = T, nrow = N1, ncol = length(tau.seq)) * result$data.stage.1$A # each row stands for a unit; each column stands for a treatment effect level
  result$data.stage.1$Y.1.impute = result$data.stage.1$Y + matrix(tau.seq, byrow = T, nrow = N1, ncol = length(tau.seq)) * (1 - result$data.stage.1$A)
  result$data.stage.2$Y.0.impute = result$data.stage.2$Y - matrix(tau.seq, byrow = T, nrow = N2, ncol = length(tau.seq)) * result$data.stage.2$A
  result$data.stage.2$Y.1.impute = result$data.stage.2$Y + matrix(tau.seq, byrow = T, nrow = N2, ncol = length(tau.seq)) * (1 - result$data.stage.2$A)
  
  # randomization test
  if(method == "CRT"){
    test.statistic.ref = CRT.exact(stage.1.A = result$data.stage.1$A,
                                   stage.2.A = result$data.stage.2$A,
                                   stage.1.X = result$data.stage.1$X,
                                   stage.2.X = result$data.stage.2$X,
                                   stage.1.Y.0.impute = result$data.stage.1$Y.0.impute,
                                   stage.1.Y.1.impute = result$data.stage.1$Y.1.impute,
                                   stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                                   stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                                   inclusion = inclusion,
                                   selected.stage.1 = selected.stage.1,
                                   selected.stage.2 = selected.stage.2,
                                   tau.seq = tau.seq)
  }else if(method == "RT.second.stage"){
    test.statistic = SATE(Y = result$data.stage.2$Y[selected.stage.2],
                          A = result$data.stage.2$A[selected.stage.2])
    if(is.na(test.statistic)){
      print("Please input a valid test statistic based on the original stage 2 data")
    }
    test.statistic.ref =  RT.second.stage.exact(stage.2.A = result$data.stage.2$A,
                                                stage.2.X = result$data.stage.2$X,
                                                stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                                                stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                                                selected.stage.2 = selected.stage.2,
                                                tau.seq = tau.seq)
  }else if(method == "RT"){
    test.statistic.ref = RT.exact(stage.1.A = result$data.stage.1$A,
                                  stage.2.A = result$data.stage.2$A,
                                  stage.1.X = result$data.stage.1$X,
                                  stage.2.X = result$data.stage.2$X,
                                  stage.1.Y.0.impute = result$data.stage.1$Y.0.impute,
                                  stage.1.Y.1.impute = result$data.stage.1$Y.1.impute,
                                  stage.2.Y.0.impute = result$data.stage.2$Y.0.impute,
                                  stage.2.Y.1.impute = result$data.stage.2$Y.1.impute,
                                  selected.stage.1 = selected.stage.1,
                                  selected.stage.2 = selected.stage.2,
                                  tau.seq = tau.seq)
  }
  p.val = p.val.function(test.statistic = test.statistic,
                         test.statistic.ref = test.statistic.ref,
                         alternative = alternative,
                         exact = T)
  if(return.all){return(list(p.val = p.val, test.statistic = test.statistic.ref, iteration = apply(!is.na(test.statistic.ref), 3, sum)))}
  return(p.val)
}


# Confidence interval based on randomization test applied to the adaptive two-stage enrichment analysis.
# Input:
# result: a list of data produced by the function "enrichment.design.data.gnr".
# tau.seq: a sequence of constant treatment effects. Default is 0.
# method: type of the randomization test. Default is "CRT". Available options:
# "RT": randomization test.
# "RT.second.stage": randomization test only using the second-stage data.
# "CRT" (default): conditional randomization test.
# sampler: sampling method for "CRT". Available options:
# "rejections.sampler" (default): rejection sampling.
# "MCMC": Markov Chain Monte Carlo method, specifically the random scan, random group, random walk Metropolis Hastings within Gibbs algorithm.   
# design: treatment assignment. Available options:
# "CRD" (default): completely randomized design.
# m.permutation: the number of permutations to be performed.
# count.max: the maximal number of candidate treatment assignments generated.
# exact: exact permutation test is adopted if True.
# alpha: significance level of the randomization test.
# left: lower bound of the constant treatment effect for binary search.
# right: upper bound of the constant treatment effect for binary search.
# eps: resolution of the binary search, the binary search is stopped if the search interval is narrower than eps.
# window.size: used by Gibbs sampling. Randomly choose a set of window.size treatments and randomly permute them. Default is 5.
# Output:
# a number indicating the lower bound of the one-sided interval [c, \inf).
# Example:
# result = enrichment.design.data.gnr(N1 = 50, N2 = 50, tau = 0, inclusion.function = function(x,...){return(1)})
# CI.enrichment.analysis(result = result, sampler = "MH")
CI.enrichment.analysis = function(result,
                                  method = "CRT",
                                  sampler = "rejection.sampler",
                                  design = "CRD",
                                  m.permutation = 100,
                                  count.max = 100 * m.permutation,
                                  exact = F,
                                  alpha = 0.05,
                                  left = -10,
                                  right = 10,
                                  eps = 0.1,
                                  alternative = "greater",
                                  window.size = 5){
  if(alternative %in% c("greater", "less")){
    CI = CI.enrichment.analysis.helper(result = result,
                                       method = method,
                                       sampler = sampler,
                                       design = design,
                                       m.permutation = m.permutation,
                                       count.max = count.max,
                                       exact = exact,
                                       alpha = alpha,
                                       left = left,
                                       right = right,
                                       eps = eps,
                                       alternative = alternative,
                                       window.size = window.size)
    return(CI)
  }
  CI.lower = CI.enrichment.analysis.helper(result = result,
                                           method = method,
                                           sampler = sampler,
                                           design = design,
                                           m.permutation = m.permutation,
                                           count.max = count.max,
                                           exact = exact,
                                           alpha = alpha/2,
                                           left = left,
                                           right = right,
                                           eps = eps,
                                           alternative = "greater",
                                           window.size = window.size)
  CI.upper = CI.enrichment.analysis.helper(result = result,
                                           method = method,
                                           sampler = sampler,
                                           design = design,
                                           m.permutation = m.permutation,
                                           count.max = count.max,
                                           exact = exact,
                                           alpha = alpha/2,
                                           left = left,
                                           right = right,
                                           eps = eps,
                                           alternative = "less",
                                           window.size = window.size)
  return(c(CI.lower, CI.upper))
}


CI.enrichment.analysis.helper = function(result,
                                         method = "CRT",
                                         sampler = "rejection.sampler",
                                         design = "CRD",
                                         m.permutation = 100,
                                         count.max = 100 * m.permutation,
                                         exact = F,
                                         alpha = 0.05,
                                         left = -10,
                                         right = 10,
                                         eps = 0.1,
                                         alternative = "greater",
                                         window.size = window.size){
  # binary search
  while (TRUE) {
    p.val = inference.enrichment.analysis(result = result,
                                          tau = left,
                                          method = method,
                                          design = design,
                                          sampler = sampler,
                                          m.permutation = m.permutation,
                                          count.max = count.max,
                                          exact = exact,
                                          alternative = alternative,
                                          window.size = window.size)
    
    if (p.val >= alpha & left > -100 & alternative == "greater") {
      left = left - 10
    }else if (p.val < alpha & left > -100 & alternative == "less") {
      left = left - 10
    }else {
      break
    }
  }
  
  while (TRUE) {
    p.val = inference.enrichment.analysis(result = result,
                                          tau = right,
                                          method = method,
                                          design = design,
                                          sampler = sampler,
                                          m.permutation = m.permutation,
                                          count.max = count.max,
                                          exact = exact,
                                          alternative = alternative,
                                          window.size = window.size)
    if (p.val < alpha & right < 100 & alternative == "greater"){
      right = right + 10
    }else if (p.val >= alpha & right < 100 & alternative == "less"){
      right = right + 10
    }else{
      break
    }
  }
  
  while (left <= (right - eps)) {
    mid = (left + right) / 2
    p.val = inference.enrichment.analysis(result = result,
                                          tau = mid,
                                          method = method,
                                          design = design,
                                          sampler = sampler,
                                          m.permutation = m.permutation,
                                          count.max = count.max,
                                          exact = exact,
                                          alternative = alternative,
                                          window.size = window.size)
    
    if ((p.val < alpha & alternative == "greater") | (p.val >= alpha & alternative == "less")){
      left = mid
    } else {
      right = mid
    }
  }
  return(mid)
}

# TODO: use data from both stages to determine the hypothesis

# TODO: add prognostic covariates; influence the test statistics, inclusion principle



# Standardized average treatment effect
# Input:
# Y: response vector.
# A: treatment assignment vector.
# Example:
# hist(replicate(400, SATE(Y = rnorm(400,0,1), A = rbinom(400,1,0.5))))
SATE = function(Y, A, ...){
  # no treated units
  if(sum(A) == 0){
    return(-mean(Y) / (sd(Y) / length(Y)))
  }
  # no control units
  if(sum(1-A) == 0){
    return(mean(Y) / (sd(Y) / length(Y)))
  }
  diff.in.means = mean(Y[which(A == 1)]) - mean(Y[which(A == 0)])
  # estimator of the standard deviation of the difference in means
  var.0 = var(Y[which(A == 0)])
  var.1 = var(Y[which(A == 1)])
  sd.estimate = sqrt(var.1 / sum(A) + var.0 / sum(1 - A))
  # if there is only one treated unit or one control unit, we estimate sd(Y(1)) and sd(Y(0)) by sd(Y)
  if(is.na(sd.estimate)){
    sd.estimate = sd(Y) * sqrt(1 / sum(A) + 1 / sum(1 - A))
  }
  result = diff.in.means / sd.estimate
  return(result)
}


# Helper: data generator.
# Input:
# n: sample size.
# Y.0: control potential outcome.
# tau: constant treatment effect.
# X: covariate vector/matrix.
# Example:
# data.gnr(n = 10, Y.0 = rep(0, 10), A = rbinom(10, 1, 0.5), tau = 1)
data.gnr = function(n, Y.0, A = rep(0, n), tau = 0, X = NULL){
  if(is.function(tau)){tau = tau(n, X)}
  Y.1 = Y.0 + tau
  Y = Y.0
  Y[which(A == 1)] = Y.1[which(A == 1)]
  result = list(X = X, A = A, Y = Y, Y.0 = Y.0, Y.1 = Y.1)
  return(result)
}


# Generate treatment indicators
# Input:
# N: number of units.
# N1: number of treated units for treatment.design = "CRD".
# prob: probability of receiving the treatment for treatment.design = "Bernoulli".
# X: binary subgroup indicator.
# treatment.design: treatment design. Acceptable values: "CRD", "Bernoulli".
# Example:
# hist(replicate(400, SATE(Y = rnorm(400,0,1), A = rbinom(400,1,0.5))))
treatment = function(N, N1 = N/2, prob = 0.5, treatment.design = "CRD"){
  A = rep(0, N)
  if(treatment.design == "CRD"){
    A[sample(N, N1)] = 1
  }else if(treatment.design == "Bernoulli"){
    A = rbinom(N, 1, prob)
  }
  return(A)
}

default.inclusion.set = function(inclusion){
  if(inclusion == 0){return(numeric(0))}
  if(inclusion == 1){return(1)}
  if(inclusion == 2){return(0)}
  if(inclusion == 3){return(c(0,1))}
  stop("Please enter a valid index corresponding to the selected inclusion rule.")
}

inclusion.set = default.inclusion.set



rejection.sampler = function(stage.1.A,
                             stage.2.A,
                             stage.1.X,
                             stage.2.X,
                             stage.1.Y.0.impute,
                             stage.1.Y.1.impute,
                             stage.2.Y.0.impute,
                             stage.2.Y.1.impute,
                             inclusion,
                             selected.stage.1,
                             selected.stage.2,
                             tau.seq,
                             m.permutation = 100,
                             count.max = count.max){
  stage.1.A.ref = stage.1.A
  stage.2.A.ref = stage.2.A
  stage.1.Y.ref = stage.1.Y.1.impute * stage.1.A.ref + stage.1.Y.0.impute * (1 - stage.1.A.ref)
  stage.2.Y.ref = stage.2.Y.1.impute * stage.2.A.ref + stage.2.Y.0.impute * (1 - stage.2.A.ref)
  count = rep(0, length(tau.seq)); count.total = 0
  test.statistic.ref = matrix(NA, nrow = m.permutation, ncol = length(tau.seq))
  while(min(count) < m.permutation){
    if(count.total > count.max){
      # print("Rejection sampling: insufficient acceptable stage 1 treatment assignments for some sharp null hypothesis. Change from CRT with rejection sampling to RT.second.stage")
      method = "CRT to RT.second.stage"
      break
    }
    count.total = count.total + 1
    # regenerate treatment assignments using CRD
    stage.1.A.ref[selected.stage.1] = sample(stage.1.A[selected.stage.1])
    stage.1.Y.ref[which(stage.1.A.ref == 0),] = stage.1.Y.0.impute[which(stage.1.A.ref == 0),]
    stage.1.Y.ref[which(stage.1.A.ref == 1),] = stage.1.Y.1.impute[which(stage.1.A.ref == 1),]
    inclusion.ref = apply(stage.1.Y.ref, 2, inclusion.function, stage.1.A.ref, stage.1.X) # inclusion rule based on the regenerated treatment assignments
    if(max(is.na(inclusion.ref)) > 0){
      stop("The inclusion rule is not well-defined under some sharp null hypothesis")
    }
    valid.assignment.index = which(inclusion.ref == inclusion) # the sharp null hypotheses under which the regenerated treatment assignment is accepted
    if(length(valid.assignment.index) == 0){next}
    count[valid.assignment.index] = pmin(m.permutation, count[valid.assignment.index] + 1)
    stage.2.A.ref[selected.stage.2] = sample(stage.2.A[selected.stage.2])
    stage.2.Y.ref[which(stage.2.A.ref == 0),] = stage.2.Y.0.impute[which(stage.2.A.ref == 0),]
    stage.2.Y.ref[which(stage.2.A.ref == 1),] = stage.2.Y.1.impute[which(stage.2.A.ref == 1),]
    test.statistic.ref[cbind(count[valid.assignment.index], valid.assignment.index)] = sapply(valid.assignment.index, function(j){SATE(Y = c(stage.1.Y.ref[selected.stage.1, j], stage.2.Y.ref[selected.stage.2, j]), A = c(stage.1.A.ref[selected.stage.1], stage.2.A.ref[selected.stage.2]))})
    if(max(is.na(test.statistic.ref[cbind(count[valid.assignment.index], valid.assignment.index)])) > 0){stop("The test statistic is not well-defined under some sharp null hypothesis")}
  }
  return(test.statistic.ref)
}


RWM = function(stage.1.A,
               stage.2.A,
               stage.1.X,
               stage.2.X,
               stage.1.Y.0.impute,
               stage.1.Y.1.impute,
               stage.2.Y.0.impute,
               stage.2.Y.1.impute,
               inclusion,
               selected.stage.1,
               selected.stage.2,
               tau.seq,
               m.permutation = 100,
               count.max = 100 * m.permutation,
               window.size = 2,
               burn.in.size = 100,
               random.initialize = T){
  # random initialization
  stage.1.A.ref = matrix(stage.1.A, nrow = length(stage.1.A), ncol = length(tau.seq))
  stage.1.Y.ref = stage.1.Y.1.impute * stage.1.A.ref + stage.1.Y.0.impute * (1 - stage.1.A.ref)
  stage.2.A.ref = stage.2.A
  stage.2.Y.ref = stage.2.Y.1.impute * stage.2.A.ref + stage.2.Y.0.impute * (1 - stage.2.A.ref)
  if(random.initialize){
    inclusion.ref = rep(NA, length(tau.seq))
    initialized = rep(0, length(tau.seq)) # the indicator of whether an acceptable treatment has been found/the initialization has been completed for a null hypothesis
    for(i in 1 : floor(count.max / m.permutation)){
      stage.1.A.ref[selected.stage.1, which(initialized == 0)] = sample(stage.1.A[selected.stage.1]) # permute units
      stage.1.Y.ref[, which(initialized == 0)] = stage.1.Y.1.impute[, which(initialized == 0)] * stage.1.A.ref[, which(initialized == 0)] + stage.1.Y.0.impute[, which(initialized == 0)] * (1 - stage.1.A.ref[, which(initialized == 0)])
      inclusion.ref[which(initialized == 0)] = sapply(which(initialized == 0), function(j){return(inclusion.function(Y = stage.1.Y.ref[, j], A = stage.1.A.ref[, j], X = stage.1.X))}) 
      if(max(is.na(inclusion.ref)) > 0){stop("The inclusion rule is not well-defined under some sharp null hypothesis")}
      initialized[which(inclusion.ref == inclusion)] = 1
      if(min(initialized) == 1){break}
    }
    if(min(initialized) == 0){
      # when a null hypothesis fails to find an appropriate initialization point, it is initialized at the observed treatment assignment
      stage.1.A.ref[selected.stage.1, which(initialized == 0)] = stage.1.A[selected.stage.1] 
      stage.1.Y.ref[, which(initialized == 0)] = stage.1.Y.1.impute[, which(initialized == 0)] * stage.1.A.ref[, which(initialized == 0)] + stage.1.Y.0.impute[, which(initialized == 0)] * (1 - stage.1.A.ref[, which(initialized == 0)])
    }
  }
  
  # MCMC
  # random scan, random group, random walk Metropolis-Hastings within Gibbs sampling algorithm
  m.permutation = m.permutation + burn.in.size
  count = rep(0, length(tau.seq)); count.total = 0
  test.statistic.ref = matrix(NA, nrow = m.permutation, ncol = length(tau.seq))
  accept.indicator = jump.distance.squared = matrix(0, nrow = m.permutation, ncol = length(tau.seq)) # record the number of accepted proposals, the squared Euclidean jump distance, per hypothesis
  group.number = floor(length(selected.stage.1) / window.size)
  start.group = (seq(1, group.number) - 1) *  window.size + 1
  end.group = seq(1, group.number) *  window.size; end.group[group.number] = length(selected.stage.1)
  permuted.index = sample(length(selected.stage.1)) # randomly permute the selected units from stage 1 in case the units are pre-ordered
  while(min(count) < m.permutation){
    count.total = count.total + 1
    # regenerate treatment assignments
    swap.group.index = sample(group.number, 1)
    # swap.index = sample(selected.stage.1[permuted.index][seq(start.group[swap.group.index], end.group[swap.group.index])], window.size, replace = F) # randomly select a group and permute
    swap.index = sample(selected.stage.1, window.size, replace = F) # randomly select window.size units and permute
    swap.index.original = sort(swap.index)
    stage.1.A.ref[swap.index.original, ] = stage.1.A.ref[swap.index, ]
    stage.1.Y.ref = stage.1.Y.1.impute * stage.1.A.ref + stage.1.Y.0.impute * (1 - stage.1.A.ref)
    inclusion.ref = sapply(seq(1, length(tau.seq)), function(j){return(inclusion.function(Y = stage.1.Y.ref[, j], A = stage.1.A.ref[, j], X = stage.1.X))}) 
    if(max(is.na(inclusion.ref)) > 0){stop("The inclusion rule is not well-defined under some sharp null hypothesis")}
    valid.assignment.index = which(inclusion.ref == inclusion)
    
    accept.indicator[count.total, valid.assignment.index] = 1 # record the number of accepted proposals
    jump.distance.squared[count.total, valid.assignment.index] = apply(stage.1.A.ref[swap.index.original, valid.assignment.index, drop = F] - stage.1.A.ref[swap.index.original[order(swap.index)], valid.assignment.index, drop = F], 2, function(x){return(mean(x^2))})
    count = pmin(m.permutation, count + 1)
    if(length(valid.assignment.index) < length(tau.seq)){
      # backtracking: when the new proposal is rejected, the Metropolis-Hastings algorithm will take the original sample as the new sample
      backtrack.index = setdiff(seq(1, length(tau.seq)), valid.assignment.index)
      stage.1.A.ref[swap.index.original, backtrack.index] = stage.1.A.ref[swap.index.original[order(swap.index)], backtrack.index]
      stage.1.Y.ref = stage.1.Y.1.impute * stage.1.A.ref + stage.1.Y.0.impute * (1 - stage.1.A.ref)
    }
    
    # second stage
    stage.2.A.ref[selected.stage.2] = sample(stage.2.A[selected.stage.2])
    stage.2.Y.ref[which(stage.2.A.ref == 0), ] = stage.2.Y.0.impute[which(stage.2.A.ref == 0), ]
    stage.2.Y.ref[which(stage.2.A.ref == 1), ] = stage.2.Y.1.impute[which(stage.2.A.ref == 1), ]
    test.statistic.ref[cbind(count, seq(1, length(tau.seq)))] = sapply(seq(1, length(tau.seq)), function(j){SATE(Y = c(stage.1.Y.ref[selected.stage.1, j], stage.2.Y.ref[selected.stage.2, j]), A = c(stage.1.A.ref[selected.stage.1, j], stage.2.A.ref[selected.stage.2]))})
    if(max(is.na(test.statistic.ref[cbind(count[seq(1, length(tau.seq))], seq(1, length(tau.seq)))])) > 0){stop("The test statistic is not well-defined under some sharp null hypothesis")}
  }
  return(list(test.statistic.ref = test.statistic.ref[-seq(1, burn.in.size), ], acceptance.rate = apply(accept.indicator[seq(burn.in.size + 1, m.permutation), ,drop = F], 2, mean), jump.distance.squared = apply(jump.distance.squared[seq(burn.in.size + 1, m.permutation),,drop = F], 2, mean)))
}



RT.second.stage = function(stage.2.A,
                           stage.2.X,
                           stage.2.Y.0.impute,
                           stage.2.Y.1.impute,
                           selected.stage.2,
                           tau.seq,
                           m.permutation = 100,
                           count.max = 100 * m.permutation){
  test.statistic.ref = matrix(NA, nrow = m.permutation, ncol = length(tau.seq))
  stage.2.A.ref = stage.2.A
  stage.2.Y.ref = stage.2.Y.0.impute
  for(i in 1 : m.permutation){
    stage.2.A.ref[selected.stage.2] = sample(stage.2.A[selected.stage.2])
    stage.2.Y.ref[which(stage.2.A.ref == 0),] = stage.2.Y.0.impute[which(stage.2.A.ref == 0),]
    stage.2.Y.ref[which(stage.2.A.ref == 1),] = stage.2.Y.1.impute[which(stage.2.A.ref == 1),]
    test.statistic.ref[i,] = sapply(seq(1, length(tau.seq)), function(j){return(SATE(Y = stage.2.Y.ref[selected.stage.2, j], A = stage.2.A.ref[selected.stage.2]))})
    if(max(is.na(test.statistic.ref[i,])) > 0){stop("The test statistic is not well-defined for some sharp null")}
  }
  return(test.statistic.ref)
}


RT = function(stage.1.A,
              stage.2.A,
              stage.1.X,
              stage.2.X,
              stage.1.Y.0.impute,
              stage.1.Y.1.impute,
              stage.2.Y.0.impute,
              stage.2.Y.1.impute,
              selected.stage.1,
              selected.stage.2,
              tau.seq,
              m.permutation = 100,
              count.max = 100 * m.permutation){
  stage.1.Y.ref = stage.1.Y.0.impute
  stage.2.Y.ref = stage.2.Y.0.impute
  test.statistic.ref = matrix(NA, nrow = m.permutation, ncol = length(tau.seq))
  for(i in 1 : m.permutation){
    stage.1.A.ref = sample(stage.1.A)
    stage.1.Y.ref[which(stage.1.A.ref == 1),] = stage.1.Y.1.impute[which(stage.1.A.ref == 1),]
    stage.1.Y.ref[which(stage.1.A.ref == 0),] = stage.1.Y.0.impute[which(stage.1.A.ref == 0),]
    stage.2.A.ref = sample(stage.2.A)
    stage.2.Y.ref[which(stage.2.A.ref == 1),] = stage.2.Y.1.impute[which(stage.2.A.ref == 1),]
    stage.2.Y.ref[which(stage.2.A.ref == 0),] = stage.2.Y.0.impute[which(stage.2.A.ref == 0),]
    test.statistic.ref[i,] = sapply(seq(1, length(tau.seq)), function(j){SATE(Y = c(stage.1.Y.ref[selected.stage.1, j], stage.2.Y.ref[selected.stage.2, j]), A = c(stage.1.A.ref[selected.stage.1], stage.2.A.ref[selected.stage.2]))})
    if(max(is.na(test.statistic.ref[i,])) > 0){stop("The test statistic is not well-defined for some sharp null")}
  }
  return(test.statistic.ref)
}



CRT.exact = function(stage.1.A,
                     stage.2.A,
                     stage.1.X,
                     stage.2.X,
                     stage.1.Y.0.impute,
                     stage.1.Y.1.impute,
                     stage.2.Y.0.impute,
                     stage.2.Y.1.impute,
                     inclusion,
                     selected.stage.1,
                     selected.stage.2,
                     tau.seq){
  stage.1.A.ref = stage.1.A
  stage.2.A.ref = stage.2.A
  stage.1.Y.ref = stage.1.Y.0.impute
  stage.2.Y.ref = stage.2.Y.0.impute
  # enumerate over all possible permutations
  stage.1.A.All = gtools::combinations(n = length(selected.stage.1), r = sum(stage.1.A[selected.stage.1]),  v = selected.stage.1)
  stage.2.A.All = gtools::combinations(n = length(selected.stage.2), r = sum(stage.2.A[selected.stage.2]),  v = selected.stage.2)
  test.statistic.ref = array(NA, dim = c(dim(stage.1.A.All)[1], dim(stage.2.A.All)[1], length(tau.seq)))
  for(i in 1 : dim(stage.1.A.All)[1]){
    stage.1.A.ref[stage.1.A.All[i, ]] = 1
    stage.1.A.ref[setdiff(selected.stage.1, stage.1.A.All[i,])] = 0
    stage.1.Y.ref[which(stage.1.A.ref == 1), ] = stage.1.Y.1.impute[which(stage.1.A.ref == 1), ]
    stage.1.Y.ref[which(stage.1.A.ref == 0), ] = stage.1.Y.0.impute[which(stage.1.A.ref == 0), ]
    inclusion.ref = sapply(seq(1,length(tau.seq)), function(j){inclusion.function(Y = stage.1.Y.ref[, j], A = stage.1.A.ref, X = stage.1.X)})
    valid.assignment.index = which(inclusion.ref == inclusion)
    if(length(valid.assignment.index) == 0){next}
    for(j in 1 : dim(stage.2.A.All)[1]){
      stage.2.A.ref[stage.2.A.All[j, ]] = 1
      stage.2.A.ref[setdiff(selected.stage.2, stage.2.A.All[j, ])] = 0
      stage.2.Y.ref[which(stage.2.A.ref == 1), ] = stage.1.Y.1.impute[which(stage.2.A.ref == 1), ]
      stage.2.Y.ref[which(stage.2.A.ref == 0), ] = stage.1.Y.0.impute[which(stage.2.A.ref == 0), ]
      test.statistic.ref[i, j, valid.assignment.index] = sapply(valid.assignment.index, function(l){SATE(Y = c(stage.1.Y.ref[selected.stage.1, l], stage.2.Y.ref[selected.stage.2, l]), A = c(stage.1.A.ref[selected.stage.1], stage.2.A.ref[selected.stage.2]))})
    }
  }
  return(test.statistic.ref)
}


RT.second.stage.exact = function(stage.2.A,
                                 stage.2.X,
                                 stage.2.Y.0.impute,
                                 stage.2.Y.1.impute,
                                 selected.stage.2,
                                 tau.seq){
  stage.2.A.ref = stage.2.A
  stage.2.Y.ref = stage.2.Y.0.impute
  # enumerate over all possible permutations
  stage.2.A.All = gtools::combinations(n = length(selected.stage.2), r = sum(stage.2.A[selected.stage.2]),  v = selected.stage.2)
  test.statistic.ref = array(NA, dim = c(1, dim(stage.2.A.All)[1], length(tau.seq)))
  for(j in 1 : dim(stage.2.A.All)[1]){
    stage.2.A.ref[stage.2.A.All[j, ]] = 1
    stage.2.A.ref[setdiff(selected.stage.2, stage.2.A.All[j, ])] = 0
    stage.2.Y.ref[which(stage.2.A.ref == 1), ] = stage.2.Y.1.impute[which(stage.2.A.ref == 1), ]
    stage.2.Y.ref[which(stage.2.A.ref == 0), ] = stage.2.Y.0.impute[which(stage.2.A.ref == 0), ]
    test.statistic.ref[1, j, ] = sapply(seq(1, length(tau.seq)), function(l){return(SATE(Y = stage.2.Y.ref[selected.stage.2, l], A = stage.2.A.ref[selected.stage.2]))})
    if(max(is.na(test.statistic.ref[1, j, ])) > 0){stop("The test statistic is not well-defined for some sharp null")}
  }
  return(test.statistic.ref)
}


RT.exact = function(stage.1.A,
                    stage.2.A,
                    stage.1.X,
                    stage.2.X,
                    stage.1.Y.0.impute,
                    stage.1.Y.1.impute,
                    stage.2.Y.0.impute,
                    stage.2.Y.1.impute,
                    selected.stage.1,
                    selected.stage.2,
                    tau.seq){
  stage.1.A.ref = stage.1.A
  stage.2.A.ref = stage.2.A
  stage.1.Y.ref = stage.1.Y.0.impute
  stage.2.Y.ref = stage.2.Y.0.impute
  # enumerate over all possible permutations
  stage.1.A.All = gtools::combinations(n = length(selected.stage.1), r = sum(stage.1.A[selected.stage.1]),  v = selected.stage.1)
  stage.2.A.All = gtools::combinations(n = length(selected.stage.2), r = sum(stage.2.A[selected.stage.2]),  v = selected.stage.2)
  test.statistic.ref = array(NA, dim = c(dim(stage.1.A.All)[1], dim(stage.2.A.All)[1], length(tau.seq)))
  for(i in 1 : dim(stage.1.A.All)[1]){
    stage.1.A.ref[stage.1.A.All[i, ]] = 1
    stage.1.A.ref[setdiff(selected.stage.1, stage.1.A.All[i,])] = 0
    stage.1.Y.ref[which(stage.1.A.ref == 1), ] = stage.1.Y.1.impute[which(stage.1.A.ref == 1), ]
    stage.1.Y.ref[which(stage.1.A.ref == 0), ] = stage.1.Y.0.impute[which(stage.1.A.ref == 0), ]
    for(j in 1 : dim(stage.2.A.All)[1]){
      stage.2.A.ref[stage.2.A.All[j, ]] = 1
      stage.2.A.ref[setdiff(selected.stage.2, stage.2.A.All[j, ])] = 0
      stage.2.Y.ref[which(stage.2.A.ref == 1), ] = stage.2.Y.1.impute[which(stage.2.A.ref == 1), ]
      stage.2.Y.ref[which(stage.2.A.ref == 0), ] = stage.2.Y.0.impute[which(stage.2.A.ref == 0), ]
      test.statistic.ref[i, j, ] = sapply(seq(1, length(tau.seq)), function(l){SATE(Y = c(stage.1.Y.ref[selected.stage.1, l], stage.2.Y.ref[selected.stage.2, l]), A = c(stage.1.A.ref[selected.stage.1], stage.2.A.ref[selected.stage.2]))})
    }
  }
  return(test.statistic.ref)
}


# Example:
# p.val.function(test.statistic = result$test.statistic,
#                test.statistic.ref = inference.CRT$test.statistic,
#                m.permutation,
#                alternative = "greater")
p.val.function = function(test.statistic, test.statistic.ref, m.permutation = 100, alternative = "greater", exact = F){
  if(!exact){
    if(alternative == "greater"){
      p.val = (1 + apply(as.matrix(test.statistic <= test.statistic.ref), 2, sum))/(1 + m.permutation)
    }else if(alternative == "less"){
      p.val = (1 + apply(as.matrix(test.statistic >= test.statistic.ref), 2, sum))/(1 + m.permutation)
    }else if(alternative == "two-sided"){
      p.val.greater = (1 + apply(as.matrix(test.statistic <= test.statistic.ref), 2, sum))/(1 + m.permutation)
      p.val.less = (1 + apply(as.matrix(test.statistic >= test.statistic.ref), 2, sum))/(1 + m.permutation)
      p.val = 2 * pmin(p.val.greater, p.val.less)
    }
    return(p.val)
  }
  # TODO: check the dimension of (test.statistic <= test.statistic.ref)
  if(alternative == "greater"){
    p.val = (1 + apply((test.statistic <= test.statistic.ref), 3, sum, na.rm = T))/(1 + apply(!is.na(test.statistic.ref), 3, sum))
  }else if(alternative == "less"){
    p.val = (1 + apply((test.statistic >= test.statistic.ref), 3, sum, na.rm = T))/(1 + apply(!is.na(test.statistic.ref), 3, sum))
  }else if(alternative == "two-sided"){
    p.val.greater = (1 + apply((test.statistic <= test.statistic.ref), 3, sum, na.rm = T))/(1 + apply(!is.na(test.statistic.ref), 3, sum))
    p.val.less = (1 + apply((test.statistic >= test.statistic.ref), 3, sum, na.rm = T))/(1 + apply(!is.na(test.statistic.ref), 3, sum))
    p.val = 2 * pmin(p.val.greater, p.val.less)
  }
  return(p.val)
}
