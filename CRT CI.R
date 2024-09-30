# Confidence interval based on randomization tests applied to the adaptive two-stage enrichment analysis

library("data.table")

# source("CRTHelper.R")
source("helper.R")

# helper
# plot the probability of rejection as a function of the constant treatment effect in the Fisher's sharp null hypothesis.
# Example:
# rejection.probability.plot(data, tau.seq = seq(1, 11))
rejection.probability.plot = function(data, tau.seq,
                                      alpha = 0.05, tau = 0){
  matplot(tau.seq, data,
          type = "l", lwd = 3, lty = 1,
          col = c("black", "blue", "red"))
  abline(h = alpha, lwd = 3, lty = 3)
  abline(v = tau, lwd = 3, lty = 3)
  legend("topright",
         legend = c("RT","RT.second.stage", "CRT"),
         col = c("black", "blue", "red"),
         lwd = 3, cex = 0.5)
}

set.seed(318)
setting = "symmetric" # "default"; "large_n_stage_1"; "Rosenblum"; "heterogeneous_subgroup"; "monotone_counter_example", "symmetric"
N1 = 50; N2 = 25
tau = 2
alpha = 0.1
exact = F
alternative = "two-sided" # "two.sided", "less", "greater"
m.permutation = 100
count.max = 100 * m.permutation
return.all = T
eps = 0.001

inclusion.set = function(inclusion){
  if(inclusion == 0){return(numeric(0))}
  if(inclusion == 1){return(1)}
  if(inclusion == 2){return(0)}
  if(inclusion == 3){return(c(0,1))}
  stop("Please enter a valid index corresponding to the selected inclusion rule.")
}

inclusion.function = function(Y, A, X, ...){
  T0 = SATE(Y = Y[X == 0], A = A[X == 0])
  T1 = SATE(Y = Y[X == 1], A = A[X == 1])
  T.inclusion = (T1 - T0) / sqrt(2)

  if(is.na(T.inclusion)){
    return(T.inclusion)
  }

  if(T.inclusion > qnorm(0.8, 0, 1)){
    result = 1
  }else if(T.inclusion < qnorm(0.2, 0, 1)){
    result = 2
  }else{
    result = 3
  }

  return(result)
}

if(setting == "large_n_stage_1"){
  N1 = 100
}else if(setting == "Rosenblum"){
  inclusion.function = function(Y, A, X, ...){
    T0 = SATE(Y = Y[X == 0], A = A[X == 0])
    T1 = SATE(Y = Y[X == 1], A = A[X == 1])
    if(min(T1, T0) > 0){
      result = 3
    }else if(T1 > 0){
      result = 1
    }else if(T0 > 0){
      result = 2
    }else{
      result = 0
    }

    return(result)
  }
}else if(setting == "heterogeneous_subgroup"){
  tau = function(n, X){
    return(X + rep(-0.5, n))
  }
}else if(setting == "symmetric"){
  inclusion.function = function(Y, A, X, ...){
    T0 = SATE(Y = Y[X == 0], A = A[X == 0])
    T1 = SATE(Y = Y[X == 1], A = A[X == 1])
    if(abs(T0) > qnorm(0.8) & abs(T1) > qnorm(0.8)){
      result = 3
    }else if(abs(T1) > qnorm(0.8)){
      result = 1
    }else if(abs(T0) > qnorm(0.8)){
      result = 2
    }else{
      result = 0
    }

    return(result)
  }
}

m = 20
record = list()
record$CRT = record$RT.second.stage = record$RT = matrix(0, nrow = m, ncol = 1 + (alternative == "two-sided"))
record$inclusion = rep(0, m)

start.time = proc.time()
for(i in 1:m){
  while(record$inclusion[i] == 0){
    result = enrichment.design.data.gnr(N1 = N1, N2 = N2, n0 = N1/2, n1 = N1/2,
                                        tau = tau,
                                        inclusion.function = inclusion.function)
    record$inclusion[i] = result$inclusion
  }
  record$CRT[i,] = CI.enrichment.analysis(result = result,
                                          m.permutation = m.permutation,
                                         method = "CRT",
                                         sampler = "rejection.sampler",
                                         alpha = alpha,
                                         count.max = count.max,
                                         exact = exact,
                                         eps = eps,
                                         alternative = alternative) # "rejection.sampler"; "MCMC"

  record$RT.second.stage[i,] = CI.enrichment.analysis(result = result,
                                                      m.permutation = m.permutation,
                                                      method = "RT.second.stage",
                                                      alpha = alpha,
                                                      count.max = count.max,
                                                      exact = exact,
                                                      eps = eps,
                                                      alternative = alternative)

  record$RT[i,] = CI.enrichment.analysis(result = result,
                                        m.permutation = m.permutation,
                                        method = "RT",
                                        alpha = alpha,
                                        count.max = count.max,
                                        exact = exact,
                                        eps = eps,
                                        alternative = alternative)
  print(paste("Current iteration:", i))
}
end.time = proc.time()
print(end.time[3] - start.time[3])

par(mfrow = c(2,2))
table(record$inclusion)
for(region in list(c(1,2,3), 1, 2, 3)){
  if(alternative == "greater"){
    data = matrix(unlist(record[c("RT","RT.second.stage", "CRT")]),
                  ncol = 3)
    colnames(data) = c("RT","RT.2", "CRT")
    boxplot(data[record$inclusion %in% region,], ylab = "CI lower bound")
    print("coverage probability")
    print(apply((data[record$inclusion %in% region, ] < (tau + eps)), 2, mean))
  }else if(alternative == "less"){
    data = matrix(unlist(record[c("RT","RT.second.stage", "CRT")]),
                  ncol = 3)
    colnames(data) = c("RT","RT.2", "CRT")
    boxplot(data[record$inclusion %in% region,], ylab = "CI upper bound")
    print("coverage probability")
    print(apply((data[record$inclusion %in% region, ] > (tau - eps)), 2, mean))
  }else if(alternative == "two-sided"){
    data = matrix(unlist(record[c("RT","RT.second.stage", "CRT")]),
                  ncol = 3)
    colnames(data) = c("RT","RT.2", "CRT")
    width = data[-seq(1, dim(data)[1]/2),,drop = F] - data[seq(1, dim(data)[1]/2),,drop = F]
    colnames(width) = c("RT","RT.2", "CRT")
    boxplot(width[record$inclusion %in% region,], ylab = "CI width")
    print("coverage probability")
    print(apply((data[seq(1, dim(data)[1]/2),][record$inclusion %in% region, ] < (tau + eps))
                * (data[-seq(1, dim(data)[1]/2),][record$inclusion %in% region, ] > (tau - eps)), 2, mean))
  }
}

boxplot(abs(record$RT.second.stage - tau), ylim = c(0,4))

