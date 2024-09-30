# Create a hypothetical two-stage enrichment trial based on the SPRINT trial.
# Subgroup: age groups: [0, 59], [60, 69], [70, 79], [80, 100].

# helper function
SATE = function(Y, A, ...){
  result = log(mean(Y[A == 1])) - log(mean(Y[A == 0]))
  result = -result # the larger the result, the larger the treatment effect
  return(result)
}

inclusion.function = function(Y, A, X, ...){
  S = rep(0, 4) # four age groups
  for(i in seq(1, 4)){
    S[i] = SATE(Y = Y[which(X == i)], A = A[which(X == i)]) 
  }
  result = which.max(S) # the larger the S, the larger the treatment effect
  return(result)
}

# Load data
# Sharing of the data is upon request because it contains sensitive information that could potentially compromise the privacy of individuals involved. 
# baselineData = read.csv("XXX")
# outcomeData = read.csv("XXX")
baselineData$outcome = outcomeData$EVENT_PRIMARY

# Create a hypothetical two-stage enrichment trial based on the SPRINT trial
set.seed(31)
N1 = 2000 
stage1.index = sample(dim(baselineData)[1], N1)
N2 = 200
stage2.index = sample(setdiff(which(baselineData$AGE >= 80), stage1.index), N2)
age.breaks = c(0, 59, 69, 79, 100)
realData = list()
realData$data.stage.1 = list(X = as.numeric(cut(baselineData$AGE, breaks = age.breaks))[stage1.index], 
                           Y = baselineData$outcome[stage1.index], 
                           A = baselineData$INTENSIVE[stage1.index])
realData$data.stage.2 =  list(X = as.numeric(cut(baselineData$AGE, breaks = age.breaks))[stage2.index], 
                              Y = baselineData$outcome[stage2.index], 
                              A = baselineData$INTENSIVE[stage2.index])
realData$inclusion = inclusion.function(Y = realData$data.stage.1$Y,
                                        A = realData$data.stage.1$A,
                                        X = realData$data.stage.1$X)
realData$hypothesis = realData$inclusion
realData$test.statistic = SATE(Y = c(realData$data.stage.1$Y[which(realData$data.stage.1$X %in% realData$hypothesis)], realData$data.stage.2$Y), 
                               A = c(realData$data.stage.1$A[which(realData$data.stage.1$X %in% realData$hypothesis)], realData$data.stage.2$A))
