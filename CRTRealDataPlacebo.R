# CRT on a hypothetical two-stage enrichment trial based on the SPRINT trial. 
# Placebo analysis

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

inclusion.set = function(inclusion){return(inclusion)}

# Load data
# Sharing of the data is upon request because it contains sensitive information that could potentially compromise the privacy of individuals involved. 
# baselineData = read.csv("XXX")
# outcomeData = read.csv("XXX")
baselineData$outcome = outcomeData$EVENT_PRIMARY

# Create a hypothetical two-stage enrichment trial based on the SPRINT trial
set.seed(318)
N1 = 2000
N2 = 200
age.breaks = c(0, 59, 69, 79, 100)
baselineData$AGE = as.numeric(cut(baselineData$AGE, breaks = age.breaks))
baselineData$AGE[baselineData$AGE == 4] = 1
baselineData$AGE[baselineData$AGE == 3] = 2
baselineData = baselineData[baselineData$INTENSIVE == 0, ]
tau = 0
m = 1000

M = 1000; count = 1
result.realData = data.frame("CRT" = numeric(M), "RT.second.stage" = numeric(M), "RT" = numeric(M), "lm" = numeric(M), "lm.second.stage" = numeric(M), "glm" = numeric(M), "glm.second.stage" = numeric(M), "inclusion" = numeric(M))
while(count <= M){
  baselineData$INTENSIVE = rbinom(length(baselineData$INTENSIVE), 1, 0.5)
  stage1.index = sample(dim(baselineData)[1], N1)
  realData = list()
  realData$data.stage.1 = list(X = baselineData$AGE[stage1.index], 
                               Y = baselineData$outcome[stage1.index], 
                               A = baselineData$INTENSIVE[stage1.index])
  realData$inclusion = inclusion.function(Y = realData$data.stage.1$Y,
                                          A = realData$data.stage.1$A,
                                          X = realData$data.stage.1$X)
  stage2.index = sample(setdiff(which(baselineData$AGE == realData$inclusion), stage1.index), N2)
  realData$data.stage.2 =  list(X = baselineData$AGE[stage2.index], 
                                Y = baselineData$outcome[stage2.index], 
                                A = baselineData$INTENSIVE[stage2.index])
  realData$hypothesis = realData$inclusion
  result.realData$inclusion[count] = realData$inclusion
  
  realData$test.statistic = SATE(Y = c(realData$data.stage.1$Y[which(realData$data.stage.1$X %in% realData$hypothesis)], realData$data.stage.2$Y), 
                                 A = c(realData$data.stage.1$A[which(realData$data.stage.1$X %in% realData$hypothesis)], realData$data.stage.2$A)) # relative ratio
  
  # CRT
  result.realData$CRT[count] = inference.enrichment.analysis(result = realData, 
                                                             m.permutation = m,
                                                             tau = tau,
                                                             method = "CRT",
                                                             sampler = "rejection.sampler",
                                                             count.max = 100 * m)
  result.realData$RT.second.stage[count] = inference.enrichment.analysis(result = realData, 
                                                                         m.permutation = m,
                                                                         tau = tau,
                                                                         method = "RT.second.stage")
  result.realData$RT[count] = inference.enrichment.analysis(result = realData, 
                                                            m.permutation = m,
                                                            tau = tau,
                                                            method = "RT")
  
  # model-based
  data.combined.ageGroup = data.frame(Y = c(realData$data.stage.1$Y[realData$data.stage.1$X == realData$inclusion], realData$data.stage.2$Y[realData$data.stage.2$X == realData$inclusion]),
                                      A = c(realData$data.stage.1$A[realData$data.stage.1$X == realData$inclusion], realData$data.stage.2$A[realData$data.stage.2$X == realData$inclusion]),
                                      stage = c(rep(1, sum(realData$data.stage.1$X == realData$inclusion)), rep(2, sum(realData$data.stage.2$X == realData$inclusion))))
  
  # lm
  model = lm(Y ~ A, data.combined.ageGroup)
  result.realData$lm[count] = summary(model)[[4]][2,4] 
  model.second.stage = lm(Y ~ A, data.combined.ageGroup[data.combined.ageGroup$stage == 2,])
  result.realData$lm.second.stage[count] = summary(model.second.stage)[[4]][2, 4]
  # glm
  model = glm(Y ~ A, data.combined.ageGroup, family = "binomial")
  result.realData$glm[count] = summary(model)[[12]][2, 4]
  model.second.stage = glm(Y ~ A, data.combined.ageGroup[data.combined.ageGroup$stage == 2,], family = "binomial")
  result.realData$glm.second.stage[count] = summary(model.second.stage)[[12]][2, 4]
  
  if(i %% 10 == 0){print(count)}
  count = count + 1
}

