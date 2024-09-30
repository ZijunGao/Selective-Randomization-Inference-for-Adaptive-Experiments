# CRT on a hypothetical two-stage enrichment trial based on the SPRINT trial. 

# helper function
source("helper.R")
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
# realData = readRDS("XXX") 

# CRT, RT.second.stage
tau = 0
m = 1000
result.realData = list()

set.seed(10)
(result.realData$CRT = inference.enrichment.analysis(result = realData, 
                                               m.permutation = m,
                                               tau = tau,
                                               method = "CRT",
                                               sampler = "rejection.sampler",
                                               count.max = 100 * m))
(result.realData$RT.second.stage = inference.enrichment.analysis(result = realData, 
                                                 m.permutation = m,
                                                 tau = tau,
                                                 method = "RT.second.stage"))
(result.realData$RT = inference.enrichment.analysis(result = realData, 
                                                   m.permutation = m,
                                                   tau = tau,
                                                   method = "RT"))


# Results
print(result.realData)

# Model-based p-value + Bonferroni correction
# Data
data.combined.ageGroup = data.frame(Y = c(realData$data.stage.1$Y[realData$data.stage.1$X == 4], realData$data.stage.2$Y[realData$data.stage.2$X == 4]),
                                    A = c(realData$data.stage.1$A[realData$data.stage.1$X == 4], realData$data.stage.2$A[realData$data.stage.2$X == 4]),
                                    stage = c(rep(1, sum(realData$data.stage.1$X == 4)), rep(2, sum(realData$data.stage.2$X == 4))))

# lm
model = lm(Y ~ A, data.combined.ageGroup)
summary(model)
summary(model)[[4]][2,4] # Bonferroni's correction; a finer partition into 10 groups, * 10
model.second.stage = lm(Y ~ A, data.combined.ageGroup[data.combined.ageGroup$stage == 2,])
summary(model.second.stage)[[4]][2, 4]
# glm
model = glm(Y ~ A, data.combined.ageGroup, family = "binomial")
summary(model)[[12]][2, 4] # Bonferroni's correction; a finer partition into 4 groups, * 4
model.second.stage = glm(Y ~ A, data.combined.ageGroup[data.combined.ageGroup$stage == 2,], family = "binomial")
summary(model.second.stage)[[12]][2, 4]
