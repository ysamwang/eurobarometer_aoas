### Example code to run analysis of eurobarometer data ###

#### Setup #####
seed <- 1
# parseCommandArgs takes in a seed and overwrites the currently stored value
parseCommandArgs()
set.seed(seed)

############
library(mixedMem)
library(batch)
library(gtools)

# Number of sub-groups
K <- 5
# dirichlet parameter for initializing theta
thet <- .6

# Total number of individuals
Total <- 11872

## Number of choices to look at for chi-squared metric
mult <- 4
### Number of parametric bootstraps to run for chi-squared metric
sim.size <- 200
## Number of choices to include in model (7 is the max across all questions)
limit <- 7



## Read in observations ##
data.drug.fix <- readRDS("~/mixedMem/euro/drugFix.rds")
data.alc.fix <- readRDS("~/mixedMem/euro/alcFix.rds")
data.aid.fix <- readRDS("~/mixedMem/euro/aidFix.rds")
N <- readRDS("~/mixedMem/euro/nijr.rds")

# Number of questions
J <- 3
# 1 repeated measure per question
Rj <- rep(1,J)
# Number of choices for each question
Vj <- c(7,10,5)

# Initalization 
alpha <- rep(.1, K)
theta <- array(0, dim = c(J, K, max(Vj)))

# All questions are rank data
dist <- rep("rank", J)

# Number of observed
Nijr <- array(0, dim = c(Total, J, max(Rj)))
Nijr[c(1:Total),c(1:J),1] <- N[1:11872, c(1,3,4), 1]
Nijr <- ifelse(Nijr  < limit, Nijr, limit)
obs <- array(-1, dim = c(Total,J, max(Rj), max(Nijr)))


# Read in observations
j <- 1
r <- 1
count <- 1
for (i in 1:11872)
{
  for(v in 1:(dim(data.drug.fix)[2]))
  {
    if(!is.na(data.drug.fix[i,v]) & (data.drug.fix[i,v] <= Nijr[count,j,r]))
    {
      obs[count,j,r,data.drug.fix[i,v]] <- (v-1) 
    }
  }
  count <- count + 1
}



j <- 2
r <- 1
count <- 1
for (i in 1:11872)
{
  for(v in 1:(dim(data.alc.fix)[2]))
  {
    if(!is.na(data.alc.fix[i,v]) & (data.alc.fix[i,v] <= Nijr[count,j,r]))
    {
      obs[count,j,r,data.alc.fix[i,v]] <- (v-1) 
    }
  }
  count <- count + 1
}

j <- 3
r <- 1
count <- 1
for(i in 1:11872)
{
  for(v in 1:(dim(data.aid.fix)[2]))
  {
    if(!is.na(data.aid.fix[i,v]) & (data.aid.fix[i,v] <= Nijr[count,j,r]))
    {
      obs[count,j,r,data.aid.fix[i,v]] <- (v-1) 
    }
  }
  count <- count + 1
}

### Get ground truth counts for comparison later 
cutoff = 5
top.drug.fix.truth <- sort(table(apply(obs[, 1, 1,c(1:mult)], MAR = 1, FUN = paste, collapse = "")), decreasing = T)
top.alc.fix.truth <- sort(table(apply(obs[, 2, 1,c(1:mult)], MAR = 1, FUN = paste, collapse = "")), decreasing = T)
top.aid.fix.truth <- sort(table(apply(obs[, 3, 1,c(1:mult)], MAR = 1, FUN = paste, collapse = "")), decreasing = T)

top.drug.fix.truth <- top.drug.fix.truth[which(top.drug.fix.truth > cutoff)]
top.alc.fix.truth <- top.alc.fix.truth[which(top.alc.fix.truth > cutoff)]
top.aid.fix.truth <- top.aid.fix.truth[which(top.aid.fix.truth > cutoff)]

### Initialize Theta ### 
r = .01
for(j in 1:J) {
  theta[j,,c(1:Vj[j])] <- gtools::rdirichlet(K, rep(thet, Vj[j]))
  
  # uniform group
  theta[j,K-1,c(1:Vj[j])] <- rep(1, Vj[j])/Vj[j]
  
  # group corresponding to the presented ordering
  a <- (1-r)/(1-r^Vj[j])
  theta[j,K,c(1:Vj[j])] <- a*r^c(0:(Vj[j]-1))
}


### Create the mixedMemModel object ###
mm <- mixedMemModel(Total = Total,
                    J = J,
                    Rj = Rj, Nijr = Nijr,
                    K = K, Vj = Vj,
                    alpha = alpha, theta = theta,
                    dist = dist, obs = obs)

### Fit the model ###
st = proc.time()
out <- mmVarFit(mm, elboTol = 1e-5,
                bMax = 3, bNaught = 1000, maxThetaIter = 2000, holdConst = c(K-2,K-1),
                maxTotalIter = 1000)
end = proc.time()

elbo <- computeELBO(out)


### Evaluate chi-squared statistic between observed rankings and parametric bootstrap rankings###

exp.drug.fix <- rep(0, length(top.drug.fix.truth))
exp.alc.fix <- rep(0, length(top.alc.fix.truth))
exp.aid.fix <- rep(0, length(top.aid.fix.truth))

countOccurence <- function(x, y) {
  sum(x == y)
}

for(i in 1:sim.size) {
  obs1 <- mixedMem::rmixedMem(Total = Total, J =  J,Rj =  Rj,Nijr =  Nijr,
                              K = K,Vj = Vj, dist =  dist, theta =  out$theta,
                              alpha = out$alpha, lambda = out$phi/rowSums(out$phi))$obs
  top.drug.fix <- apply(obs1[, 1, 1,c(1:mult)], MAR = 1, FUN = paste, collapse = "")
  top.alc.fix <- apply(obs1[, 2, 1,c(1:mult)], MAR = 1, FUN = paste, collapse = "")
  top.aid.fix <- apply(obs1[, 3, 1,c(1:mult)], MAR = 1, FUN = paste, collapse = "")
  
  exp.drug.fix <- exp.drug.fix + apply(matrix(names(top.drug.fix.truth), ncol = 1), MAR = 1, FUN = countOccurence, top.drug.fix)
  exp.alc.fix <- exp.alc.fix + apply(matrix(names(top.alc.fix.truth), ncol = 1), MAR = 1, FUN = countOccurence, top.alc.fix)
  exp.aid.fix <- exp.aid.fix + apply(matrix(names(top.aid.fix.truth), ncol = 1), MAR = 1, FUN = countOccurence, top.aid.fix)
}


exp.drug.fix <- exp.drug.fix / sim.size
exp.alc.fix <- exp.alc.fix / sim.size
exp.aid.fix <- exp.aid.fix / sim.size


chi1 <- sum((top.drug.fix.truth - exp.drug.fix)^2 / exp.drug.fix)
chi2 <- sum((top.alc.fix.truth - exp.alc.fix)^2 / exp.alc.fix)
chi3 <- sum((top.aid.fix.truth - exp.aid.fix)^2 / exp.aid.fix)


### Write result ###
result <- data.frame(seed = seed, elbo = elbo,
                     alpha = max(out$alpha), K = K, 
                     chi1 = chi1, chi2 = chi2, chi3 = chi3, chiSum = chi1 + chi2 + chi3,
                     t.init = thet,
                     time = (end-st)[1])

result_name = paste("~/mixedMem/euro/results/res",seed, ".csv", sep = "")
write.csv(result, result_name, row.names = F)