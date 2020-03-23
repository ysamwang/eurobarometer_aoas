args <- commandArgs(TRUE)

eval(parse(text=args[[1]]))

library(mixedMem)
library(gtools)

  set.seed(seed)
  opt.model <- readRDS("~/mixedMem/euro/optFixed2_q3HOfull.rds")
  
  
  Total <- 11872
  sub.sample <- sample.int(Total, size = Total, replace = T)
  J <- 3
  Vj <- c(7, 10, 5)
  K <- 5
  alpha <- opt.model$alpha
  theta <- opt.model$theta
  dist <- rep("rank", J)
  Nijr<- opt.model$Nijr[sub.sample, , , drop = FALSE]
  obs <- opt.model$obs[sub.sample, , , , drop = FALSE]
  phi <- opt.model$phi[sub.sample, , drop = FALSE]
  
  
  mm <- mixedMemModelVarInf(Total = Total, J = J, Nijr = Nijr, 
                      K = K, Vj = Vj, phi = phi, 
                      alpha = alpha, theta = theta,
                      dist = dist, obs = obs)
  st = proc.time()
  out <- mmVarInfFit(mm, elboTol = 1e-5,
                  bMax = 2, bNaught = 1000, bMult = 5000,  maxThetaIter = 2000, holdConst = c(K-2,K-1),
                  maxTotalIter = 1000)
  end = proc.time()
  
  elbo <- computeELBO(out)
  
  
  res <- matrix(c(seed, elbo, (end - st)[3], out$alpha, out$theta), nrow = 1)
  
  result_name = paste("~/mixedMem/euro/bootstrap/res_old_init/bs",seed, ".csv", sep = "")
  write.csv(res, result_name, row.names = F)

