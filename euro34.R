# This file cleans the raw data from the eurobarometer website into the data which is ultimately used
# in the paper.


library(foreign)
library(mixedMem)
setwd("C:/Users/Sam/Dropbox/1 RA - MMS/derivations/euroBarometer")
dat <- read.dta("ZA1961_v1-0-1.dta")


# Take out people who didn't rank anything
drug.fix.comp <- (!is.na(dat$v74) | !is.na(dat$v76) | !is.na(dat$v78) | !is.na(dat$v80)|
               !is.na(dat$v82) |!is.na(dat$v84)|!is.na(dat$v86)) 
sum(drug.fix.comp)
alc.consq.comp <- (!is.na(dat$v139) | !is.na(dat$v140) | !is.na(dat$v141)| !is.na(dat$v142)|
                !is.na(dat$v143) | !is.na(dat$v144)| !is.na(dat$v145))
sum(alc.consq.comp)
alc.fix.comp <- (!is.na(dat$v146)| !is.na(dat$v148) | !is.na(dat$v150)| !is.na(dat$v152) |
              !is.na(dat$v154) | !is.na(dat$v156) | !is.na(dat$v158)| !is.na(dat$v160)|
            !is.na(dat$v162) | !is.na(dat$v164))
sum(alc.fix.comp)
aid.fix.comp <- (!is.na(dat$v208)|!is.na(dat$v210)|!is.na(dat$v212)|!is.na(dat$v214)|
              !is.na(dat$v216))
sum(aid.fix.comp)

indices <- drug.fix.comp & alc.consq.comp & alc.fix.comp & aid.fix.comp
sum(indices)


# Take out people who repeated rankings
check.repeat <- function( x )
{
  x <- na.omit(x)
  return(all(sort(as.numeric(substr(x,1,1)))==c(1:length(x))))
}

data <- data.frame(dat[indices,])
data.drug.fix <- data[,76+c(0:6)*2]
data.alc.consq <- data[,c(141:147)]
data.alc.fix <- data[,148+c(0:9)*2]
data.aid.fix <- data[,210+ c(0:4)*2]
rep1 <- apply(data.drug.fix, MAR = 1, FUN = check.repeat)
rep2 <- apply(data.alc.consq, MAR = 1, FUN = check.repeat)
rep3 <- apply(data.alc.fix, MAR = 1, FUN = check.repeat)
rep4 <- apply(data.aid.fix, MAR = 1, FUN = check.repeat)
repeats <- rep1 & rep2 & rep3 & rep4
data <- data[repeats,]
data.drug.fix <- data.drug.fix[repeats,]
data.alc.consq <- data.alc.consq[repeats,]
data.alc.fix <- data.alc.fix[repeats,]
data.aid.fix <- data.aid.fix[repeats,]

make.num <- function(x)
{
  x <- as.numeric(x)
  if(min(x, na.rm = T)>1){
    x <- x-1
  }
  return(x)
}
data.drug.fix1 <- matrix(0, nrow =  dim(data.drug.fix)[1], ncol = dim(data.drug.fix)[2])
data.alc.consq1 <- matrix(0, nrow =  dim(data.drug.fix)[1], ncol = dim(data.alc.consq)[2])
data.alc.fix1 <- matrix(0, nrow =  dim(data.drug.fix)[1], ncol = dim(data.alc.fix)[2])
data.aid.fix1 <- matrix(0, nrow =  dim(data.drug.fix)[1], ncol = dim(data.aid.fix)[2])


for(i in 1:dim(data.drug.fix)[1])
{
  data.drug.fix1[i,] <- make.num(data.drug.fix[i,])
  data.alc.consq1[i,] <- make.num(data.alc.consq[i,])
  data.alc.fix1[i,] <- make.num(data.alc.fix[i,])
  data.aid.fix1[i,] <- make.num(data.aid.fix[i,])
}


data.combined <-data.frame(data.drug.fix1, data.alc.consq1, data.alc.fix1, data.aid.fix1)
names(data.combined) <- c(paste("drugFix",c(1:dim(data.drug.fix)[2]), sep = ""),
                          paste("alcConsq",c(1:dim(data.alc.consq)[2]), sep = ""),
                          paste("alcFix",c(1:dim(data.alc.fix)[2]), sep = ""),
                          paste("aidFix",c(1:dim(data.aid.fix)[2]), sep = ""))

write.table(data.combined,"cleanEuro34.txt")

