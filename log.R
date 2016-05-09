####################################
### Load Libraries and Functions ###
####################################
## R version: 3.2.3
source("./src.R") # source functions 


#############################
### Initialize Parameters ###
#############################
## grid resolution in arbitrary spatial units (in this case km)
res <- 10

## matrix initialization
boosts <- matrix(1,100,100)
seeds <- matrix(NA,100,100)

## adding seeds
seeds[50,25] <- 0
seeds[50,75] <- 0

## base speeds
speeds <- c(1,3)


####################################
### Run Fast Marching Algorithms ###
####################################
## simple case, with no boost factor
cdist1 <- MFM(seeds, speeds, res, boosts)
plot(cdist1)

## case with barrier to movement in-between two sources
boosts[25:75,45:55] <- 0
cdist2 <- MFM(seeds, speeds, res, boosts)
plot(cdist2)

## case with randomly located sources and corresponding speeds
boosts <- matrix(1,100,100)
n.seeds <- 5
seeds <- matrix(NA,100,100)
seeds[matrix(c(ceiling(runif(n.seeds)*100),ceiling(runif(n.seeds)*100)),ncol=2)] <- 0
speeds <- ceiling(runif(n.seeds)*5)
cdist3 <- MFM(seeds, speeds, res, boosts)
plot(cdist3)
