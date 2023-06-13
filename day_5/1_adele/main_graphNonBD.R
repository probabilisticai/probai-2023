rm(list=ls())

###################################
# Installing required  R packages #
###################################

# install.packages("BiocManager")
# BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
#
# install.packages(c("igraph", "pcalg", "dagitty"), dependencies=TRUE)
# install.packages("causaleffect", dependencies=TRUE)

# install.packages("devtools", dependencies=TRUE)
# devtools::install_github("adele/PAGId", dependencies=TRUE)

##############################################

library(pcalg)
library(dagitty)
library(MASS)
library(matrixcalc)
library(causaleffect)
library(igraph)
library(PAGId)

source("helperFunctions.R")

#########################################
# Constructing the true Causal Diagram  #
#########################################

# Getting the true causal diagram
trueDAG.out <- getDAG2bR42()
trueDAG.dag <- trueDAG.out$adag # a dagitty object

# Note: this plots the DAG including the U's variables
plot(trueDAG.dag)

# Getting the ADMG corresponding to the true DAG
trueDAG.amat <- trueDAG.out$amat # the adjacency matrix following the pcalg notation
trueDAG.gNEL <- as(t(trueDAG.out$amat), "graphNEL")
trueADMG.ig <- igraph_from_graphNel(trueDAG.gNEL, latents(trueDAG.dag))

# Note: this plots the causal diagram with bidirected edges indicating unmeasured confounders
plot(trueADMG.ig)

# Getting the minimal conditional independence relationships over the observed variables,
# implied by the true ADMG -- obtained by minimal d-separators
trueImpliedCI <- dagitty::impliedConditionalIndependencies(trueDAG.dag, type = "missing.edge")
trueImpliedCI


################################################################
# Checking if P(y1,y2,y3,y4,y5|do(x1,x2)) is identifiable      #
# from the ADMG/MAG using the generalized adjustment criterion #
################################################################

x <- 1:2
y <- 3:7

mag <- dagitty::toMAG(trueDAG.dag)
amat.mag <- dagitty2amat(mag, type="mag")

adj <- adjustment(amat = amat.mag,
                  amat.type = "mag", x = x, y = y,
                  set.type = "all")
adj # There is no valid adjustment set for X={X1,X2} and Y={Y1,Y2,Y3,Y4,Y5}


########################################################
# Identifying causal effect P(y |do(x)) from the ADMG  #
########################################################

x <- c("x1", "x2")
y <- c("y1", "y2", "y3", "y4", "y5")
z <- NULL


exprDAG <- causaleffect::causal.effect(y=y, x=x, z=z,
                                       G = trueADMG.ig,
                                       expr = TRUE,
                                       simp = TRUE)
print(paste0("ID from DAG: ", exprDAG))

# to see all steps applied by the ID algorithm,
# use expr = FALSE and steps = TRUE.
retDAG <- causaleffect::causal.effect(y=y, x=x, z=z,
                                      G = trueADMG.ig,
                                      expr = FALSE,
                                      simp = TRUE,
                                      steps = TRUE)
print(retDAG)



#####################################
# Simulating data from the true DAG #
#####################################

N = 10000 # try different values of N
#seed = ceiling(runif(1, 0, 10000))
seed = 3632
set.seed(seed)

# The following function generates a random linear Gaussian SCM
# compatible with the true DAG and then samples from it:
dat <- generateDatasetFromDAG(adag = trueDAG.dag, N=N)
head(dat)


#########################################################
# Testing some expected dependencies and independencies #
#########################################################

# For Gaussian data, partial correlation tests (Fisher's Z test)
# can be used to assess conditional independencies
# For other data types, appropriate conditional independence tests should be used.
indepTest <- pcalg::gaussCItest
alpha <- 0.05
suffStat <- list(C = cor(dat), n = N)

X1ind <- 1
X2ind <- 2
Y1ind <- 3
Y2ind <- 4
Y3ind <- 5
Y4ind <- 6
Y5ind <- 7

# some dependencies:
p_x1y4.y3 <- indepTest(X1ind, Y4ind, Y3ind, suffStat)
p_x1y5.y1 <- indepTest(X1ind, Y5ind, Y1ind, suffStat)
p_x2y3.y2 <- indepTest(X2ind, Y3ind, Y2ind, suffStat)
p_y3y1.y5 <- indepTest(Y3ind, Y1ind, Y5ind, suffStat)
p_y3y5 <- indepTest(Y3ind, Y5ind, c(), suffStat)

ret1 <- all(c(p_x1y4.y3, p_x1y5.y1, p_x2y3.y2, p_y3y1.y5, p_y3y5) < alpha)
ret1

# minimal indepedencies
p_x1y2.x2 <- indepTest(X1ind, Y2ind, X2ind, suffStat)
p_x1y4 <- indepTest(X1ind, Y4ind, c(), suffStat)
p_x1y5 <- indepTest(X1ind, Y5ind, c(), suffStat)
p_x2y3 <- indepTest(X2ind, Y3ind, c(), suffStat)
p_x2y4 <- indepTest(X2ind, Y4ind, c(), suffStat)
p_x2y5 <- indepTest(X2ind, Y5ind, c(), suffStat)

p_y1y2.x2 <- indepTest(Y1ind, Y2ind, X2ind, suffStat)
p_y1y3.x1y4 <- indepTest(Y1ind, Y3ind, c(X1ind, Y4ind), suffStat)
p_y1y3.x1y5 <- indepTest(Y1ind, Y3ind, c(X1ind, Y5ind), suffStat)
p_y1y4.y5 <- indepTest(Y1ind, Y4ind, Y5ind, suffStat)

p_y2y4 <- indepTest(Y2ind, Y4ind, c(), suffStat)
p_y2y5 <- indepTest(Y2ind, Y5ind, c(), suffStat)
p_y3y5.y4 <- indepTest(Y3ind, Y5ind, Y4ind, suffStat)

ret2 <- all(c(p_x1y2.x2, p_x1y4, p_x1y5, p_x2y3, p_x2y4, p_x2y5, p_y1y2.x2,
      p_y1y3.x1y4, p_y1y3.x1y5, p_y1y4.y5, p_y2y4, p_y2y5, p_y3y5.y4) >= alpha)
ret2

###############################
# Estimating a PAG using FCI  #
###############################

# Get the PAG using the true ADMG as an oracle of cond. indep. relationships
truePAG <- getTruePAG(trueDAG.dag, verbose = TRUE)
plot(truePAG)


# For Gaussian data, partial correlation tests (Fisher's Z test)
# can be used to assess conditional independencies
# For other data types, appropriate conditional independence tests should be used.
indepTest <- pcalg::gaussCItest
alpha <- 0.05
suffStat <- list(C = cor(dat), n = N)
estPAG <- pcalg::fci(suffStat,
                     indepTest = indepTest,
                     labels= colnames(dat), alpha = alpha,
                     verbose = TRUE)
plot(estPAG)

# True and estimated PAGs are identical
ret3 <- all(truePAG@amat - estPAG@amat == 0)
ret3

# Getting the minimal conditional independence relationships
# implied by the estimated PAG
# This can be obtained by computing the minimal d-separators in
# any MAG that is member of the PAG
amat.mag <- pcalg::pag2magAM(estPAG@amat, 1)
plotAG(amat.mag)
amag <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")
estImpliedCI <- dagitty::impliedConditionalIndependencies(amag)
estImpliedCI


###########################################################
# Checking if P(y1,y2,y3,y4,y5|do(x1,x2)) is identifiable #
# from the PAG using generalized adjustment criterion.    #
###########################################################

#P_{x1,x2}(y1,y2,y3,y4,y5)
Xinds <- c(1,2)
Yinds <- 3:7

x <- Xinds
y <- Yinds

# adjacency matrix of the estimated PAG following the pcalg notation
estPAG.amat <- estPAG@amat

adj <- adjustment(amat = estPAG.amat, amat.type = "pag", x = x, y = y, set.type = "all")
adj # There is no valid adjustment set for X={A, F} and Y={Y}


# Testing if the particular set {} is valid for adjustment:
z <- c()
gac_z <- gac(estPAG.amat, x, y, z, type = "pag")
gac_z$gac # the set is not valid for adjustment

# the third condition of the gac is violated, i.e.,
# not all proper definite status non-causal paths from x to y are blocked by z
gac_z$res



########################################################
# Identifying P(y1,y2,y3,y4,y5|do(x1,x2)) from the PAG #
# using the CIDP algorithm                             #
########################################################

x <- c("x1", "x2")
y <- c("y1", "y2", "y3", "y4", "y5")
z <- NULL


# adjacency matrix of the estimated PAG following the pcalg notation
estPAG.amat <- estPAG@amat

retPAG <- CIDP(estPAG.amat, x, y, z)
retPAG$id
print(paste0("ID from PAG: ", retPAG$Qexpr[[retPAG$query]]))

# This shows the steps taken by the CIDP algorithm
# by substitution and simplification, we will get the same
# identification formula from the ADMG
retPAG$Qexpr


