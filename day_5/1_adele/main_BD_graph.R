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

#############################################
# Constructing the Structural Causal Model  #
# and corresponding Causal Diagram          #
#############################################

# Getting the true SCM - follow the BD model
# A -> B -> C;
# B <- D <-> C
trueSCM <- getBDSEM()

# Getting the true causal diagram - follow the BD model
trueDAG.out <- getBDGraph()
trueDAG.dag <- trueDAG.out$adag # a dagitty object

# Note: this plots the DAG including the U's variables
plot(trueDAG.dag)

# Getting the ADMG corresponding to the true DAG
trueDAG.amat <- trueDAG.out$amat # the adjacency matrix following the pcalg notation
trueDAG.gNEL <- as(t(trueDAG.out$amat), "graphNEL")
trueADMG.ig <- igraph_from_graphNel(trueDAG.gNEL, latents(trueDAG.dag))

# Note: this plots the causal diagram with bidirected edges indicating unmeasured confounders
plot(trueADMG.ig)


###################################################
# Checking the conditional independence relations #
# implied by the true Causal Diagram              #
###################################################

# Getting the minimal conditional independence relationships over the observed variables,
# implied by the true ADMG -- obtained by minimal d-separators
trueImpliedCI <- dagitty::impliedConditionalIndependencies(trueDAG.dag, type = "missing.edge")
print(trueImpliedCI)


##########################################
# Checking if P(c|do(b)) is identifiable #
# from the ADMG/MAG using the            #
# generalized adjustment criterion       #
##########################################

Aind <- 1
Bind <- 2
Cind <- 3
Dind <- 4

x <- Bind
y <- Cind

mag <- dagitty::toMAG(trueDAG.dag)
amat.mag <- dagitty2amat(mag, type="mag")
plotAG(amat.mag)

adj <- adjustment(amat = amat.mag,
                  amat.type = "mag", x = x, y = y,
                  set.type = "all")
print(adj) # The valid adjustment sets for X={B} and Y={C} are {D} and {A,D}


#######################################################
# Identifying causal effect P(c|do(b)) from the ADMG  #
#######################################################

y = "C"
x = "B"
z = c()
exprDAG <- causaleffect::causal.effect(y=y, x=x, z=z,
                                       G = trueADMG.ig,
                                       expr = TRUE,
                                       simp = TRUE)
print(paste0("ID from DAG: ", exprDAG))

# Note: As expected, P(c|do(b)) is identifiable
# through adjustment.

# to see all steps applied by the ID algorithm,
# use expr = FALSE and steps = TRUE.
retDAG <- causaleffect::causal.effect(y=y, x=x, z=z,
                                      G = trueADMG.ig,
                                      expr = FALSE,
                                      simp = TRUE,
                                      steps = TRUE)
#print(retDAG)



##################################################
# Simulating data from the true SCM and true DAG #
##################################################

N = 10000 # try different values of N
#seed = ceiling(runif(1, 0, 10000))
seed = 831
set.seed(seed)
dat <- generateDatasetFromSEM(trueSCM$beta, trueSCM$lat, N)
head(dat)
dat <- dat[,c("A", "B", "C", "D")]

# Note: we could also simulate data from the true DAG
# The following function generates a random linear Gaussian SCM
# compatible with the true DAG and then samples from it:
dat2 <- generateDatasetFromDAG(adag = trueDAG.dag, N=N)
head(dat2)


#########################################################
# Testing some expected dependencies and independencies #
#########################################################

# For Gaussian data, partial correlation tests (Fisher's Z test)
# can be used to assess conditional independencies
# For other data types, appropriate conditional independence tests should be used.
indepTest <- pcalg::gaussCItest
alpha <- 0.05
suffStat <- list(C = cor(dat), n = N)

Aind <- 1
Bind <- 2
Cind <- 3
Dind <- 4

# dependencies:
p_AC.B <- indepTest(Aind, Cind, Bind, suffStat)
p_AC <- indepTest(Aind, Cind, c(), suffStat)
p_AD.B <- indepTest(Aind, Dind, Bind, suffStat)

print(all(c(p_AC, p_AC.B, p_AD.B) < alpha))

# indepedencies
p_AC.BD <- indepTest(Aind, Cind, c(Bind, Dind), suffStat)
p_AD <- indepTest(Aind, Dind, c(), suffStat)

print(all(c(p_AD, p_AC.BD) >= alpha))


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
print(all(truePAG@amat - estPAG@amat == 0))

###################################################
# Checking the conditional independence relations #
# implied by the estimated PAG                    #
###################################################

# Getting the minimal conditional independence relationships
# implied by the estimated PAG
# This can be obtained by computing the minimal d-separators in
# any MAG that is member of the PAG
amat.mag <- pcalg::pag2magAM(estPAG@amat, 1)
plotAG(amat.mag)
amag <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")

estImpliedCI <- dagitty::impliedConditionalIndependencies(amag)
print(estImpliedCI)


######################################################
# Identifying causal effect P(c|do(b)) from the PAG  #
# using generalized adjustment criterion.            #
######################################################

Aind <- 1
Bind <- 2
Cind <- 3
Dind <- 4

y = Cind
x = Bind

# adjacency matrix of the estimated PAG following the pcalg notation
estPAG.amat <- estPAG@amat

adj <- adjustment(amat = estPAG.amat, amat.type = "pag", x = x, y = y, set.type = "all")
print(adj) # The valid adjustment sets are  {D} and {A, D}



z <- c() # empty set
gac_empty <- gac(estPAG.amat, x, y, z, type = "pag")
print(gac_empty$gac) # the empty set is not valid for adjustment

# the third condition of the gac is violated, i.e.,
# not all proper definite status non-causal paths from x to y are blocked by z
print(gac_empty$res)

z <- adj[[1]] # set {D}
gac_d <- gac(estPAG.amat, x, y, z, type = "pag")
print(gac_d$gac) # the set {D} is valid for adjustment

# all conditions of the gac are satisfied
print(gac_d$res)



######################################################
# Identifying causal effect P(c|do(b)) from the PAG  #
# using the CIDP algorithm                           #
######################################################

y = "C"
x = "B"
z = c()

# adjacency matrix of the estimated PAG following the pcalg notation
estPAG.amat <- estPAG@amat

retPAG <- CIDP(estPAG.amat, x, y, z)
print(retPAG$id)
print(paste0("ID from PAG: ", retPAG$Qexpr[[retPAG$query]]))

# This shows the steps taken by the CIDP algorithm
# by substitution and simplication, we will get the same adjustment formula
print(retPAG$Qexpr)


