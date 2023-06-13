library(pcalg)
library(dagitty)
library(MASS)
library(matrixcalc)
library(causaleffect)
library(igraph)
library(PAGId)

####################
# Helper Functions #
####################

# A <- f(Ua) = Ua
# D <- f(Ucd, Ud) = beta_D.Ucd * Ucd + Ud
# B <- f(A, D, Ub) = beta_B.A * A - beta_B.D * D + Ub
# C <- f(Ucd, B, Uc) = beta_C.B * B - beta_C.Ucd * Ucd + Uc
# betaList is a list with entries for:
# beta_D.Ucd, beta_B.A, beta_B.D, beta_C.Ucd, and beta_C.B
# if betaList is NULL, then coefficients are randomly selected
getBDSEM <- function(betaList=NULL) {
  allvars <- c("A", "B", "C", "D", "Ucd")
  p <- length(allvars)
  beta <- matrix(0, p, p)
  colnames(beta) <- rownames(beta) <- allvars
  topolOrd <- c( "Ucd", "A", "D", "B", "C")

  if (is.null(betaList)) {
    betaList <- as.list(sample(c(-1, 1), p, replace = T) * runif(p, 0.3, 0.7))
    names(betaList) <- c("beta_D.Ucd", "beta_B.A", "beta_B.D", "beta_C.Ucd", "beta_C.B")
  }

  beta["Ucd","D"] <- betaList$beta_D.Ucd; # coeff of Ucd in the function for d
  beta["A", "B"] <- betaList$beta_B.A; # coeff of a in the function for b
  beta["D", "B"] <- betaList$beta_B.D; # coeff of d in the function for b
  beta["Ucd", "C"] <- betaList$beta_C.Ucd; # coeff of Ucd in the function for c
  beta["B", "C"] <- betaList$beta_C.B; # coeff of b in the function for c

  # Note that beta is triangular, implying that the SEM is recursive (or acyclic)
  beta <- beta[topolOrd, topolOrd]
  lat <- c("Ucd")

  return(list(beta=beta, lat=lat))
}

# The SEM is constructed as
# Y = beta Y + eps, where
# beta is a matrix of the coefficients for all variables (V and U)
generateDatasetFromSEM <- function(beta, lat, N) {
  p <- ncol(beta)
  ident <- diag(1, p, p)
  colnames(ident) <- rownames(ident) <- colnames(beta)
  # Or, similarly,
  # Y = (I - Beta)^{-1} eps
  # For Gaussian Y and errors eps independent from each other (we have all the Us),
  # we just need to simulate Y from a multivariate distribution
  # of mean zero and cov Sigma=I
  IminusBinv <- ginv(ident - beta)
  Sigma = t(IminusBinv) %*% ident %*% IminusBinv

  valR <- matrixcalc::is.symmetric.matrix(Sigma) &&
    matrixcalc::is.positive.definite(Sigma, tol=1e-8)
  if (!valR) {
    stop("This SEM generates a non-positive definite covariance matrix. Try another SEM.")
  }

  dat <-  MASS::mvrnorm(N, rep(0, p), Sigma, empirical = FALSE)
  dat <- as.data.frame(dat)
  colnames(dat) <- colnames(beta)
  head(dat)

  uvars <- which(colnames(dat) %in% lat)
  dat <- dat[,-uvars]

  return(dat)
}

# Randomly generate a linear SEM following a dagitty DAG, adag,
# and then draw samples from it.
generateDatasetFromDAG <- function(adag, N, ntries=30) {
  done <- FALSE
  tries <- 0
  obs.dat <- NULL
  while (!done && tries <= ntries) {
    done <- tryCatch(
      {
        obs.dat <- dagitty::simulateSEM(adag, b.lower = -0.6, b.upper = 0.6, N=N)
        R <- cor(obs.dat)
        valR <- matrixcalc::is.symmetric.matrix(R) &&
          matrixcalc::is.positive.definite(R, tol=1e-8)
        valR
      }, error=function(cond) {
        message(cond)
        FALSE
      })
    tries <- tries + 1
  }
  return(obs.dat)
}

dagittyOracleCI <- function(x, y, S, suffStat) {
  g <- suffStat$g
  labels <- names(g)
  if (dagitty::dseparated(g, labels[x], labels[y], labels[S])) {
    return(1)
  } else {
    return(0)
  }
}

# receives a dagitty g of type "dag" or "mag" and returns
# the true PAG as an pcalg fci object
getTruePAG <- function(g, verbose = FALSE) {
  indepTest <- dagittyOracleCI
  if (graphType(g) == "dag") {
    g <- dagitty::toMAG(g)
  }
  suffStat <- list(g=g)
  truePag <- pcalg::fci(suffStat,
                        indepTest = indepTest,
                        labels= names(suffStat$g), alpha = 0.9999,
                        verbose = verbose)
  return(truePag)
}

igraph_from_graphNel <- function(graphN, latNodes){
  igraph_dag <- igraph::igraph.from.graphNEL(graphN, weight = FALSE)
  for (n in latNodes) {
    adj_list <- graphN@edgeL[[n]]$edges
    if (length(adj_list) == 2) {
      igraph_dag <- igraph::add_edges(igraph_dag, c(adj_list[1], adj_list[2], adj_list[2], adj_list[1]))
      igraph_dag <- igraph::set.edge.attribute(graph = igraph_dag,
                                               name ="description",
                                               index = c(length(igraph::E(igraph_dag))-1, length(igraph::E(igraph_dag))), value = "U")
    }
  }
  for (n in latNodes){
    igraph_dag <- igraph::delete_vertices(igraph_dag, n)
  }
  return(igraph_dag)
}


# A -> B -> C; B <- D <- Ucd -> C
getBDGraph <- function() {
  allvars <- c("A", "B", "C", "D", "Ucd")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["A","B"] <- 0; amat["B","A"] <- 1; # a -> b
  amat["B","C"] <- 0; amat["C","B"] <- 1; # b -> c
  amat["D","B"] <- 0; amat["B","D"] <- 1; # d -> b
  amat["Ucd","C"] <- 0; amat["C","Ucd"] <- 1; # Ucd -> c
  amat["Ucd","D"] <- 0; amat["D","Ucd"] <- 1; # Ucd -> c

  lat <- c("Ucd")
  adag <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(adag) <- lat

  return(list(adag=adag, amat=amat, lat=lat))
}

# DAG in Fig. 2b) at https://causalai.net/r42.pdf
getDAG2bR42 <- function(ret_dagg = TRUE) {
  allvars <- c("x1", "x2", "y1", "y2", "y3", "y4", "y5", "ux1x2", "uy2y3", "uy4y5", "ux1y3", "ux2y1")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["x2", "y2"] <- 0; amat["y2", "x2"] <- 1; # x2 -> y2
  amat["x1", "y1"] <- 0; amat["y1", "x1"] <- 1; # x1 -> y1
  amat["y4", "y3"] <- 0; amat["y3", "y4"] <- 1; # y4 -> y3
  amat["y5", "y1"] <- 0; amat["y1", "y5"] <- 1; # y5 -> y1
  amat["y5", "y4"] <- 0; amat["y4", "y5"] <- 1; # y5 -> y4
  amat["ux1x2", "x1"] <- 0; amat["x1", "ux1x2"] <- 1; # ux1x2 -> x1
  amat["ux1x2", "x2"] <- 0; amat["x2", "ux1x2"] <- 1; # ux1x2 -> x2
  amat["uy2y3", "y2"] <- 0; amat["y2", "uy2y3"] <- 1; # uy2y3 -> y2
  amat["uy2y3", "y3"] <- 0; amat["y3", "uy2y3"] <- 1; # uy2y3 -> y3
  amat["uy4y5", "y4"] <- 0; amat["y4", "uy4y5"] <- 1; # uy4y5 -> y4
  amat["uy4y5", "y5"] <- 0; amat["y5", "uy4y5"] <- 1; # uy4y5 -> y5
  amat["ux1y3", "x1"] <- 0; amat["x1", "ux1y3"] <- 1; # ux1y3 -> x1
  amat["ux1y3", "y3"] <- 0; amat["y3", "ux1y3"] <- 1; # ux1y3 -> y3
  amat["ux2y1", "x2"] <- 0; amat["x2", "ux2y1"] <- 1; # ux2y1 -> x2
  amat["ux2y1", "y1"] <- 0; amat["y1", "ux2y1"] <- 1; # ux2y1 -> y1
  # plot(as(t(amat), "graphNEL"))

  lat <- c("ux1x2", "uy2y3", "uy4y5", "ux1y3", "ux2y1")

  adag <- pcalg::pcalg2dagitty(amat, colnames(amat),type="dag")

  dagitty::latents(adag) <- lat
  return(list(adag=adag, amat=amat, lat=lat))
}

# returns an pcalg amat (adjacency matrix) of type amat.pag (same type for MAGs),
# where:
# 0: No edge
# 1: Circle
# 2: Arrowhead
# 3: Tail
dagitty2amat <- function(adagg, type="mag") {
  edg <- dagitty:::edges(adagg)
  node_names <- dagitty:::names.dagitty(adagg)
  ans_mat <- matrix(
    data = 0, nrow = length(node_names),
    ncol = length(node_names),
    dimnames = list(node_names, node_names)
  )

  diredg <- subset(edg, e == "->")

  ans_mat[as.matrix(diredg[c("w", "v")])] <- 3
  ans_mat[as.matrix(diredg[c("v", "w")])] <- 2

  bidiredg <-  subset(edg, e == "<->")
  ans_mat[as.matrix(bidiredg[c("w", "v")])] <- 2
  ans_mat[as.matrix(bidiredg[c("v", "w")])] <- 2

  return(ans_mat)
}

# type can be "pag" or "dag"
getMAG <- function(amat, type="pag") {
  if (type == "pag") {
    amat.mag <- pcalg::pag2magAM(amat, 1)
    #plotAG(amat.mag)
    magg <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")
    #plot(magg)
  } else {
    adagg <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
    magg <- dagitty::toMAG(adagg)
    amat.mag <- dagitty2amat(magg, type="mag")
  }
  return(list(amat.mag = amat.mag, magg=magg))
}
