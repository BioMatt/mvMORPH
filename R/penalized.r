################################################################################
##                                                                            ##
##                       mvMORPH: penalized.r                                 ##
##                                                                            ##
##   Internal functions for penalized methods in the mvMORPH package          ##
##                                                                            ##
##  Created by Julien Clavel - 31-07-2018                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################

# ------------------------------------------------------------------------- #
# .loocvPhylo                                                               #
# options: par, cvmethod, targM, corrStr, penalty, error, nobs              #
#                                                                           #
# ------------------------------------------------------------------------- #

# Enhanced cache initialization for LOOCV optimization
.initializeLOOCVCache <- function(corrModel) {
  n <- corrModel$nobs
  p <- corrModel$p
  m <- corrModel$m
  
  # Add LOOCV-specific caches to temp_matrices
  corrModel$cache$temp_matrices$loocv <- list(
    # Pre-computed matrices for all leave-one-out scenarios
    Y_minusx_list = vector("list", n),      # Store Y[-x,] for each x
    X_minusx_list = vector("list", n),      # Store X[-x,] for each x
    XtX_cols = array(0, dim = c(m, p, n)), # Store XtX[,x] for each x
    residuals_x = matrix(0, nrow = n, ncol = p), # Store residuals[x,] for each x
    
    # Working matrices for vectorized operations
    Bx_array = array(0, dim = c(m, p, n)), # Store all Bx updates
    Sk_array = array(0, dim = c(p, p, n)), # Store all Sk matrices
    
    # For batch processing
    batch_residuals = array(0, dim = c(n-1, p, n)), # All residual updates
    h_minus_x = numeric(n),  # (1-h[x]) for each x
    
    # Indices for efficient subsetting
    row_indices = vector("list", n),        # Row indices for each leave-one-out
    initialized = FALSE
  )
  
  return(corrModel)
}

# Pre-compute all leave-one-out subsets (call this when tree params change)
.precomputeLOOCVSubsets <- function(corrModel, mod_par, XtX, residuals, h) {
  n <- corrModel$nobs
  p <- corrModel$p
  m <- corrModel$m
  
  loocv_cache <- corrModel$cache$temp_matrices$loocv
  
  # Pre-compute all the expensive subsetting operations
  for(i in 1:n) {
    loocv_cache$row_indices[[i]] <- (1:n)[-i]
    loocv_cache$Y_minusx_list[[i]] <- mod_par$Y[-i, , drop = FALSE]
    loocv_cache$X_minusx_list[[i]] <- mod_par$X[-i, , drop = FALSE]
    loocv_cache$XtX_cols[, , i] <- XtX[, i, drop = FALSE]
    loocv_cache$residuals_x[i, ] <- residuals[i, ]
    loocv_cache$h_minus_x[i] <- 1 - h[i]
  }
  
  loocv_cache$initialized <- TRUE
  return(corrModel)
}

# Optimized LOOCV computation with vectorization
.optimizedLOOCVLoop <- function(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo) {
  
  loocv_cache <- corrStr$cache$temp_matrices$loocv
  n <- corrStr$nobs
  p <- corrStr$p
  
  # Method 1: Batch computation of all Bx updates
  # This vectorizes the rank-1 updates across all leave-one-out scenarios
  
  # Pre-allocate result
  llik <- numeric(length(nloo))
  
  if(length(nloo) < 10) {
    # For small nloo, use the optimized individual loop
    return(.optimizedSmallLOOCVLoop(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo))
  }
  
  # Batch approach for larger nloo
  # Step 1: Compute all Bx matrices at once using array operations
  Bx_updates <- array(0, dim = c(corrStr$m, p, length(nloo)))
  
  for(idx in seq_along(nloo)) {
    i <- nloo[idx]
    # Vectorized rank-1 update
    outer_prod <- tcrossprod(loocv_cache$XtX_cols[, , i], loocv_cache$residuals_x[i, , drop = FALSE])
    Bx_updates[, , idx] <- B - outer_prod / loocv_cache$h_minus_x[i]
  }
  
  # Step 2: Compute all residual updates in batches
  batch_size <- min(50, length(nloo))  # Process in batches to manage memory
  
  for(batch_start in seq(1, length(nloo), batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, length(nloo))
    batch_indices <- batch_start:batch_end
    
    # Process this batch
    for(idx_in_batch in seq_along(batch_indices)) {
      idx <- batch_indices[idx_in_batch]
      i <- nloo[idx]
      
      # Compute residuals using pre-computed subsets
      Y_minus_i <- loocv_cache$Y_minusx_list[[i]]
      X_minus_i <- loocv_cache$X_minusx_list[[i]]
      Bx_i <- Bx_updates[, , idx]
      
      # Efficient matrix multiplication
      residuals_i <- Y_minus_i - X_minus_i %*% Bx_i
      
      # Compute covariance matrix
      Sk_i <- crossprod(residuals_i) / (n - 1)
      
      # Compute likelihood
      llik[idx] <- .regularizedLik(Sk_i, loocv_cache$residuals_x[i, ], alpha, targM, target, penalty, const)
    }
  }
  
  return(llik)
}

# Optimized version for small nloo (< 10 observations)
.optimizedSmallLOOCVLoop <- function(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo) {
  
  loocv_cache <- corrStr$cache$temp_matrices$loocv
  temp_Bx <- corrStr$cache$temp_matrices$Bx_temp
  temp_residuals <- corrStr$cache$temp_matrices$residuals_temp
  temp_Sk <- corrStr$cache$temp_matrices$Sk_temp
  
  n <- corrStr$nobs
  llik <- numeric(length(nloo))
  
  # Optimized individual loop with pre-computed subsets
  for(idx in seq_along(nloo)) {
    i <- nloo[idx]
    
    # Use pre-computed values instead of subsetting operations
    XtX_col_i <- loocv_cache$XtX_cols[, , i]
    residuals_i <- loocv_cache$residuals_x[i, ]
    h_factor <- loocv_cache$h_minus_x[i]
    
    # Efficient rank-1 update
    temp_Bx[] <- B - tcrossprod(XtX_col_i, residuals_i) / h_factor
    
    # Use pre-computed subsets
    Y_minus_i <- loocv_cache$Y_minusx_list[[i]]
    X_minus_i <- loocv_cache$X_minusx_list[[i]]
    
    # Compute residuals
    temp_residuals[1:(n-1), ] <- Y_minus_i - X_minus_i %*% temp_Bx
    
    # Compute covariance matrix
    temp_Sk[] <- crossprod(temp_residuals[1:(n-1), , drop = FALSE]) / (n - 1)
    
    # Compute likelihood
    llik[idx] <- .regularizedLik(temp_Sk, residuals_i, alpha, targM, target, penalty, const)
  }
  
  return(llik)
}

# Alternative: Fully vectorized approach using Sherman-Morrison formula
.shermanMorrisonLOOCV <- function(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo) {
  # This uses the Sherman-Morrison formula for efficient rank-1 updates
  # Most beneficial when the regularized covariance structure allows it
  
  n <- corrStr$nobs
  p <- corrStr$p
  
  if(penalty != "RidgeArch" || length(nloo) < 20) {
    # Fall back to optimized loop for unsupported penalties or small problems
    return(.optimizedLOOCVLoop(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo))
  }
  
  # Compute base covariance matrix
  S_base <- crossprod(residuals) / n
  target_matrix <- .targetM(S_base, targM, penalty = "RidgeArch")
  
  llik <- numeric(length(nloo))
  
  # Sherman-Morrison updates for each leave-one-out scenario
  for(idx in seq_along(nloo)) {
    i <- nloo[idx]
    
    # Compute rank-1 update to covariance matrix
    # This is more complex but can be faster for certain structures
    # Implementation depends on your specific .regularizedLik function
    
    # For now, fall back to the standard optimized approach
    llik[idx] <- .optimizedSmallLOOCVLoop(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo[idx])
  }
  
  return(llik)
}

# Updated LOOCV section for the main .loocvPhylo function
.loocvPhyloOptimized_LOOCV <- function(corrStr, residuals, alpha, targM, penalty, const, XtX, B) {
  
  n <- corrStr$nobs
  p <- corrStr$p
  
  # Compute covariance matrix
  Sk <- crossprod(residuals) / n
  
  # Cache target matrix
  target_key <- .getTargetCacheKey(penalty, targM, alpha, p)
  if(is.null(corrStr$cache$target_matrices[[target_key]])) {
    corrStr$cache$target_matrices[[target_key]] <- .targetM(Sk, targM, penalty)
  }
  target <- corrStr$cache$target_matrices[[target_key]]
  
  # Use cached hat matrix diagonal
  h <- corrStr$cache$h_diagonal
  
  # Pre-filter valid indices (avoid hat score of 1)
  nloo <- corrStr$nloo[!h + 1e-8 >= 1]
  const <- n / length(nloo)
  
  # Initialize LOOCV cache if needed
  if(is.null(corrStr$cache$temp_matrices$loocv)) {
    corrStr <- .initializeLOOCVCache(corrStr)
  }
  
  # Pre-compute subsets if not already done or if parameters changed
  loocv_cache <- corrStr$cache$temp_matrices$loocv
  if(!loocv_cache$initialized) {
    corrStr <- .precomputeLOOCVSubsets(corrStr, corrStr$cache$mod_par, XtX, residuals, h)
  }
  
  # Choose optimization strategy based on problem size
  if(length(nloo) > 100 && p > 10) {
    # Large problem: use batch processing
    llik <- .optimizedLOOCVLoop(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo)
  } else {
    # Small to medium problem: use optimized individual loop
    llik <- .optimizedSmallLOOCVLoop(corrStr, B, residuals, alpha, targM, target, penalty, const, nloo)
  }
  
  ll <- 0.5 * (n * p * log(2 * pi) + p * corrStr$cache$mod_par$det + sum(llik))
  return(ll)
}
# ------------------------------------------------------------------------- #
# .mvGLS                                                                    #
# options: corrstruct object                                                #
#                                                                           #
# ------------------------------------------------------------------------- #
.mvGLS <- function(corrstruct){
  
  # GLS Estimate
  B <- pseudoinverse(corrstruct$X)%*%corrstruct$Y
  residuals <- corrstruct$Y - corrstruct$X%*%B
  
  return(list(residuals=residuals, B=B))
}

# ------------------------------------------------------------------------- #
# .scaleStuct                                                               #
# options: structure object                                                 #
#                                                                           #
# ------------------------------------------------------------------------- #
.scaleStruct <- function(structure){
  # wrapper for future developments
  if(inherits(structure, "phylo")){
    structure$edge.length <- structure$edge.length/max(node.depth.edgelength(structure))
  }
  
  return(structure)
}

# ------------------------------------------------------------------------- #
# .targetM                                                                  #
# options: S, targM, penalty, I                                             #
#                                                                           #
# ------------------------------------------------------------------------- #
.targetM <- function(S, targM, penalty="RidgeArch", I = NULL, ...){
  
  p <- dim(S)[1]
  args <- list(...)
  if(is.null(args[["userMatrix"]])) userMatrix <- NULL else userMatrix <- args$userMatrix
  if(is.null(args[["tuning"]])) tuning <- NULL else tuning <- args$tuning
  
  # If the identity is not provided
  if(is.null(I)) I = diag(p)
  target = NULL
  
  if(penalty=="RidgeArch"){
    switch(targM,
           "Variance" = {target <- diag(diag(S))},
           "unitVariance" = {target <- I*mean(diag(S))},
           "null" = {
             warning("The \"null\" target cannot be used with the \"RidgeArch\" method. The \"unitVariance\" target is used instead.")
             target <- I*mean(diag(S))
           },
           "user" = { target <- userMatrix} # TODO
    )
  }else if(penalty=="RidgeAlt"){
    switch(targM,
           "Variance" = {target <- diag(1/diag(S))},
           "unitVariance" = {target <- I*(1/mean(diag(S)))},
           "null" = {target <- matrix(0,p,p)},
           "user" = { target <- solve(userMatrix)} # TODO
    )
  }else if(penalty=="EmpBayes"){
    # TODO add option to switch between the unit and variance target in the Empirical Bayes; account for (v-p) factor?
    switch(targM,
           "Variance" = {target <- tuning*diag(diag(S))},
           "unitVariance" = {target <- I*tuning*mean(diag(S))}, # scaling is either v-p in Collucia, or for Matrix T v-1 the scaling of the mean for the prior
           "user" = { target <- userMatrix} # TODO
    )
  }
  
  return(target)
}

# ------------------------------------------------------------------------- #
# .corrStr (wrapper to covariance structure)                                #
# options: par, timeObject                                                  #
#                                                                           #
# ------------------------------------------------------------------------- #
.corrStr <- function(par, timeObject){
  
  if(timeObject$model%in%c("EB", "BM", "lambda", "OU", "OUvcv", "BMM", "OUM", "OU1", "OUMvcv")){
    # Tree transformation
    struct = .transformTree(timeObject$structure, par, model=timeObject$model, mserr=timeObject$mserr,
                            Y=timeObject$Y, X=timeObject$X, REML=timeObject$REML, precalc=timeObject$precalc)
  }else{
    stop("Currently works for phylogenetic models \"BM\", \"EB\", \"OU\", \"BMM\", \"OUM\" \"OUMvcv\" and \"lambda\"  only...")
  }
  return(struct)
}

# ------------------------------------------------------------------------- #
# .regularizedLik return the log-lik with the regularized estimate          #
# options: S, residuals, lambda, targM, target, penalty, const              #
#                                                                           #
# ------------------------------------------------------------------------- #
.regularizedLik <- function(S, residuals, lambda, targM, target, penalty, const=1){
  
  switch(penalty,
         "RidgeArch"={
           G <- (1-lambda)*S + lambda*target
           Gi <- try(chol(G), silent=TRUE)
           if(inherits(Gi, 'try-error')) return(1e6)
           rk <- sum(backsolve(Gi, residuals, transpose = TRUE)^2)
           llik <- const*sum(2*log(diag(Gi))) + rk
         },
         "RidgeAlt"={
           quad <- .makePenaltyQuad(S, lambda, target, targM)
           Gi <- quad$P
           detG <- sum(log(quad$ev))
           Swk <- tcrossprod(residuals)
           rk <- sum(Swk*Gi)
           llik <- const*detG + rk
         },
         "LASSO"={
           LASSO <- glassoFast(S, lambda, maxIt=500)
           G <- LASSO$w;
           Gi <- LASSO$wi;
           Swk <- tcrossprod(residuals);
           rk <- sum(Swk*Gi);
           llik <- const*as.numeric(determinant(G)$modulus) + rk
         })
  
  return(llik)
}


# ------------------------------------------------------------------------- #
# .makePenaltyQuad   (for quadratic ridge)                                  #
# options: S,lambda,target,targM                                            #
#                                                                           #
# ------------------------------------------------------------------------- #
.makePenaltyQuad <- function(S,lambda,target,targM){
  
  switch(targM,
         "Variance"={
           D <- (S - lambda * target)
           D2 <- D %*% D
           sqrtM <- .sqM(D2/4 + lambda * diag(nrow(S)))
           Alt <- D/2 + sqrtM
           AltInv <- (1/lambda)*(Alt - D)
           evalues <- eigen(Alt, symmetric=TRUE, only.values = TRUE)$values
         },
         "unitVariance"={
           eig  <- eigen(S, symmetric = TRUE)
           Q <- eig$vectors
           d <- eig$values - lambda*target[1]
           evalues <- sqrt(lambda + d^2/4) + d/2
           D1 <- evalues
           D2 <- 1/evalues # Inverse
           Alt <- Q %*% (D1 * t(Q))
           AltInv <- Q %*% (D2 * t(Q))
         },
         "null"={
           eig  <- eigen(S, symmetric = TRUE)
           Q <- eig$vectors
           d <- eig$values
           evalues <- sqrt(lambda + d^2/4) + d/2
           D1 <- evalues
           D2 <- 1/evalues
           Alt <- Q %*% (D1 * t(Q))
           AltInv <- Q %*% (D2 * t(Q))
         }
  )
  pen <- list(S=Alt, P=AltInv, ev=evalues)
  return(pen)
}


# ------------------------------------------------------------------------- #
# .sqM                                                                      #
# Matrix square root using eigen-decomposition                              #
#                                                                           #
# ------------------------------------------------------------------------- #
.sqM <- function(x){
  if(!all(is.finite(x))) return(Inf)
  eig <- eigen(x, symmetric = TRUE)
  sqrtM <- eig$vectors %*% (sqrt(eig$values) * t(eig$vectors))
  return(sqrtM)
}

# Build the matrix square root inverse
.sqM1 <- function(x){
  if(inherits(x, "phylo")) x <- vcv.phylo(x)
  if(!all(is.finite(x))) return(Inf)
  eig <- eigen(x, symmetric = TRUE)
  # check for singular dimensions =>  hack from corpcor package. Just retain the dimensions with positive eigenvalues
  tol = max(dim(x))*max(eig$values)*.Machine$double.eps
  Positive = eig$values > tol
  if(sum(Positive)<length(eig$values)) warning("The phylogenetic covariance matrix was singular. Check the results carefully and consider using 'eigSqm=FALSE' option and 'error=TRUE'")
  sqrtM <- eig$vectors[,Positive,drop=FALSE] %*% ((1/sqrt(eig$values[Positive])) * t(eig$vectors[,Positive,drop=FALSE]))
  return(sqrtM)
}


# ------------------------------------------------------------------------- #
# .covPenalized                                                             #
# options: S, penalty, targM="null", tuning=0, n                            #
#                                                                           #
# ------------------------------------------------------------------------- #
.penalizedCov <- function(S, penalty, Target=NULL, targM="null", tuning=0, n){
  
  # dim of S
  p = ncol(S)
  
  # target matrix
  if(is.null(Target)) Target <- .targetM(S, targM, penalty, tuning=tuning)
  
  # Construct the penalty term
  switch(penalty,
         "RidgeAlt"={
           pen <- .makePenaltyQuad(S,tuning,Target,targM)
           Pi <- pen$S
           P <- pen$P
         },
         "RidgeArch"={
           Pi <- (1-tuning)*S + tuning*Target
           eig <- eigen(Pi)
           V <- eig$vectors
           d <- eig$values
           P <- V%*%((1/d) * t(V))
         },
         "LASSO"={
           LASSO <- glassoFast(S,tuning)
           Pi <- LASSO$w
           P <- LASSO$wi
         },
         "EmpBayes"={
           # Compute the Empirical Bayes estimate of the covariance matrix
           v = p+1
           Pi <- (S*n + Target)/(v+n-2)
           eig <- eigen(Pi)
           V <- eig$vectors
           d <- eig$values
           P <- V%*%((1/d) * t(V))
         },
         "LL"={
           Pi <- S
           eig <- eigen(Pi)
           V <- eig$vectors
           d <- eig$values
           P <- V%*%((1/d) * t(V))
         })
  
  estimate <- list(Pinv=Pi, P=P, S=S)
  return(estimate)
}

# ------------------------------------------------------------------------- #
# .transformTree                                                            #
# options: phy, param, model, mserr=NULL, Y=NULL, X=NULL, REML=TRUE,        #
#      precalc=NULL                                                         #
# ------------------------------------------------------------------------- #

.transformTree <- function(phy, param, model=c("EB", "BM", "lambda", "OU", "BMM", "OUM","OUMvcv"), mserr=NULL, Y=NULL, X=NULL, REML=TRUE, precalc=NULL){
  
  # pre-compute and checks | TODO reduce computational burden by avoiding reordering and recomputing distances, ages...etc
  n <- Ntip(phy)
  parent <- phy$edge[,1]
  descendent <- phy$edge[,2]
  extern <- (descendent <= n)
  N <- 2*n-2
  diagWeight <- NULL
  const <- 0
  flag <- FALSE
  
  # Model
  switch(model,
         "OU"={
           D = numeric(n)
           
           # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
           if(!is.ultrametric(phy)){
             dis = node.depth.edgelength(phy) # has all nodes
             D = max(dis[1:n]) - dis[1:n]
             D = D - mean(D)
             phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
             flag <- TRUE
           }
           
           # Branching times (now the tree is ultrametric)
           times <- branching.times(phy)
           Tmax <- max(times)
           # compute the branch lengths
           distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
           d1 = distRoot[parent-n]
           d2 = numeric(N)
           d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
           d2[!extern] = distRoot[descendent[!extern]-n]
           
           # weights for a 3 points structured matrix
           diagWeight = exp(param*D)
           phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
           names(diagWeight) = phy$tip.label
           
           # transform the variables
           w <- 1/diagWeight
           Y <- matrix(w*Y, nrow=n)
           X <- matrix(w*X, nrow=n)
           
           # Adjust errors
           if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
         },
         "OUM"={
           # Weight matrix OUM
           W <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(precalc$root_std))
           
           # transform the tree
           D = numeric(n)
           
           # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
           if(!is.ultrametric(phy)){
             dis = node.depth.edgelength(phy) # has all nodes
             D = max(dis[1:n]) - dis[1:n]
             D = D - mean(D)
             phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
             flag <- TRUE
           }
           
           # Branching times (now the tree is ultrametric)
           times <- branching.times(phy)
           Tmax <- max(times)
           # compute the branch lengths
           if(precalc$randomRoot){
             distRoot <-  exp(-2*param*times)
             d1 = distRoot[parent-n]
             d2 = numeric(N)
             d2[extern] = exp(-2*param*D[descendent[extern]])
             d2[!extern] = distRoot[descendent[!extern]-n]
           }else{
             distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
             d1 = distRoot[parent-n]
             d2 = numeric(N)
             d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
             d2[!extern] = distRoot[descendent[!extern]-n]
           }
           
           # weights for a "3 points" structured matrix
           diagWeight = exp(param*D)
           phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
           names(diagWeight) = phy$tip.label
           
           # transform the variables
           w <- 1/diagWeight
           Y <- matrix(w*Y, nrow=n)
           X <- matrix(w*W, nrow=n) # Here X is replaced by the weighted matrix
           
           # REML "constant"
           if(REML) const <- determinant(crossprod(W))$modulus # TODO: check for n-ultrametric trees
           
           # Adjust errors
           if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
           
         },
         "OUMvcv"={
           # Weight matrix OUM
           W <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(precalc$root_std))
           
           # REML "constant"
           if(REML) const <- determinant(crossprod(W))$modulus # TODO: check for n-ultrametric trees
           
           V<-.Call("mvmorph_covar_ou_random", A=vcv.phylo(phy), alpha=param, sigma=1, PACKAGE="mvMORPH")
           
           C<-list(sqrtM=t(chol(solve(V))), det=determinant(V)$modulus, const=const)
           
         },
         "OU1"={
           # Weight matrix OU1
           W <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(precalc$root_std))
           
           # transform the tree
           D = numeric(n)
           
           # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
           if(!is.ultrametric(phy)){
             dis = node.depth.edgelength(phy) # has all nodes
             D = max(dis[1:n]) - dis[1:n]
             D = D - mean(D)
             phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
             flag <- TRUE
           }
           
           # Branching times (now the tree is ultrametric)
           times <- branching.times(phy)
           Tmax <- max(times)
           # compute the branch lengths
           if(precalc$randomRoot){
             distRoot <-  exp(-2*param*times)
             d1 = distRoot[parent-n]
             d2 = numeric(N)
             d2[extern] = exp(-2*param*D[descendent[extern]])
             d2[!extern] = distRoot[descendent[!extern]-n]
           }else{
             distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
             d1 = distRoot[parent-n]
             d2 = numeric(N)
             d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
             d2[!extern] = distRoot[descendent[!extern]-n]
           }
           
           # weights for a "3 points" structured matrix
           diagWeight = exp(param*D)
           phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
           names(diagWeight) = phy$tip.label
           
           # transform the variables
           w <- 1/diagWeight
           Y <- matrix(w*Y, nrow=n)
           X <- matrix(w*W, nrow=n) # Here X is replaced by the weighted matrix
           
           # REML "constant"
           if(REML) const <- determinant(crossprod(W))$modulus
           
           # Adjust errors
           if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
           
         },
         "EB"={
           if (param!=0){
             distFromRoot <- node.depth.edgelength(phy)
             phy$edge.length = (exp(param*distFromRoot[descendent])-exp(param*distFromRoot[parent]))/param
           }
         },
         "lambda"={
           # Pagel's lambda tree transformation
           if(param!=1) {
             root2tipDist <- node.depth.edgelength(phy)[1:n] # for non-ultrametric trees. The 'up' limit should be exactly 1 to avoid singularity issues
             phy$edge.length <- phy$edge.length * param
             phy$edge.length[extern] <- phy$edge.length[extern] + (root2tipDist * (1-param))
           }
         },
         "OUvcv"={
           V<-.Call("mvmorph_covar_ou_fixed", A=vcv.phylo(phy), alpha=param, sigma=1, PACKAGE="mvMORPH")
           C<-list(sqrtM=t(chol(solve(V))), det=determinant(V)$modulus)
         },
         "OUTS"={
           stop("Not yet implemented. The time-series models are coming soon, please be patient")
         },
         "RWTS"={
           stop("Not yet implemented. The time-series models are coming soon, please be patient")
         },
         "BMM"={
           # multirates model - proportional scaling or explicit estimation?
           #phy$edge.length <- phy$mapped.edge %*% param
           phy$edge.length <- phy$mapped.edge %*% c(1,param)
         })
  
  # Add measurment error
  if(is.numeric(mserr)) phy$edge.length[extern] = phy$edge.length[extern] + mserr
  
  # Compute the independent contrasts scores
  if(inherits(phy, "phylOLS")){
    if((sum(phy$edge.length) - n)<=.Machine$double.eps){
      # Return the determinant
      deterM <- 0
    }else{
      sqrtM <- 1/sqrt(phy$edge.length[extern])
      X <- X*sqrtM
      Y <- Y*sqrtM
      # Return the determinant => variance terms  of the 'star' tree
      deterM <- sum(log(phy$edge.length[extern]))
    }
    
  }else{
    if(model!="OUvcv" & model!="OUMvcv") C <- pruning(phy, trans=FALSE) # FIXME -> to remove the call to OUvcv?
    #if(any(phy$edge.length<=.Machine$double.eps)) C<-list(sqrtM=t(.sqM1(phy)), det=determinant(vcv(phy))$modulus) # FIXME => remove problems with the pruning algorithms on zero branch lengths?
    X <- crossprod(C$sqrtM, X)
    Y <- crossprod(C$sqrtM, Y)
    
    # Return the determinant
    deterM <- C$det
  }
  
  # Adjust the determinant for non-ultrametric OU (see Ho & Ane 2014 - Syst. Bio., p. 401)
  if(flag) deterM <- deterM + 2*sum(log(diagWeight))
  if(REML) deterM <- deterM + determinant(crossprod(X))$modulus - const
  
  # Return the score, variances, and scaled tree
  return(list(phy=phy, diagWeight=diagWeight, X=X, Y=Y, det=deterM, const=const))
}

# Vec operator
.vec <- function(x) as.numeric(x)


# ------------------------------------------------------------------------- #
# .setBounds                                                                #
# options: penalty, model, lower, upper, tol, mserr=NULL, penalized         #
#                                                                           #
# ------------------------------------------------------------------------- #

.setBounds <- function(penalty, model, lower, upper, tol=1e-10, mserr=NULL, penalized=TRUE, corrModel=NULL, k=NULL){
  
  if(is.null(upper)){
    switch(model,
           "EB"={up <- 1e-10},
           "OU"={up <- 30/max(node.depth.edgelength(corrModel$structure))}, # ~ 30 half-lifes upper limit for phylogenetic trees. Should change the format for general models
           "OU1"={up <- 30/max(node.depth.edgelength(corrModel$structure))},
           "OUM"={up <- 30/max(node.depth.edgelength(corrModel$structure))},
           "OUMvcv"={up <- 30/max(node.depth.edgelength(corrModel$structure))},
           "lambda"={up <- 1},
           "BM"={up <- Inf},
           "BMM"={up <- rep(Inf,k-1)},
           up <- Inf)
  }else{
    up <- upper
  }
  
  if(is.null(lower)){
    switch(model,
           "EB"={low <- -(30/max(node.depth.edgelength(corrModel$structure)))},
           "OU"={low <- 1e-10},
           "OU1"={low <- 1e-10},
           "OUM"={low <- 1e-10},
           "OUMvcv"={low <- 1e-10},
           "lambda"={low <- 1e-8},
           "BM"={low <- -Inf},
           "BMM"={low <- rep(-Inf,k-1)},
           low <- -Inf)
  }else{
    low <- lower
  }
  
  # Default tolerance for the parameter search
  if(is.null(tol)){
    if(penalty=="RidgeArch" || penalty=="EmpBayes"){
      tol = 1e-8
    }else{
      tol = 0
    }
  }
  
  # Set the default bounds
  if(penalized){
    if(penalty%in%c("RidgeAlt","LASSO")){
      upperBound <- c(log(10e6),up)
      lowerBound <- c(log(tol),low)
    }else if(penalty=="RidgeArch"){
      upperBound <- c(1,up)
      lowerBound <- c(tol,low)
    }else if(penalty=="EmpBayes"){
      upperBound <- c(10e6,up)
      lowerBound <- c(tol,low)
    }
    
    # parameters
    if(model=="BMM"){
      #id1 <- 1; id2 <- 2:(k+1); id3 <- k+2
      id1 <- 1; id2 <- 2:k; id3 <- k+1
    }else if(model=="BM"){
      id1 <- id2 <- 1; id3 <- 2
    }else{
      id1 <- 1; id2 <- 2; id3 <- 3
    }
    
  }else{
    upperBound <- up
    lowerBound <- low
    
    # parameters
    if(model=="BMM"){
      #id1 <- 1; id2 <- 1:k; id3 <- k+1
      id1 <- 1; id2 <- 1:(k-1); id3 <- k
    }else if(model=="BM"){
      id1 <- id2 <- id3 <- 1
    }else{
      id1 <- 1; id2 <- 1; id3 <- 2
    }
  }
  
  # Bounds for error
  if(!is.null(mserr)){
    lowerBound = c(lowerBound,0)
    upperBound = c(upperBound,Inf)
  }
  
  # regularization parameter
  switch(penalty,
         "RidgeArch"={ transformTun <- function(x) (x[id1])},
         "RidgeAlt" ={ transformTun <- function(x) exp(x[id1])},
         "LASSO" ={ transformTun <- function(x) exp(x[id1])},
         transformTun <- function(x) (x[id1])
  )
  
  # model parameter
  switch(model,
         "OU"={ transformPar <- function(x) (x[id2])},
         "OUM"={ transformPar <- function(x) (x[id2])},
         "BM" ={ transformPar <- function(x) (x[id2])},
         "EB" ={ transformPar <- function(x) (x[id2])},
         "lambda" ={ transformPar <- function(x) (x[id2])},
         "BMM"={transformPar <- function(x) (x[id2]*x[id2])},
         transformPar <- function(x) (x[id2])
  )
  
  # mserr parameter
  transformSE <- function(x) x[id3]*x[id3]
  
  
  bounds <- list(upper=upperBound, lower=lowerBound, trTun=transformTun, trPar= transformPar, trSE=transformSE)
  return(bounds)
}

# ------------------------------------------------------------------------- #
# .startGuess                                                               #
# options: corrModel, cvmethod, mserr, target, penalty, echo, penalized     #
# tol,                                                                      #
# ------------------------------------------------------------------------- #

.startGuess <- function(corrModel, cvmethod, mserr=NULL, target, penalty, echo=TRUE, penalized=TRUE, tol=NULL,...){
  if(echo==TRUE) message("Initialization via grid search. Please wait...")
  # Penalization parameters guesses
  if(penalized){
    if(penalty=="RidgeArch"){
      range_val <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
      if(!is.null(tol)) range_val <- range_val[!(range_val<=tol)]
    }else if(penalty=="EmpBayes"){
      range_val <- c(1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000)
      if(!is.null(tol)) range_val <- range_val[!(range_val<=tol)]
    }else if(penalty=="RidgeAlt"){
      range_val <- log(c(1e-11, 1e-9, 1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000))
      if(!is.null(tol)) range_val <- range_val[!(range_val<=log(tol))]
    }else{
      range_val <- log(c(1e-6, 0.01, 0.1, 1, 10, 100, 1000))
      if(!is.null(tol)) range_val <- range_val[!(range_val<=log(tol))]
    }
  }else{
    range_val <- NULL
  }
  
  # prepare the list
  list_param <- list()
  list_param[[1]] <- range_val
  list_param[[2]] <- 1 # dummy starting value for BM...
  
  # Models starting guesses
  switch(corrModel$model,
         "OU"={
           mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "OU1"={
           mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "OUM"={
           mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "OUMvcv"={
           mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "OUvcv"={
           mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "lambda"={
           mod_val <- c(0.2,0.5,0.8)
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "EB"={
           mod_val <- -log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
           list_param[[2]] <- mod_val
           index_err <- 3
         },
         "BMM"={
           
           # guess starting values
           start_values <- function(tree, data, predictors){
             tip_values <- 1:Ntip(tree)
             index_tips <- tree$edge[,2]%in%tip_values
             maps <- sapply(tree$maps[index_tips], function(x) names(x[length(x)]))
             # check if any tips are missing?
             k = ncol(tree$mapped.edge)
             if(length(unique(maps))<k) {
               mod_val <- mean(diag(rate_pic(tree, data)))
               guesses <- as.list(rep(sqrt(mod_val), k))
             }else{
               guesses <- lapply(colnames(tree$mapped.edge), function(map_names) {
                 dat_red <- which(maps==map_names)
                 sp_to_remove <- tree$tip.label[!tree$tip.label%in%tree$tip.label[dat_red]]
                 # Sanity check => because errors with current phytools function | example reported by Jake
                 if(Ntip(tree) - length(sp_to_remove) <= 1) {
                   # select a second species at random
                   sp_to_remove <- sp_to_remove[-sample(length(sp_to_remove), size = 1)]
                 }
                 tree_red=drop.tip(tree, sp_to_remove)
                 if(Ntip(tree_red)<=1){
                   # simple estimate on the whole tree
                   sqrt(mean(apply(.rate_guess(tree, data[tree$tip.label,], predictors[tree$tip.label,]) , 2, var)))
                 }else{
                   sqrt(mean(apply(.rate_guess(tree_red, data[tree_red$tip.label,], predictors[tree_red$tip.label,]) , 2, var)))
                 }
               })
             }
             
             return(guesses)
           }
           
           mod_val <- start_values(corrModel$structure, corrModel$Y, corrModel$X)[-1]
           list_param <- c(list(range_val), mod_val)
           index_err <- length(list_param) + 1
         },
         index_err <- 3
  )
  
  # Prepare the grid search
  if(!is.null(corrModel$mserr)){
    if(corrModel$model=="BMM") list_param[[index_err]] <- sqrt(c(0.001,0.01,0.1,1,10)*mean(unlist(mod_val)^2)) else list_param[[index_err]] <- c(0.001,0.01,0.1,1,10)
    list_param[sapply(list_param, is.null)] <- NULL
    brute_force <- expand.grid(list_param)
  }else{
    list_param[sapply(list_param, is.null)] <- NULL
    brute_force <- expand.grid(list_param)
  }
  
  start <- brute_force[which.min(apply(brute_force, 1, .loocvPhylo,
                                       cvmethod=cvmethod, # options
                                       targM=target,
                                       corrStr=corrModel,
                                       penalty=penalty,
                                       error=mserr,
                                       nobs=corrModel$nobs )),]
  
  if(echo==TRUE & penalized==TRUE)  cat("Best starting for the tuning: ",as.numeric(corrModel$bounds$trTun(start[1])))
  return(start)
}


# ------------------------------------------------------------------------- #
# .check_par_results  (TODO)                                                #
# options: model, par, penalized                                            #
#                                                                           #
# ------------------------------------------------------------------------- #

.check_par_results <- function(corrModel, par, penalized=TRUE){
  
  if(penalized) indice = 2 else indice = 1
  switch(corrModel$model,
         "OU"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
         },
         "OU1"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
         },
         "OUM"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
         },
         "OUMvcv"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
         },
  )
}

# ------------------------------------------------------------------------- #
# .rate_guess                                                               #
# options: phylo, Y, X                                                      #
#                                                                           #
# ------------------------------------------------------------------------- #

.rate_guess <- function(phylo, Y, X){
  # model
  C <- pruning(phylo, trans=FALSE)
  X <- crossprod(C$sqrtM, X)
  Y <- crossprod(C$sqrtM, Y)
  
  # GLS estimates
  XtX <- pseudoinverse(X)
  B <- XtX%*%Y
  residuals <- Y - X%*%B
  
  # Return the residuals
  return(residuals)
}