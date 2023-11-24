
# GPForth.df is the main GP algorithm for orthogonal rotation.
# GPFoblq.df is the main GP algorithm for oblique rotation.
# Gf computes the numerical derivative and is needed by GPForth.df and GPFoblq.df.
# For both algorithms is required: a loadings matrix A. Optional  
# a initial rotation matrix Tmat. By default this is the identity matrix.
# Optional: the rotation method to be used. Between quation marks have to
# be the last part of the name of the ff function, e.g. for ff.varimax
# the argument is "varimax". Identical arguments can be used for oblique
# rotation. Some rotation criteria (including simplimax, pst, fss, 
# cf,...) require one or more additional arguments. 


GPForth.df <- function(A, Tmat = diag(ncol(A)), normalize = FALSE, eps = 1e-05, 
  maxit = 1000, method="varimax", methodArgs = NULL){
# begin of function Gf
  Gf <- function(Tmat, A, Method, methodArgs){
   k <- nrow(Tmat)
   G <- Z <- matrix(0,k,k)
   ep <- 0.0001
   for (r in 1:k){
     for (s in 1:k){
        dT <- Z
        dT[r,s] <- ep
        p1 <- do.call(Method, append(list(A %*% (Tmat + dT)), methodArgs))
        p2 <- do.call(Method, append(list(A %*% (Tmat - dT)), methodArgs))
        G[r,s] <- (p1$f - p2$f)/(2 * ep)
      }
    }
  G
  }
# end of function Gf
 if((!is.logical(normalize)) || normalize) {
   W <- NormalizingWeight(A, normalize=normalize)
   normalize <- TRUE
   A <- A/W
  }
  k <- nrow(Tmat)
  if (1 >= k) 
    stop("rotation does not make sense for single factor models.")
 Method <- paste("ff", method, sep = ".")
 al <- 1
 Table <- NULL
 L <- A %*% Tmat
 for (iter in 0:maxit){
   f <- do.call(Method, append(list(L), methodArgs))
   G <- Gf(Tmat, A, Method, methodArgs)
   M <- crossprod(Tmat,G)
   S <- (M + t(M))/2
   Gp <- G-Tmat %*% S
      s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table,c(iter,f$f,log10(s),al))
   if (s < eps)
         break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al * Gp
     UDV <- svd(X)
     Tmatt <- UDV$u %*% t(UDV$v)
     L <- A %*% Tmatt
     ft <- do.call(Method, append(list(L), methodArgs))
     if (ft$f < (f$f-.5*s^2*al))
       break
     al <- al/2
    }
    Tmat <- Tmatt
  }
  L <- A %*% Tmat
  convergence <- (s < eps)
  if ((iter == maxit) & !convergence) 
        warning("convergence not obtained in GPForth.df. ", maxit, 
            " iterations used.")
  if(normalize) L <- L * W
  dimnames(L) <- dimnames(A)
  r <- list(loadings = L, Th = Tmat, Table = Table, method = f$Method, 
        orthogonal = TRUE, convergence = convergence, G = G)
  colnames(r$Table) <- c("iter", "f", "log10(s)", "alpha")
  class(r) <- "GPArotation"
  r
}

 

GPFoblq.df <- function (A, Tmat = diag(ncol(A)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000, method = "quartimin", methodArgs = NULL){ 
  Gf <- function(Tmat, A, Method, methodArgs){ # begin function Gf
   k <- nrow(Tmat)
   ep <- .0001
   G <- Z <- matrix(0,k,k)
   for (r in 1:k){
     for (s in 1:k){
        dT <- Z
        dT[r,s] <- ep
        p1 <- do.call(Method, append(list(A %*% t(solve(Tmat + dT))), methodArgs))
        p2 <- do.call(Method, append(list(A %*% t(solve(Tmat - dT))), methodArgs))
        G[r,s] <- (p1$f - p2$f)/(2*ep)
      }
    }
  G
  }  # end of function Gf
  if((!is.logical(normalize)) || normalize) {
    W <- NormalizingWeight(A, normalize=normalize)
    normalize <- TRUE
    A <- A/W
  }
  k <- nrow(Tmat)
  if (1 >= k) 
     stop("rotation does not make sense for single factor models.")
 Method <- paste("ff", method, sep = ".")
 al <- 1
 Table <- NULL
 L <- A %*% t(solve(Tmat))
 for (iter in 0:maxit){
   f <- do.call(Method, append(list(L), methodArgs))
   G <- Gf(Tmat, A, Method, methodArgs)
    Gp <- G-Tmat %*% diag(apply(Tmat*G,2,sum))
    s <- sqrt(sum(diag(crossprod(Gp))))
    Table <- rbind(Table,c(iter,f$f,log10(s),al))
    if (s < eps)
       break
     al <- 2*al
     for (i in 0:10){
       X <- Tmat-al*Gp
       v <- 1/sqrt(apply(X^2,2,sum))
       Tmatt <- X %*% diag(v)
       L <- A %*% t(solve(Tmatt))
       ft <- do.call(Method, append(list(L), methodArgs))
       if (ft$f < (f$f-.5*s^2*al))
          break
       al <- al/2
     } 
     Tmat <- Tmatt
   }
   L <- A %*% t(solve(Tmat))
   convergence <- (s < eps)
   if ((iter == maxit) & !convergence) 
        warning("convergence not obtained in GPFoblq.df. ", maxit, 
            " iterations used.")
   if(normalize) L <- L * W
    dimnames(L) <- dimnames(A)
    r <- list(loadings = L, Phi = t(Tmat) %*% Tmat, Th = Tmat, 
        Table = Table, method = ft$Method, orthogonal = FALSE, 
        convergence = convergence, G = G)
    colnames(r$Table) <- c("iter", "f", "log10(s)", "alpha")
    class(r) <- "GPArotation"
    r
}

ff.quartimax <- function(L){
  f = -sum(L^4) / 4
  list(f = f, Method = "DF-Quartimax")
}

ff.varimax <- function(L){
  QL <- sweep(L^2,2,apply(L^2,2,mean),"-")
  f= -sqrt(sum(diag(crossprod(QL))))^2/4
  list(f = f,
  		Method = "DF-Varimax")
}

ff.quartimin <- function(L){
  L2 <- L^2
  k <- ncol(L)
  M <- matrix(1,k,k)-diag(k)
  f <- sum(L2 * (L2 %*% M))/4
  list(f = f,
  		Method = "DF-Quartimin")
} 

ff.target <- function(L,Target = NULL){
  # Needs Target matrix, e.g.  Target <-matrix(c(rep(9,4),rep(0,8),rep(9,4)),8) 
  f <- sum((L-Target)^2, na.rm = TRUE)
  list(f = f,
  		Method = "DF-Target rotation")
}

ff.pst <- function(L,W = NULL,Target = NULL){
  # Needs weight matrix W with 1's at specified values, 0 otherwise
  # e.g. W = matrix(c(rep(1,4),rep(0,8),rep(1,4)),8). 
  # When W has only 1's this is procrustes rotation
  # Needs a Target matrix Target with hypothesized factor loadings.
  # e.g. Target = matrix(0,8,2)
  Btilde <- W * Target
  f <- sum((W*L-Btilde)^2, na.rm = TRUE)
  list(f = f,
  		Method = "DF-PST")
}

ff.oblimax <- function(L){
  f <- -(log(sum(L^4))-2*log(sum(L^2)))
  list(f = f,
  		Method = "DF-Oblimax")
}

ff.entropy <- function(L){
  f <- -sum(L^2 * log(L^2 + (L^2==0)))/2
  list(f = f,
  		Method = "DF-Entropy")
}
 
cubimax.df <-  function(A, Tmat = diag(ncol(A)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000){
    GPForth.df(A, Tmat = Tmat, method = "cubimax", normalize = normalize, 
        eps = eps, maxit = maxit)
}
 
ff.cubimax <- function(L){
  f <- -sum(diag(t(L^2) %*% abs(L)))
  list(f = f,
  		Method = "DF-Cubimax")
}

ff.simplimax <- function(L,k=nrow(L)){
  # k: Number of close to zero loadings
  Imat <- sign(L^2 <= sort(L^2)[k])
  f <- sum(Imat*L^2)
  list(f = f,
  		Method = "DF-Simplimax")
}


fssT.df <-  function(A, Tmat = diag(ncol(A)), kij = 2, normalize = FALSE, eps = 1e-05, 
    maxit = 1000){
    GPForth.df(A, Tmat = Tmat, method = "fss", normalize = normalize, 
        eps = eps, maxit = maxit, methodArgs = (kij = kij) )
}
fssQ.df <-  function(A, Tmat = diag(ncol(A)), kij = 2, normalize = FALSE, eps = 1e-05, 
    maxit = 1000){
    GPFoblq.df(A, Tmat = Tmat, method = "fss", normalize = normalize, 
        eps = eps, maxit = maxit, methodArgs = (kij = kij) )
}

ff.fss <- function(L, kij=2){
  m <- ncol(L)
  p <- nrow(L)
  zm <- m + kij
  Imat <- matrix(0, p, m)
  for (j in 1:m){
    Imat[abs(L[,j]) <= sort(abs(L[,j]))[zm],j] <- 1 }
  for (i in 1:(m-1)){
    for (j in (i+1):m){
      nz <- sum( (Imat[,i] + Imat[,j]) ==1)
      while (nz < zm && sum(Imat[ ,c(i,j)]) < m * 2){
	    tbc <- c(abs(L[,i]), abs(L[,j]))
	    tbcs <- sort(tbc [c(Imat[,i], Imat[,j])==0])[1]
	    Imat[abs(L) == tbcs] <- 1
	    nz <- sum( (Imat[,i] + Imat[,j]) ==1)
      }
    }
  }
  Method <- paste("DF-Forced Simple Structure (kij = ",kij,")", sep="")
  f <- sum(Imat*L^2)
  list(f = f, Imat = Imat,
       Method = Method)
}


ff.bentler <- function(L){
  L2 <- L^2
  M <- crossprod(L2)
  D <- diag(diag(M))
  f <- -(log(det(M))-log(det(D)))/4
  list(f = f,
  		Method = "DF-Bentler")
}

ff.geomin <- function(L,delta=0.01){
  k <- ncol(L)
  L2 <- L^2 + delta
  pro <- exp(apply(log(L2),1,sum)/k) #apply(L2,1,prod)^(1/k)
  f <- sum(pro)
  list(f = f,
  		Method = "DF-Geomin")
}

ff.infomax <- function(L){
  k <- ncol(L)
  S <- L^2
  s <- sum(S)
  s1 <- apply(S, 1, sum)
  s2 <- apply(S, 2, sum)
  E <- S/s
  e1 <- s1/s
  e2 <- s2/s
  Q0 <- sum(-E * log(E))
  Q1 <- sum(-e1 * log(e1))
  Q2 <- sum(-e2 * log(e2))
  f <- log(k) + Q0 - Q1 - Q2
  list(f = f,
  		Method = "DF-Infomax")
}


ff.cf <- function(L,kappa=0){
  k <- ncol(L)
  p <- nrow(L)
  # kappa <- 0 # Quartimax 
  # kappa <- 1/p # Varimax
  # kappa <- m/(2*p) # Equamax
  # kappa <- (m-1)/(p+m-2) # Parsimax
  # kappa <- 1 # Factor parsimony
  # Method <- paste("Crawford-Ferguson:k=",kappa,sep="")
  N <- matrix(1,k,k)-diag(k)
  M <- matrix(1,p,p)-diag(p)
  L2 <- L^2
  f1 <- (1-kappa)*sum(diag(crossprod(L2,L2 %*% N)))/4
  f2 <- kappa*sum(diag(crossprod(L2,M %*% L2)))/4
  f <- f1 + f2
  list(f = f,
  		Method = "DF-Crawford-Ferguson")
}
