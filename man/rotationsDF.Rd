\name{rotationsDF}
    \alias{rotationsDF}
	\alias{cubimax.df}
	\alias{fssQ.df}
	\alias{fssT.df}

\title{Rotations}
\usage{
	cubimax.df(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000)
	fssQ.df(A, Tmat=diag(ncol(A)), kij=2, normalize=FALSE, eps=1e-5, maxit=1000)
	fssT.df(A, Tmat=diag(ncol(A)), kij=2, normalize=FALSE, eps=1e-5, maxit=1000)
}
\arguments{
    \item{A}{an initial factor loadings matrix to be rotated.}
    \item{Tmat}{initial rotation matrix.}
    \item{kij}{minimum additional number of forced simple structure loadings in a pair of factors.}
    \item{normalize}{parameter passed to optimization routine (\code{GPForth.df} or \code{GPFoblq.df}).}
    \item{eps}{parameter passed to optimization routine (\code{GPForth.df} or \code{GPFoblq.df}).}
    \item{maxit}{parameter passed to optimization routine (\code{GPForth.df} or \code{GPFoblq.df}).}
}
\value{A list (which includes elements used by \code{factanal}) with:
   \item{loadings}{Lh from \code{GPForth.df} or \code{GPFoblq.df}.} 
   \item{Th}{Th from \code{GPForth.df} or \code{GPFoblq.df}.}  
   \item{Table}{Table from \code{GPForth.df} or \code{GPFoblq.df}.}
   \item{method}{A string indicating the rotation objective function.}
   \item{orthogonal}{A logical indicating if the rotation is orthogonal.}
   \item{convergence}{Convergence indicator from \code{GPForth.df} or \code{GPFoblq.df}.}
   \item{Phi}{t(Th) \%*\% Th.  The covariance matrix of the rotated factors.
   This will be the identity matrix for orthogonal
   rotations so is omitted (NULL) for the result from GPForth.df.}
   }
\description{
   Optimize factor loading rotation objective.
}
\details{
 The functions listed here optimize a rotation objective. They can be used directly or the
 function name can be passed to factor analysis functions like \code{factanal}.  
 
 Available rotations are 
\tabular{lll}{
 \code{cubimax.df}   \tab orthogonal  \tab  Cubimax	\cr
 \code{fssQ.df} 	  \tab oblique	   \tab  Forced Simple Structure (see Vignette)   \cr
 \code{fssT.df} 	  \tab orthogonal  \tab  Forced Simple Structure (see Vignette)   \cr
 }
 The argument \code{kij} for Forced Simple Structure is the minimum number 
 of forced simple structure loadings in a pair of factors, in addition to 
 the number of factors itself. Meaningful values are integers (1, ..., items - factors )

 }
\examples{
  data(ability.cov)
  x <- factanal(factors = 3, covmat = ability.cov, rotation="none")
  fssT.df(x$loadings, kij = 2)
  fssQ.df(x$loadings, kij = 4)

  # 3 different methods
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <- factanal(factors = 3, covmat=NetherlandsTV, rotation="none")
  #
  fa.varimax <- GPForth.df(loadings(fa.unrotated), method = "varimax", normalize = TRUE)
  fa.cubimax <- cubimax.df(loadings(fa.unrotated), normalize = TRUE)
  fa.quartimax <- GPForth.df(loadings(fa.unrotated), method = "quartimax", normalize = TRUE)
  print(cbind(loadings(fa.varimax), loadings(fa.cubimax), loadings(fa.quartimax)), digits = 2)   
  }
\seealso{
   \code{\link{GPForth.df}}, 
   \code{\link{GPFoblq.df}}, 
   \code{\link{ff.cubimax}}, 
   \code{\link{ff.fss}}, 
   \code{\link[stats]{factanal}} 
 }
\references{
  Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

    Jennrich, R.I. (2004) Derivative free gradient projection algorithms for rotation,
    \emph{Psychometrika}, \bold{69}(3), 475--480.
  }

\author{Coen A. Bernaards and Robert I. Jennrich}
\concept{rotation}
\keyword{multivariate}


