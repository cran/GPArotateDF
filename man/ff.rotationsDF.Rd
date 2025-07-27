\name{ff.rotationsDF}
	\alias{ff.bentler}
	\alias{ff.cf}
	\alias{ff.cubimax} 
	\alias{ff.entropy}
	\alias{ff.fss} 
	\alias{ff.geomin} 
	\alias{ff.infomax} 
	\alias{ff.oblimax}
	\alias{ff.pst}
	\alias{ff.quartimax}
	\alias{ff.quartimin} 
	\alias{ff.simplimax} 
	\alias{ff.target} 
	\alias{ff.varimax} 

\title{Rotations}
\usage{
	ff.bentler(L)
	ff.cf(L, kappa=0)
	ff.cubimax(L)
	ff.entropy(L)
	ff.geomin(L, delta=0.01)
	ff.infomax(L) 
	ff.oblimax(L)
	ff.pst(L, W=NULL, Target=NULL)
	ff.quartimax(L)
	ff.quartimin(L)
	ff.simplimax(L, k=nrow(L))
	ff.fss(L, kij=2)
	ff.target(L, Target=NULL) 
	ff.varimax(L)
}
\arguments{
    \item{L}{a factor loading matrix}
    \item{kappa}{see details.}
    \item{delta}{constant added to Lambda^2 in objective calculation.}
    \item{Target}{rotation target for objective calculation.}
    \item{W}{weighting of each element in target.}
    \item{k}{number of close to zero loadings.}
    \item{kij}{minimum additional number of forced simple structure loadings in a pair of factors.}
}

\value{
   \item{f}{criterion function value.}
   \item{method}{A string indicating the rotation objective function.}
   }
\description{
   Optimize factor loading rotation objective.
}
\details{
 These functions are used to optimize a rotation objective. The name need to be included
 in a call to \code{GPForth.df} or \code{GPFoblq}. Calling the functions itself computes the values
 but no rotation is performed.
 
 Functions listed here are all exported through NAMESPACE, and primarily serve as examples 
 for programming new rotation methods.  New rotation methods can be programmed with a name
 \code{ff.newmethod}. The inputs are the matrix L, and optionally any additional arguments. The
 output should be a list with elements 
\tabular{ll}{
  \code{f}       \tab the value of the criterion at L.\cr
  \code{Method}  \tab a string indicating the criterion.\cr         
  }
  Please note that the function value \code{f} has to be minimized. If the rotation criterion
  is supposed to maximize, then use the negative of the criterion to miniize. 
 Functions which are available are 
\tabular{lll}{
\code{ff.bentler}   \tab orthogonal or oblique \tab  Bentler's invariant pattern simplicity criterion\cr
\code{ff.cf}		   \tab orthogonal or oblique \tab  Crawford-Ferguson family	   \cr
\code{ff.cubimax}   \tab orthogonal    \tab  	\cr
\code{ff.entropy}   \tab orthogonal \tab  minimum entropy			\cr
\code{ff.fss} 	   \tab orthogonal or oblique	   \tab  Forced Simple Structure (see Vignette)   \cr
\code{ff.geomin}    \tab orthogonal or oblique  \tab 				   \cr
\code{ff.infomax}   \tab orthogonal or oblique   \tab 				   \cr
\code{ff.oblimax}   \tab oblique	   \tab         				\cr
\code{ff.pst}	   \tab orthogonal or oblique \tab  partially specified target rotation	\cr
\code{ff.quartimax} \tab orthogonal \tab 					\cr
\code{ff.quartimin} \tab oblique	   \tab         				\cr
\code{ff.simplimax} \tab oblique	   \tab         				\cr
\code{ff.target}    \tab orthogonal or oblique \tab  target rotation			\cr
\code{ff.varimax}   \tab orthogonal \tab 					\cr
 }

 The argument \code{kappa} parameterizes the family for the Crawford-Ferguson 
 method. If \code{m} is the number of factors and \code{p} is the number of 
 items then \code{kappa} values having special names are \code{0=}Quartimax, 
 \code{1/p=}Varimax, \code{m/(2*p)=}Equamax, \code{(m-1)/(p+m-2)=}Parsimax, \code{1=}Factor parsimony.
 
 For the argument \code{kij} for Forced Simple Structure see \code{\link{rotationsDF}}.
 
 }
\examples{
  # Example from Mulaik, S.A. (2018) Fundamentals of Common Analysis in 
  # Irwing, P., Booth, T., Hughes, D. J.  (20180314). The Wiley Handbook of 
  # Psychometric Testing, 2 Volume Set: A Multidisciplinary Reference on Survey, 
  # Scale and Test Development. [VitalSource Bookshelf 10.2.0]. 
  # Retrieved from vbk://9781118489703
  # See Part II, Table 8.5 (p.225) and Table 8.6 (p.235)

  ff.absolmin <- function(L)
  {
    f <- sum(abs(L))
    list(f = f, Method = "DF-Absolmin")
  }
  absolmin.df <- function(A, Tmat = diag(ncol(A)), normalize = FALSE, eps = 1e-5, maxit = 1000){
   GPFoblq.df(A, Tmat=Tmat, method = "absolmin", normalize = normalize, eps = eps, maxit = maxit)	
  }
  # absolmin.df(fa.unrotated$loadings, normalize = TRUE)
  
  # See GPArotation vignette for a bug in factanal with oblique rotation.

  data("Griffith")  
  fa.unrotated <- factanal(covmat = Griffith, factors = 6, rotation = "none")
  GPFoblq.df(fa.unrotated$loadings, method = "quartimin", normalize = TRUE)
  
  absolmin.df(fa.unrotated$loadings, normalize = TRUE)
  
  data("Harman", package="GPArotation")
  qHarman  <- GPForth.df(Harman8, Tmat=diag(2), method="quartimax")
  
  # define a new function as ff.newname orthogonal rotation for use with factanal
  ff.expomax <- function(L)
  {
    f <- -sum(expm1(abs(L)))
    list(f = f, Method = "DF-Expomax")
  }
  GPForth.df(Harman8, method ="expomax")

  expomax.df <- function(L, Tmat = diag(ncol(L)), normalize = FALSE, eps = 1e-5, maxit = 1000){
   GPForth.df(L, Tmat=Tmat, method = "expomax", normalize = normalize, eps= eps, maxit = maxit)	
  }
  expomax.df(Harman8, normalize = TRUE)
  factanal(factors = 2, covmat = ability.cov, rotation = "expomax.df", 
      control = list(rotate =c(normalize = TRUE)))
 }
\seealso{
   \code{\link{GPForth.df}}, 
   \code{\link{GPFoblq.df}}, 
   \code{\link{fssQ.df}},
   \code{\link{fssT.df}},
   \code{\link{cubimax.df}},
   \code{\link{rotationsDF}},
   \code{\link[stats]{factanal}} 
 }
\references{
  Jennrich, R.I. (2004) Derivative free gradient projection algorithms for rotation,
  Psychometrika: \bold{69}(3), 475--480.

  }

\author{Coen A. Bernaards and Robert I. Jennrich}
\concept{rotation}
\keyword{multivariate}


