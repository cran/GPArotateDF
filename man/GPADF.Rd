\name{GPADF}
\alias{GPForth.df}
\alias{GPFoblq.df}

\title{Rotation Optimization}
\usage{
    GPForth.df(A, Tmat=diag(ncol(A)), normalize = FALSE, eps=1e-5, 
       maxit=1000, method="varimax", methodArgs=NULL)
    GPFoblq.df(A, Tmat=diag(ncol(A)), normalize = FALSE, eps=1e-5, 
       maxit=1000, method="quartimin", methodArgs=NULL)
    }
\arguments{
    \item{A}{initial factor loadings matrix for which the rotation 
             criterian is to be optimized.}
    \item{Tmat}{initial rotation matrix.}
    \item{normalize}{see details.}
    \item{eps}{convergence is assumed when the norm of the gradient is smaller
         than eps.}
    \item{maxit}{maximum number of iterations allowed in the main loop.}
    \item{method}{rotation objective criterian.}
    \item{methodArgs}{a list ofmethodArgs arguments passed to the rotation objective}
}

\description{
   Derivative free gradient projection rotation optimization routine used by various rotation
   objective.
   }
\value{A GPArotation object which is a list with elements 
  \item{loadings}{The rotated loadings, one column for each factor. If randomStarts
  		were requested then this is the rotated loadings matrix with the 
  		lowest criterion value.}
  \item{Th}{The rotation matrix, loadings \%*\% t(Th) = A.}
  \item{Table}{A matrix recording the iterations of the rotation optimization.}
  \item{method}{A string indicating the rotation objective function.}
  \item{orthogonal}{A logical indicating if the rotation is orthogonal.}
  \item{convergence}{A logical indicating if convergence was obtained.}
  \item{Phi}{t(Th) \%*\% Th.  The covariance matrix of the rotated factors.
  This will be the identity matrix for orthogonal
  rotations so is omitted (NULL) for the result from GPForth.df.}
  \item{G}{The gradient of the objective function at the rotated loadings.}
  }
\details{
 Derivative free gradient projection rotation optimization routines can be used to
 rotate a loadings matrix. The rotation criteria in the \code{GPArotation} package
 require a derivative to operate. In certain cases, the derivative is complex
 or non-existent. The derivative free gradient projection method provides a numerical
 alternative to the \code{GPArotation} package. 
 The functions in the package \code{GPArotateDF} follow most of the functionality
 and logic as in the \code{GPArotation} package. Please consult the documentation 
 in \code{GPArotation} for further details. 
 
 The argument \code{method} can be used to specify a string indicating
 the rotation objective. \code{GPFoblq} defaults to \code{"quartimin"} 
 and \code{GPForth} defaults to \code{"varimax"}. Available rotation objective functions
 include \code{"ff.bentler"}, \code{"ff.cf"}, \code{"ff.cubimax"}, \code{"ff.entropy"}, 
 \code{"ff.fss"}, \code{"ff.geomin"}, \code{"ff.infomax"}, \code{"ff.oblimax"},
 \code{"ff.pst"}, \code{"ff.quartimax"},\code{"ff.quartimin"}, \code{"ff.simplimax"}, 
 \code{"ff.target"}, and \code{"ff.varimax"}. 
 Most of the rotation criteria are avaible in the GPArotation pacakage 
 except for \code{cubimax} and \code{Forced Simple Structure}.

The rotation criteria are in the functions prefixed by "ff." that are used 
in the actual function call. The \code{ff.*} function call
 would typically not be used directly, but are needed for rotation. Since
 these are illustrative of computation, these are all exported
 from the package namespace. 
 New criteria for use with derivative free GP rotation do require a function of the type
 \code{ff.newCriterionName} that provides value for complexity \code{f}, and name of method.

 Some rotation criteria (including \code{"simplimax"}, \code{"pst"}, 
 \code{"target"}, \code{"cf"}, \code{"fss"}) require one or more additional arguments. 
 Check GPArotation documentation for details or see \code{\link{ff.fss}}. 
 
 The argument normalize gives an indication of if and how any normalization should 
 be done before rotation, and then undone after rotation. If normalize is FALSE 
 (the default) no normalization is done. If normalize is TRUE then Kaiser normalization 
 is done. (So squared row entries of normalized A sum to 1.0. This is sometimes 
 called Horst normalization.) If normalize is a vector of length equal to the 
 number of indicators (= number of rows of A) then the colums are divided by 
 normalize before rotation and multiplied by normalize after rotation. 
 If normalize is a function then it should take A as an argument and return a vector 
 which is used like the vector above.

}
\seealso{
	\code{\link{cubimax.df}}
	\code{\link{fssQ.df}}
	\code{\link{fssT.df}}
	\code{\link{ff.bentler}}, 
	\code{\link{ff.cf}}, 
	\code{\link{ff.cubimax}}, 
	\code{\link{ff.entropy}},
	\code{\link{ff.fss}}, 
	\code{\link{ff.geomin}}, 
	\code{\link{ff.infomax}}, 
	\code{\link{ff.oblimax}}, 
	\code{\link{ff.pst}},
	\code{\link{ff.quartimax}},
	\code{\link{ff.quartimin}}, 
	\code{\link{ff.simplimax}}, 
	\code{\link{ff.target}}, 
	\code{\link{ff.varimax}}
}
\examples{
  # GPRSorth and rotation name 
  data("Harman", package = "GPArotation")
  GPForth.df(Harman8, method = "quartimax")
  GPForth.df(Harman8, method = "cubimax")
  GPForth.df(Harman8, method = "varimax")
  GPFoblq.df(Harman8, method = "quartimin")

  # displaying results of factor analysis rotation output
  origdigits <- options("digits")
  Abor.unrotated <- factanal(factors = 2, covmat = ability.cov, rotation = "none")
  Abor <- GPFoblq.df(loadings(Abor.unrotated), method = "quartimin")
  Abor
  print(Abor)
  print(Abor, Table = TRUE)
  print(Abor, digits = 2)
  summary(Abor)
  options(digits = origdigits$digits)
  }
\author{Coen A. Bernaards and Robert I. Jennrich
		with some R modifications by Paul Gilbert.}
\references{
   Jennrich, R.I. (2004). Derivative free gradient projection algorithms for rotation. 
    \emph{Psychometrika}, \bold{69}, 475--480.
       
   Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.
   }
\concept{rotation}
\keyword{multivariate}
