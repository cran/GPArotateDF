\name{00.GPArotateDF}
\alias{GPArotateDF}
\alias{GPArotateDF-package}
\alias{00.GPArotateDF-package}
\alias{GPArotateDF.Intro}
\docType{package}
\title{Derivative-Free Gradient Projection Algorithms for Factor Rotation}
\description{Derivative-Free GPA Rotation for Factor Analysis

The GPArotateDF package contains functions for the rotation of factor loadings
matrices without the need for deriving gradients. Both orthogonal and oblique rotation 
algorithms are available. Additionally, a number of rotation criteria are provided.
The GP algorithms minimize the rotation criterion function, and provide the 
corresponding rotation matrix. For oblique rotation, the covariance / correlation matrix
of the factors is also provided. The derivartive-free rotation method implemented in this package 
is described in Jennrich (2004). The GPArotation package is based on Bernaards and Jennrich (2005). 

\tabular{ll}{
Package: \tab GPArotateDF\cr
Depends: \tab R (>= 2.0.0)\cr
License: \tab GPL Version 2.\cr
Imports: \tab GPArotation \cr
}

Index of functions:

Derivative-Free Gradient Projection Rotation Algorithms \cr
\tabular{ll}{
\code{\link{GPForth.df}} \tab Orthogonal rotation function \cr
\code{\link{GPFoblq.df}} \tab Oblique rotation function \cr
}

Rotations\cr
\tabular{ll}{
\code{\link{cubimax.df}} \tab Cubimax rotation  \cr
\code{\link{fssQ.df}} \tab Oblique Forced Simple Structure rotation \cr
\code{\link{fssT.df}} \tab Orthogonal Forced Simple Structure rotation \cr
}

ff routines to compute value of the criterion \cr
\tabular{ll}{
\code{\link{ff.bentler}}\tab Bentler's Invariant Pattern Simplicity ff \cr
\code{\link{ff.cf}} \tab Crawford-Ferguson Family ff \cr
\code{\link{ff.cubimax}} \tab Cubimax ff  \cr
\code{\link{ff.entropy}} \tab Minimum Entropy ff \cr
\code{\link{ff.geomin}}\tab Geomin ff \cr
\code{\link{ff.infomax}} \tab Infomax ff \cr
\code{\link{ff.oblimax}} \tab Oblimax ff \cr
\code{\link{ff.pst}}\tab Partially Specified Target ff \cr
\code{\link{ff.quartimax}}\tab Quartimax ff \cr
\code{\link{ff.quartimin}}\tab Quartimin ff \cr
\code{\link{ff.simplimax}}\tab Simplimax ff \cr
\code{\link{ff.fss}}\tab Forced Simple Structure ff \cr
\code{\link{ff.target}}\tab Target ff \cr
\code{\link{ff.varimax}}\tab Varimax ff \cr
}

Utility functions\cr
\tabular{ll}{
\code{\link{NormalizingWeight}} \tab Kaiser normalization (not exported from NAMESPACE) \cr
}


}



\author{Coen A. Bernaards and Robert I. Jennrich
}
\references{
   Bernaards, C.A. and Jennrich, R.I. (2005) Gradient Projection Algorithms 
     and Software for Arbitrary Rotation Criteria in Factor Analysis. 
     \emph{Educational and Psychological Measurement}, \bold{65}, 676--696.

   Jennrich, R.I. (2004). Derivative free gradient projection algorithms for rotation.
       \emph{Psychometrika}, \bold{69}, 475--480.
}
\keyword{ package }

\seealso{
   \code{\link{GPForth.df}},
   \code{\link{GPFoblq.df}}
}
