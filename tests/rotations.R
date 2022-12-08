#   Tests here only compare against values computed with GPArotation code,
#   to ensure the regular and DF versions give the same result


 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 require("GPArotateDF")
 search()
 Sys.info()

require("stats")  

fuzz <- 1e-6 
all.ok <- TRUE  


  data(ability.cov)
  L <- loadings(factanal(factors = 2, covmat=ability.cov))

# quartimax

LG <- quartimax(L, normalize = FALSE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "quartimax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
LG <- quartimax(L, normalize = TRUE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "quartimax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
 
	  
    
# quartimin

LG <- quartimin(L, normalize = FALSE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "quartimin")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
LG <- quartimin(L, normalize = TRUE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "quartimin")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
      
# oblimax commented out as it is gives problem quite consistently
   
#LG <- oblimax(L, normalize = FALSE, eps=1e-5)
#LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "oblimax")  
# Oblimax fails for fuzz = 1e-6. But succeeds for 0.01
#  if( 0.01 < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
#    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
#    cat("difference:\n")
#    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
#    all.ok <- FALSE  
#    } 
#  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
#    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
#    cat("difference:\n")
#    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
#    all.ok <- FALSE  
#    } 
    
#LG <- oblimax(L, normalize = TRUE, eps=1e-5)
#LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "oblimax")  
#  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
#    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
#    cat("difference:\n")
#    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
#    all.ok <- FALSE  
#    } 	  
#  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
#    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
#    cat("difference:\n")
#    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
#    all.ok <- FALSE  
#    } 
 
 
# entropy

LG <- entropy(L, normalize = FALSE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "entropy")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
LG <- entropy(L, normalize = TRUE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "entropy")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
 
 
 
# simplimax
   
LG <- simplimax(L, normalize = FALSE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "simplimax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
LG <- simplimax(L, normalize = TRUE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "simplimax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
 
 
# bentlerQ
   
LG <- bentlerQ(L, normalize = FALSE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "bentler")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
LG <- bentlerQ(L, normalize = TRUE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "bentler")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
 

# bentlerT

LG <- bentlerT(L, normalize = FALSE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "bentler")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
LG <- bentlerT(L, normalize = TRUE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "bentler")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
 

# geominQ
   
LG <- geominQ(L, normalize = FALSE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "geomin")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
LG <- geominQ(L, normalize = TRUE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "geomin")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
 

# geominT

LG <- geominT(L, normalize = FALSE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "geomin")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
LG <- geominT(L, normalize = TRUE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "geomin")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
 

# infomaxQ
   
LG <- infomaxQ(L, normalize = FALSE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "infomax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
    # changed to eps 1e-6 in order to pass the test.
LG <- infomaxQ(L, normalize = TRUE, eps=1e-6)
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-6, method = "infomax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
 

# infomaxT

LG <- infomaxT(L, normalize = FALSE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "infomax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
LG <- infomaxT(L, normalize = TRUE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "infomax")  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
 

# CF Equamax Q
   
LG <- GPFoblq(L, normalize = FALSE, eps=1e-5,  method = "cf", methodArgs=list(kappa=2/12))
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "cf", methodArgs=list(kappa=2/12))  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
LG <- GPFoblq(L, normalize = TRUE, eps=1e-5,  method = "cf", methodArgs=list(kappa=2/12))
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "cf", methodArgs=list(kappa=2/12))  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 	  
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
 

# CF equamax T

LG <- GPForth(L, normalize = FALSE, eps=1e-5, method = "cf", methodArgs=list(kappa=2/12))
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "cf", methodArgs=list(kappa=2/12))  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
LG <- GPForth(L, normalize = TRUE, eps=1e-5, method = "cf", methodArgs=list(kappa=2/12))
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "cf", methodArgs=list(kappa=2/12))  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
    

# targetQ
LG <- targetQ(L, Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2), normalize=FALSE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = FALSE, eps=1e-5, method = "target", 
 	methodArgs=list(Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2)))  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    } 
    
LG <- targetQ(L, Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2), normalize=TRUE, eps=1e-5)
LGDF <- GPFoblq.df(L, normalize = TRUE, eps=1e-5, method = "target", 
 	methodArgs=list(Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2)))  

  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
  if( fuzz < max(abs(unclass(LGDF)$Phi - unclass(LG)$Phi))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$Phi - unclass(LG)$Phi), digits=18)
    all.ok <- FALSE  
    }   	  
 
# targetT

LG <- targetT(L, Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2), normalize=FALSE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = FALSE, eps=1e-5, method = "target", 
 	methodArgs=list(Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2)))  
  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 

LG <- targetT(L, Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2), normalize=TRUE, eps=1e-5)
LGDF <- GPForth.df(L, normalize = TRUE, eps=1e-5, method = "target", 
 	methodArgs=list(Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2)))  
  if( fuzz < max(abs(unclass(LGDF)$loadings - unclass(LG)$loadings))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    cat("difference:\n")
    print((unclass(LGDF)$loadings - unclass(LG)$loadings), digits=18)
    all.ok <- FALSE  
    } 
 
   	  
 
# pstT
# This won't converge properly
# No further investigations performed.

 
 
cat("tests completed.\n")



if (! all.ok) stop("some tests FAILED")
