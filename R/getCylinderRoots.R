#' Characteristic root of an inhomogeneous truncated shear beam
#'
#' Returns the eigenvalue \code{an} of mode \code{no} for a Gazetas &
#' Dakoulas (1985) inhomogeneous truncated shear beam parameterised by
#' inhomogeneity ratio \code{mo} and truncation ratio \code{lo}. The
#' value is interpolated from the precomputed `CylinderRoots` table
#' using the requested model (linear glm, non-linear glm, decision
#' tree, or random forest).
#'
#' @param mo Double. Inhomogeneity ratio of the shear-modulus profile.
#'   \code{mo = 0} is a homogeneous column. Range 0 to 0.95.
#' @param lo Double. Truncation ratio (top-radius / bottom-radius);
#'   \code{lo = 0} is a fully truncated cone. Range 0 to 0.495.
#' @param no Integer. Mode number, 1..8. Default 1 (fundamental).
#' @param model Character. Interpolation model: \code{"lm"} (linear
#'   GLM), \code{"nlm"} (non-linear GLM with interactions, default),
#'   \code{"dt"} (decision tree), or \code{"rf"} (random forest).
#' @param extrapolate Logical. If \code{TRUE} (default), allow
#'   \code{mo}/\code{lo} outside the table support and warn; if
#'   \code{FALSE}, clamp to the table maxima.
#' @param OSF Double. Outlier scale factor used by the random-forest
#'   model to shrink the training set around the query point. Default
#'   0.10.
#'
#' @return Double. Characteristic root \code{an} for the requested
#'   mode and parameters.
#' @export getCylinderRoots
#'
#' @import data.table
#' @import stats
#' @import rpart
#' @import randomForest
#'
#' @examples
#' getCylinderRoots(mo = 0.45, lo = 0.20, no = 1, model = "nlm")

getCylinderRoots <- function(mo,lo,no=1,model="nlm",extrapolate=TRUE,OSF=0.10){

  # Ensure the lazy-loaded dataset is bound in this scope. Without
  # this call, bare-name references below resolve via the package
  # namespace under `library(newmark)` but can fail to resolve when
  # the package is loaded via `pkgload::load_all()` (the path pkgdown
  # takes when install = FALSE), producing
  #   "object 'CylinderRoots' not found".
  utils::data("CylinderRoots", package = "newmark", envir = environment())

  OK <- exists("CylinderRoots") & !is.null(CylinderRoots)
  stopifnot(OK)
  . <- NULL
  # Check ranges
  if(!(no %in% seq(min(CylinderRoots$n), max(CylinderRoots$n)) )){
    stop(sprintf("no=%d .Wrong mode number",no))}
  DATA <- CylinderRoots[n==no,-c("n")] |> unique()
  # Check if you were lucky
  if(nrow(DATA[l==lo & m==mo])>0){
    VALUE <- DATA[l==lo & m==mo]$an
    return(VALUE)
  }
  # Check ranges
  LoMax <- max(DATA$l)
  if(lo>LoMax & extrapolate==FALSE){
    lo <- min(lo,LoMax)
  }
  
  MoMax <- max(DATA$m)
  if(mo>MoMax & extrapolate==FALSE){
    mo <- min(mo,MoMax)
  }
  
  # Linear interpolation
  if(model=="lm"){
    MODEL <- glm(an ~l+m+l2+m2,data=DATA[,.(an,l,m,l2=l^2,m2=m^2)])
    NEWDATA <- list(l=lo,m=mo,l2=lo^2,m2=mo^2)
    VALUE <- predict(MODEL,newdata=NEWDATA)
    
  }
  # Non-Linear interpolation. Mixed effects. Best 
  if(model=="nlm"){
    MODEL <- glm(an ~l+m+l2+m2+l*m+l2*m2,data=DATA[,.(an,l,m,l2=l^2,m2=m^2)])
    NEWDATA <- list(l=lo,m=mo,l2=lo^2,m2=mo^2)
    VALUE <- predict(MODEL,newdata=NEWDATA)
    
  }
  
  # Decision Trees
  if(model=="dt"){
    MODEL <- rpart(an ~ l+m,data=DATA)
    NEWDATA <- list(l=lo,m=mo)
    VALUE <- predict(MODEL,newdata=NEWDATA)
    }
  
  
  # Random Forest
  if(model=="rf"){
    # Reduce the dataset
    DATA <- DATA[ between(m,(1-OSF)*mo,(1+OSF)*mo) & between(l,(1-OSF)*lo,(1+OSF)*lo)]
    MODEL <- randomForest(an ~ l+m,data=DATA,importance=FALSE,proximity=FALSE)
    NEWDATA <- list(l=lo,m=mo)
    VALUE <- predict(MODEL,newdata=NEWDATA)
  }
  if(lo>LoMax & mo<=MoMax & extrapolate==TRUE){
    warning(sprintf("lo=%g is higher than LoMax=%g. The rooth an=%g is a predicted value that may not satisfy the boundary conditions of characteristics (cylinder) equation. Use this value with proper judgement.",lo,LoMax,VALUE))
  }
  if(mo>MoMax & lo<=LoMax & extrapolate==TRUE){
    warning(sprintf("mo=%g is higher than MoMax=%g. The root an=%g is a predicted value that may not satisfy the boundary conditions of characteristics (cylinder) equation. Use this root value with proper judgement.",mo,MoMax,VALUE))
  }
  if(mo>MoMax & lo>LoMax & extrapolate==TRUE){
    warning(sprintf("mo=%g and lo=%g  are higher than MoMax=%g and LoMax=%g, respectively. The root an=%g is a predicted value that may not satisfy the boundary conditions of characteristics (cylinder) equation. Use this root value with proper judgement.",mo,lo,MoMax,LoMax,VALUE))
  }
  return(VALUE)
}
