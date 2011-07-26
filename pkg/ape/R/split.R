setClass("split",representation(x="matrix",n="integer",labels="vector",freq="vector"))

validSplit <-function(object)
 {
  if(is.raw(object@x[1,1]) & nrow(object@x)==ceiling(object@n/8) & length(object@labels)==object@n & length(object@freq)==ncol(object@x) & is.numeric(object@freq[1]))
   {    
    TRUE
   }
 else{str=""
  if(!is.raw(object@x[1,1])) str="x must contain elements of type raw"
  if(!length(object@labels)==object@n) str="labels must be of length n"
  if(!nrow(object@x)==ceiling(object@n/8)) str="x must have ceiling(n/8) rows"
  if(!length(object@freq)==ncol(object@x)) str="freq must have ncol(x) elements"
  if(!is.numeric(object@freq[1])) str="freq must contain numeric elements"
      paste(str)
     }
 }

setValidity("split",validSplit)

