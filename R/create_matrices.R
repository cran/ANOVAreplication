#' Create (in)equality constraint matrices
#'
#' Parses a character string describing an informative hypothesis,
#' and returns (in)equality constraint matrices
#'
#' Informative hypotheses specified as a character string by "hyp" should
#' adhere to the following simple syntax: \itemize{
#' \item The hypothesis consists of a (series of) (in)equality
#' constraint(s). Every single (in)equality constraint is of the form "R1*mu1 +
#' R2*mu2+... = r", where capital Rs refer to numeric scaling constants, must
#' refer to the names of parameters in the model, and the lower case r refers
#' to a constant. Standard mathematical simplification rules apply; thus,
#' "R1*mu1 = R2*mu2" is equivalent to "R1*mu1 - R2*mu2 = 0".  \item Multiple
#' unrelated constraints within one hypothesis can be chained by "&". Thus,
#' "a=b&c=d" means that H1: a=b AND c=d.  \item Multiple related constraints
#' within one hypothesis can be chained by repeating the (in)equality operators
#' "=", "<", or ">". Thus, "a<b<c" means that H1: a < b AND b < c.  \item
#' Parameters can be grouped by placing them in a parenthesized, comma
#' separated list. Thus, "(a,b)>c" means that H1: a > c AND b > c.  Similarly,
#' "(a,b)>(c,d)" means that H1: a > c AND b > c AND b > c AND b > d.  }
#'
#' @aliases create_matrices create_matrices
#' @param varnames A character (vector of characters), containing names of
#' variables used in the hypotheses.  %Object of class \code{\link{lm}}, from
#' which the model parameters are extracted.
#' @param hyp A character string, containing a Bain hypothesis (see Details).
#' @return A pair of named matrices for every hypothesis specified in the
#' \code{hyp} argument; one matrix named ERr, specifying equality constraints,
#' and one matrix named IRr, specifying inequality constraints.
#' @author Caspar van Lissa
#' @keywords internal utilities
create_matrices <- function(varnames, hyp){
  if(is.null(varnames)) stop("Please input proper linear model object")
  hyp <- gsub("\\s", "", hyp)
  if(grepl("[><=]{2,}", hyp)) stop("Do not use combined comparison signs e.g., '>=' or '==', and use '&' to add a related constraint instead of ','")

  hyp_list <- strsplit(hyp, ";")[[1]]                         #mz
  hyp_list <- lapply(hyp_list, function(x){ strsplit(x, "&")[[1]]})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_compound_constraints))})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_parentheses))})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, flip_inequality)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, constraint_to_equation)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, order_terms)})

  hyp_list <- unlist(lapply(hyp_list, function(x){
    ERr <- x[grep("=", x)]
    IRr <- x[grep("[<>]", x)]
    if(length(ERr) == 0){
      ERr <- NULL
      E <- 0                                                    #mz
    } else {
      ERr <- t(sapply(ERr, constraint_to_row, varnames = varnames))
      E <- nrow(ERr)                                            #mz
      #colnames(ERr)[ncol(ERr)]<- "="                           #mz
    }
    if(length(IRr) == 0){
      IRr <- NULL
    } else {
      IRr <- t(sapply(IRr, constraint_to_row, varnames = varnames))
      #colnames(IRr)[ncol(IRr)]<- ">"                           #mz
    }

    if(E==0){
      R <- IRr[,-ncol(IRr)]}else{                               #mz
        if(length(IRr) == 0){R <- ERr[,-ncol(ERr)]}else{        #mz
          R <- rbind(ERr[,-ncol(ERr)],IRr[,-ncol(IRr)])}}       #mz
    if(is.vector(R)==TRUE){R<-matrix(R,nrow=1)}
    #mz
    r <- c(ERr[,ncol(ERr)],IRr[,ncol(IRr)])                     #mz
    list(Amat = R, difmin = r, E=E)                             #mz

  }), recursive = FALSE)

  #names(hyp_list) <- paste0(names(hyp_list), rep(1:(length(hyp_list)/2), each = 2))     #mz
  hyp_list
}







