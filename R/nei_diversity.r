#' Calculates Nei diversity, given a set of genotypes (gt)
#' and a vector of weight values (equal length to the 
#' number of individuals 
#' 
#' @param gt  - genotype matrix (individuals, loci)
#' @param w   - vector of weights
#' @return Nei diversity, a float value 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

nei_diversity <- function(gt, w=NULL) {

   # estimate allele frequencies
   calc_afs <- function(gt, w=NULL) {

      if (is.null(w)) {

         p  <- colSums(gt,na.rm=TRUE) / (colSums(!is.na(gt))*2)

      } else {

         if (nrow(gt) == length(w)) {
            wm <- matrix(rep(w,ncol(gt)),ncol=ncol(gt),byrow=FALSE)
            gtnn <- gt
            gtnn[ is.na(gt) ] <- 0
            gtnn[ !is.na(gt) ] <- 1 
            p  <- colSums(gt*wm,na.rm=TRUE) / colSums( gtnn * wm *2)
         } else {
            cat("   weights supplied: must have length equal to number of rows in gt \n")
         }
      }
      return(p)
   }

   p <- calc_afs(gt, w)
   #p <- p[which(is.na(p)==FALSE)] #remove na values to calulate nei diversity 

   # calculate diversity from frequencies
   Hs <- 0
   for (i in 1:ncol(gt)) {

      Hs <- Hs + ( 1-(p[i])^2 - (1-p[i])^2  )

   }

   # finish and return
   nei <- 1 / ncol(gt) * Hs 
   return(nei)
}

