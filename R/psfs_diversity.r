#' Calculates projected sfs, given a set of genotypes (gt),
#' a value for down-projection (m), and a vector of weight 
#' values (equal length to the number of individuals) 
#' 
#' @param gt  - genotype matrix (individuals, loci)
#' @param w   - vector of weights
#' @param m   - projection value
#' @return psfs - a float value 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

psfs_diversity <- function(gt, m, w=NULL, pMAC_mode=FALSE, Nmat=NULL, ncpu=NULL, unlim_m=FALSE) {

   if (is.null(m)) {
      m = 20
   }

   # estimate allele frequencies
   if (is.null(w)) {
      wgt = gt
   } else {

      if (nrow(gt) == length(w)) {
         iw <- rep(1:nrow(gt),w)
         wgt <- gt[iw,]
      } else {
         cat("   weights supplied: must have length equal to number of rows in gt \n")
      }
   }

   if (pMAC_mode) {
      projected_SFS <- project_SFS_from_MAC(gt_MAC=wgt, gt_N=Nmat, m=m, ncpu=ncpu, unlim_m=unlim_m)
   } else {
      projected_SFS <- project_SFS_from_genotypes(gt_SNP=wgt, m=m, ncpu=ncpu, unlim_m=unlim_m)
   }

   fixedSFS <- projected_SFS[1] + projected_SFS[m+1]
   ppSFS <- (sum(projected_SFS) - fixedSFS)/sum(projected_SFS)

     
  if (is.na(ppSFS)){
      ppSFS <- 0
  }

   return(ppSFS)
}

