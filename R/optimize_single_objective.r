#' Performs simulated annealing to find mixture of individuals 
#' that maximize a function (e.g. a measure of genetic diversity) 
#' 
#' @param gt                - a genotype matrix (samples, loci)
#' @param initial_weights   - a vector of initial weights
#' @param weights_max       - maximum weights for each individual
#' @param weights_min       - minimum weights for each individual
#' @param max_t             - maximum temperature 
#' @return a list of output 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


optimize_single_objective <- function( gt=NULL, sm=NULL, measure=NULL, max_steps=10000, N_t=NULL, initial_weights=NULL, weights_max=NULL, weights_min=NULL, max_t=1, q=NULL, m=NULL, p_depends_delta=TRUE, disp=0, pMAC_mode=FALSE, Nmat=NULL, ncpu=NULL, unlim_m = NULL) {

  proceed = TRUE
  if (is.null(gt) & is.null(sm)) {
    cat("   Provide either initial_weights or a number of target individuals (N_t) \n")
    proceed = FALSE
  }
  if (is.null(unlim_m)) {
    unlim_m <- FALSE
  }
  if (!is.null(gt) & !is.null(sm)) {
    cat("   Provide either a genotype or dimilarity matrix \n")
    proceed = FALSE
  }
  N_g = NULL
  if (!is.null(gt)) {
    N_g = nrow(gt)
    v = gt
  }
  if (!is.null(sm)) {
    N_g = nrow(sm)
    v = sm
  }
  if (is.null(N_t) & is.null(initial_weights)) {
    cat("   Provide either initial_weights or a number of target individuals (N_t) \n")
    proceed = FALSE
  }
  if (!is.null(N_t) & !is.null(initial_weights)) {
    if (N_t != sum(initial_weights)) {
      cat("   Conflict in number of target individuals and initial weights provided \n")
      proceed = FALSE
    }
  }
  if (!is.null(N_t) & is.null(initial_weights)) {
    cat("   Generating initial weights  \n")
    initial_weights <- propose_initial_weights(N_g, N_t)
  }
  if (measure == "psfs") {
    require(sfsCalcs)
  }
  summary <- NULL
  summary <- generate_measure_value(v = v, measure = measure, 
                                    w = initial_weights, q = q, m = m, disp = disp, pMAC_mode = pMAC_mode, 
                                    Nmat = Nmat, ncpu = ncpu, unlim_m = unlim_m)
  if (is.null(summary)) {
    proceed = FALSE
    cat("   Initial weights returned null value for measure \n")
  }
  if (!proceed) {
    cat("   Something has gone wrong. Reconsider data or settings. \n")
  }
  else {
    weight <- mat.or.vec(max_steps, N_g)
    value <- mat.or.vec(max_steps, 1)
    accept <- mat.or.vec(max_steps, 1)
    chain <- list(weight = weight, value = value, accept = accept)
    chain$value[1] <- summary
    chain$weight[1, ] <- initial_weights
    s <- 2
    weights <- initial_weights
    while (s <= max_steps) {
      proposed_weights <- propose_new_weights(weights, 
                                              w_max = weights_max, w_min = weights_min)
      proposal_summary <- generate_measure_value(v, measure = measure, 
                                                 w = proposed_weights, q = q, m = m, disp = disp, 
                                                 pMAC_mode = pMAC_mode, Nmat = Nmat, ncpu = ncpu, 
                                                 unlim_m = unlim_m)
      if (is.na(proposal_summary)) {
        proposal_summary <- 0
      }
      temp <- temp_scheduler(s, max_steps, max_t = max_t)
      
      #if kinship isn't nulll, then we run an additioanl kinship check to ensure no klinship between proposed samples to be optimised
      if(is.null(kinall)){
        accept_proposal <- proposal_accept_reject(summary = summary,proposal_summary = proposal_summary, t=temp, p_depends_delta = p_depends_delta)
        } else {
          
          #check for any kinship here and reject if there is:
          sampstotestforkin <- rownames(gt_sw_comp2)[which(proposed_weights>0)]
          kincombos <- combn(sampstotestforkin,2)
          finalkinval <- 0
          for (u in (length(kincombos)/2)){
            finalkinval <- as.numeric(finalkinval + kinall$value[which(kinall$Var1==kincombos[,u][1]&kinall$Var2==kincombos[,u][2])])
          }
          
          if (finalkinval>0){
            accept_proposal <- FALSE
          } else{
            accept_proposal <- proposal_accept_reject(summary = summary,proposal_summary = proposal_summary, t=temp, p_depends_delta = p_depends_delta)
          }
      }
      
      
      if (accept_proposal) {
        weights <- proposed_weights
        summary <- proposal_summary
      }
      chain$accept[s] <- accept_proposal
      chain$value[s] <- summary
      chain$weight[s, ] <- weights
      s <- s + 1
      pfreq <- floor(max_steps/20)
      if (s%%pfreq == 0) {
        cat("   up to step", s, "\n")
      }
    }
    return(chain)
  }
}
