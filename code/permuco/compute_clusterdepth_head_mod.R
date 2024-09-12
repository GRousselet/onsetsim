# modified from compute_clusterdepth_head, GAR, University of Glasgow, January 2024
# set alternative to greater as we use t2 values.
# cluster matrix calculated separately as we use the permutation distribution to form clusters; 
# otherwise assumes a unique threshold for all time points, which implies constant variance -- here var.equal = FALSE, as it should always be. 

############################################################################################################
#' Cluster-depth correction (from the head only)
#'
#' @description Compute the cluster-depth test correction (from the head) given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param threshold A scalar that represents the threshold to create the clusters.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @export
#' @family multcomp
compute_clusterdepth_head_mod <- function(distribution, cluster){

  alternative <- "greater"
  
  depth_head <- get_clusterdepth_head(cluster, border = "ignore")
  distr_head <- depth_distribution_head(distribution, depth_head)

  pvalue <- rep(NA,ncol(cluster))
  max_cl_size <- max(table(cluster[1,cluster[1,]!=0]))

  for(cli  in seq_len(max(cluster[1,]))){
    sample <- which(cluster[1,]==cli)
    di <- distr_head#[,seq_len(max_cl_size), drop=FALSE]
    di <- rbind(c(distribution[1, sample],rep(0,ncol(di)-length(sample))),di)
    pvalue[sample] <- compute_troendle(di,
                                       alternative = alternative)$main[seq_along(sample),2]
  }

  list(main = cbind(statistic = distribution[1,], pvalue = pvalue,cluster_id = cluster[1,]))

}