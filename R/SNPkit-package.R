#' SNPkit: S4 tools for reading and organizing genetic data
#'
#' Utilities for reading, cleaning, summarizing, and modeling SNP genotype data.
#'
#' @keywords internal
#' @docType package
#' @name SNPkit
#' @author
#' Vinícius Junqueira \email{junqueiravinicius@hotmail.com}
#' Roberto Higa \email{roberto.higa@embrapa.br}
#' Fernando Flores Cardoso \email{fernando.cardoso@embrapa.br}
#' Marcos Jun Iti Yokoo \email{marcos.yokoo@embrapa.br}
"_PACKAGE"

dummy_imports <- function() {
  MASS::isoMDS
  anticlust::fast_anticlustering
  dplyr::select
  ggplot2::ggplot
  ggplot2::aes
  ggplot2::geom_point
  ggplot2::labs
  ggplot2::theme_minimal
  ggplot2::theme
  ggplot2::element_rect
  ggplot2::ggsave
  grDevices::dev.off
  graphics::hist
  graphics::par
  graphics::text
  magrittr::`%>%`
  reshape2::acast
  snpStats::col.summary
  snpStats::row.summary
  snpStats::snp.pre.multiply
  snpStats::snp.post.multiply
  stats::dist
  stats::hclust
  stats::pchisq
  stats::prcomp
  stats::sd
  utils::read.table
  utils::write.table
  Rcpp::evalCpp
  invisible(NULL)
}

if (FALSE) {
  dummy_imports()
}
