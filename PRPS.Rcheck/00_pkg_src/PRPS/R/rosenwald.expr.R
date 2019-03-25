#' rosenwald expression data, which is copied from "LPS" R package
#' All information in the following are copied from https://rdrr.io/cran/LPS/man/rosenwald.html
#'
#' This dataset contains 60 Diffuse Large B-Cell Lymphomas analysed on Lymphochip microarrays, 
#' as published by Rosenwald et al. The "Germinal Center B-cell like" and "Activated B-Cell like" subtypes, 
#' as determined by hierarchical clustering, were predicted by a LPS approach in Wright et al.
#' To minimize package size, values were rounded at 3 decimals and only 60 DLBCL from the 240 series were randomly
#'  selected (40 from the "Training" set, 20 from the "Validation" set), excluding "Type III" sub-types.
#'
#' @format A numeric matrix of expression values, with probes in rows and samples in columns. 
#' Both dimensions are named, probes by there "UNIQID" and samples by there "LYM numbers". Many NA values are present.
#' \describe{
#'   \item{rosenwald.expr}{60 Diffuse Large B-Cell Lymphomas analysed on Lymphochip microarrays}
#' }
#' @source \url{https://rdrr.io/cran/LPS/man/rosenwald.html}
