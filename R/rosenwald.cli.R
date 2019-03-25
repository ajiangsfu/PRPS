#' rosenwald information data, which is copied from "LPS" R package
#' All information in the following are copied from https://rdrr.io/cran/LPS/man/rosenwald.html
#'
#' This dataset contains 60 Diffuse Large B-Cell Lymphomas analysed on Lymphochip microarrays, 
#' as published by Rosenwald et al. The "Germinal Center B-cell like" and "Activated B-Cell like" subtypes, 
#' as determined by hierarchical clustering, were predicted by a LPS approach in Wright et al.
#' To minimize package size, values were rounded at 3 decimals and only 60 DLBCL from the 240 series were randomly
#'  selected (40 from the "Training" set, 20 from the "Validation" set), excluding "Type III" sub-types.
#'
#' @format A data.frame with a row for each sample, and 4 factor columns described below. 
#' Rows are named by samples "LYM numbers", in the same order than rosenwald.expr
#' \describe{
#'   \item{set}{the "Training" or "Validation" set the sample comes from.}
#'   \item{group}{the DLBCL sub-type that is to be predicted ("GCB" or "ABC").}
#'   \item{follow.up}{follow-up of the patient, in years.}
#'   \item{status}{status of the patient at the end of the follow-up ("Dead" or "Alive").}
#' }
#' @source \url{https://rdrr.io/cran/LPS/man/rosenwald.html}
