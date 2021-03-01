#' @title calculates proteins aa compositions and returns the results as a data frame
#' @description calculates proteins aa compositions and returns the results as a data frame.
#'
#' @param proteins a data frame loaded with the readProteinsFasta2df function
#' @param proteins_subset an array of protein names to subset from the dataset provided as the proteins parameter (default = NULL)
#'
#' @return a datat_frame composed of the protName as first column and the composition of each AA in the next columns
#'
#' @importFrom protr extractAAC
#' @export
#'

calculateIndividualAAComposition2df <- function(proteins,proteins_subset = NULL){
  #### Prepare the set of proteins to process based on if it needs to be subset or not :
  if (is.null(proteins_subset)) { proteins_to_process <- proteins }
  else {proteins_to_process <- proteins[proteins$protName %in% proteins_subset]}
  temp_df <- NULL
  #### iterate on every row (== protein sequence) individually
  for (i in proteins_to_process$protSeq) {
    #### remove unrecognized characters ("*","X","U") :
    temp_stars_removed <- gsub("*","",as.character(i),fixed = TRUE)
    temp_stars_X_removed <- gsub("X+","", temp_stars_removed)
    temp_stars_X_U_removed <- gsub("U+","",temp_stars_X_removed)
    #### calculate the AA compositions using extractAAC() function from the protr package
    temp <- extractAAC(temp_stars_X_U_removed)
    #### convert to data.frame and append all the data
    temp_df <- rbind(temp_df,data.frame(t(temp)))
  }
  ### post_process the dataframe to get rid of the column with the fasta sequence :
  temp_proteins_df <- cbind(proteins,temp_df)
  ### drop the protein sequence from the table
  temp_proteins_df <- temp_proteins_df[, names(temp_proteins_df) != "protSeq"]
  # round all the scores in the dataframe :
  temp_proteins_df[,-1] <- round(temp_proteins_df[,-1], digits = 3)
  return(temp_proteins_df)
}
