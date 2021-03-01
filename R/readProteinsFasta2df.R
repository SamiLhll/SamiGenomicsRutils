#' @title reads a proteins fasta file as a data_frame
#' @description reads a proteins fasta file as a data_frame for further protein wise analyses
#'
#' @param fasta_file a path to the FASTA file
#'
#' @return a datat_frame composed of two columns protName with the preotein names and protSeq containing the fasta sequence
#'
#' @importFrom protr readFASTA
#' @export
#'


readProteinsFasta2df <- function(fasta_file){
  ### get protnames as list :
  temp_command <- paste0("grep \">\" ",fasta_file," | awk '{print $1}' | sed -e 's/>//g'")
  prot_names <- as.character(system(command = temp_command,
                                    intern = TRUE))
  ### read fasta file with protr package :
  protr_fasta_proteins <- readFASTA(fasta_file)
  ### convert to df :
  protein_sequences <- NULL
  for (i in prot_names) {
    protein_sequences <- c(protein_sequences,
                           toString(protr_fasta_proteins[i][[1]]))
  }
  protein_fasta_df <- data.frame(protName = prot_names,
                                 protSeq = protein_sequences)
  ### return dataframe :
  return(protein_fasta_df)
}
