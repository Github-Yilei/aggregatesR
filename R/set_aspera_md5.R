#' @title set_aspera_md5
#' @description Message Digest Algorithm 5(MD5) is a cryptographic hash algorithm that can be used to
#' verify the integrity of files.
#' Here, we using set_aspera_md5 to prepare
#' ENA(\url{https://www.ebi.ac.uk/ena/browser/home}) md5 data.
#'
#' @param path the absolute path of the output file.
#' @param md5_tsv the ENA md5 file.
#'
#' @export
#'
#' @examples path <- 'C:/Users/Admin/Desktop/'
#' @examples md5_tsv <- system.file("extdata", "filereport_read_run_PRJNA648176_md5.txt", package = "aggregatesR" )
#' @examples set_aspera_md5(path, md5_tsv)

set_aspera_md5 <- function(path, md5_tsv){
  raw_data <- read.delim(md5_tsv, header = TRUE)
  raw_md5 <- stringr::str_split(raw_data$fastq_md5, ";")
  clean_md5 <- data.frame(unlist(raw_md5),
                          paste0(rep(raw_data$run_accession,
                                     each=2),
                                 c("_1.fastq.gz", "_2.fastq.gz")))
  # set file name
  #temp <- gsub("_tsv.txt", "", ena_tsv)
  #file_name <- gsub("filereport_read_run_", "", temp)
  write.table(clean_md5, quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              file = paste0(path,
                            "aspera_md5.txt"))
}
