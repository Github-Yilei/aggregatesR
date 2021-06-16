#' @title set_aspera_fq
#' @description setting the command lines of downloading fastq data from
#' ENA(\url{https://www.ebi.ac.uk/ena/browser/home}) by using
#' IBM Aspera.
#' @param path the absolute path of the output file.
#' @param ena_tsv the ENA file.
#' @param ascp the absolute path of the ascp package.
#' @param step the jobs of each file.
#'
#' @param openssh the openssh of ascp package(at etc file).
#' @param fastq_cmd the command line for downloading fastq data(please checking it by hand).
#'
#' @export
#'
#' @examples path <- 'C:/Users/Admin/Desktop/'
#' @examples ena_tsv <- system.file("extdata", "filereport_read_run_PRJNA648176_tsv.txt", package = "aggregatesR" )
#' @examples ascp <- '~/miniconda3/pkgs/aspera-cli-3.9.1-0'
#' @examples set_aspera_fq(path, ena_tsv, ascp, 10)

set_aspera_fq <- function(path, ena_tsv, ascp, step){
  raw_data <- read.delim(ena_tsv, header = TRUE)
  openssh <- paste0(ascp, '/etc/asperaweb_id_dsa.openssh')
  fastq_cmd <- paste(ascp,
                     '/bin/ascp -i',
                     openssh,
                     '--overwrite=diff -P33001 -T -l 30m era-fasp@')
  # fastq_aspera/fastq_ftp
  raw_links <- stringr::str_split(raw_data$fastq_aspera, ";")
  links <- unlist(raw_links)
  links <- paste0(fastq_cmd, links, " ./")

  dd2 <- data.frame()
  index <- seq(1, length(links)/step, 1)

  for (i in index){
    sub_links <- data.frame(links[(1 + step*(i-1)):(i*step)])
    names(sub_links) <- "links"
    df <- data.frame(links= c(
      paste0("echo start download reads ", 1 + step*(i-1),
             "-", i*step),
      paste0("mkdir " , (1 + step*(i-1)), "-", (i*step)),
      paste0("cd ", (1 + step*(i-1)), "-", (i*step)))
    )
    dd <- rbind(df, sub_links)
    dd2 <- rbind(dd2, dd)
  }


  left_links <- data.frame(links[(max(index)* step + 1): length(links)])
  names(left_links) <- "links"

  df <- data.frame(links= c(
    paste0("echo start download reads ",
           max(index)* step + 1, "-", length(links)),
    paste0("mkdir " , max(index)* step + 1,
           "-", length(links)),
    paste0("cd ", max(index)* step + 1,
           "-", length(links))
  ))

  my_df <- rbind(dd2, df, left_links)

  # set file name
  temp <- gsub("_tsv.txt", "", ena_tsv)
  file_name <- gsub(".*filereport_read_run_", "", temp)
  write.table(my_df, quote = FALSE, row.names = FALSE,
              col.names = FALSE,
              file = paste0(path, file_name, "_links.txt"))
}
