#' @name getMotifFromMEME
#'
#' @title Extract Motif Information from MEME Software.
#' @description
#' \code{getMotifFromMEME} Extract motif information from the MEME software results.
#'
#' @param data A txt or xml file from MEME software.
#' @param format txt or xml.
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_pad
#' @importFrom dplyr select mutate
#' @importFrom utils read.table
#' @importFrom tidyr pivot_longer unite
#' @importFrom XML xmlParse xmlRoot xmlSize xmlToList


#'
#' @examples
#'
#' filepath <- system.file("examples", "meme.txt", package = "aggregatesR")
#' file <- getMotifFromMEME(data = filepath, format = "txt")
#' motif_seq <- getMotifFromMEME(file, format = "txt")
#'
#' filepath <- system.file("examples", "meme.xml", package = "aggregatesR")
#' file <- getMotifFromMEME(data = filepath, format = "xml")
#' gene_info <- getMotifFromMEME(file, format = "xml")
#'
#'
#' @export
#'


getMotifFromMEME <- function(data, format = "txt") {
  if (format == "txt") {
    df <- read.table(data, sep = ",") %>% dplyr::rename(raw = V1)

    motif_seq_start <- which(grepl("Motif [0-9]* in BLOCKS format", df$raw)) + 3
    motif_seq_end <- which(grepl("^//$", df$raw)) - 1

    motif_seq <- NULL
    for (i in c(1:length(motif_seq_start))){
      tmp <- df[motif_seq_start[i]:motif_seq_end[i], ] %>% as.data.frame()
      colnames(tmp) <- "raw"

      tmp %>%
        separate(col = raw, sep = " *\\(.*\\) ", into = c("seq_name", "remaind"), remove = FALSE) %>%
        separate(col = remaind, into = c("motif_seq"), sep = "  ", extra = "drop") %>%
        mutate(motif_num = paste0("Motif_", i)) -> tmp

      motif_seq <- rbind(motif_seq, tmp)
    }

    return(motif_seq )

  } else {
    df <- XML::xmlParse(file = data)
    df.root <- XML::xmlRoot(df)

    # gene motif data
    gene_motif <- df.root[[4]] %>%
      XML::xmlToList() %>% unlist(recursive = FALSE)

    gene_motif <- do.call(rbind, lapply(gene_motif, data.frame))
    gene_motif$V2 <- rownames(gene_motif)
    colnames(gene_motif) <- c("V1", "V2")
    gene_motif$V2 <- gsub("scanned_sites\\..*\\.", "", gene_motif$V2)
    gene_motif$V2 <- gsub("[0-9]*", "", gene_motif$V2)

    ## seq index
    motif_seq_start <- which(grepl("sequence_", gene_motif$V2))
    motif_seq_end <- which(grepl("num_sites", gene_motif$V2))

    gene_info <- NULL
    for (i in c(1: length(motif_seq_start))){
      if (i == 1){
        tmp_motif <- gene_motif[1: motif_seq_start[i]- 1, ] %>%
          group_by(V2) %>%
          mutate(index = row_number()) %>%
          pivot_wider(names_from = V2, values_from = V1) %>%
          select(-index)

        tmp_seq <- gene_motif[motif_seq_start[i]:motif_seq_end[i], ] %>%
          pivot_wider(names_from = V2, values_from = V1) %>%
          rename(seq_pvalue = pvalue)

        gene_info <- cbind(tmp_motif, tmp_seq)
      } else {
        tmp_motif <- gene_motif[(motif_seq_end[i-1] + 1): (motif_seq_start[i]- 1), ] %>%
          group_by(V2) %>%
          mutate(index = row_number()) %>%
          pivot_wider(names_from = V2, values_from = V1) %>%
          select(-index)

        tmp_seq <- gene_motif[motif_seq_start[i]:motif_seq_end[i], ] %>%
          pivot_wider(names_from = V2, values_from = V1) %>%
          rename(seq_pvalue = pvalue)

        tmp <- cbind(tmp_motif, tmp_seq)
        gene_info <- rbind(gene_info, tmp)
      }
    }

    tmp <- df.root[[1]] %>% XML::xmlToList()
    seq_info <- tmp[2:(XML::xmlSize(tmp) - 1)] %>% bind_rows() %>%
      select(-alphabet_array) %>%
      na.omit() %>%
      rename(sequence_id = id) %>%
      rename(seq_name = name)


    motif_info <- df.root[[3]] %>%
      XML::xmlToList() %>% sapply("[[", ".attrs") %>%
      t() %>% as.data.frame() %>%
      rename(motif_id = id)

    gene_info <- merge(gene_info, seq_info, by = 'sequence_id')

    gene_info <- merge(gene_info, motif_info, by = 'motif_id') %>%
      dplyr::select(!c("sequence_id", "num_sites", "weight",
                       "name", "alt", "sites", "elapsed_time"))

    return(gene_info)

  }
}



#' @name combinedMotifseq
#'
#' @description Combine all motif seq to 1.
#' @param data motif_seq from \code{getMotifFromMEME}.
#' @return A single line of all motifs.
#'
#' @examples
#' data <- paste0(file, "meme.txt")
#' motif_seq <- getMotifFromMEME(data = data, format = "txt")
#' combined_seq <- combinedMotifseq(motif_seq)
#'
#'
#' @export
#'

combinedMotifseq <- function(data){
  df <- select(data, c("motif_seq", "motif_num")) %>%
    group_by(motif_num) %>%
    mutate(index = row_number()) %>%
    pivot_wider(names_from = 'motif_num', values_from = 'motif_seq') %>%
    select(-index)

  combineSeq <- function(seq_cols){
    seq_len = str_count(seq_cols)[!is.na(seq_cols)][1]
    ifelse(is.na(seq_cols),
           seq_cols[is.na(seq_cols)] <- stringr::str_pad("-", seq_len, "left", pad = "-"),
           seq_cols)}

  combinedSeq <- as.data.frame(apply(df, 2, function(x) combineSeq(x))) %>%
    unite(col, sep = "", remove = TRUE)

  return(combinedSeq)
}




