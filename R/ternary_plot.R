#' clean data
#' @title data_clean
#' @description buildding a stanard table for ggtern
#' @param otu a tabbed text files in which OTUs are rows and samples are columns,
#' one entry in the table is usually a number of reads, or a frequency in the range 0.0 to 1.0.
#' @param design a data table in which fixing the header of sample as 'sampleID' and the header of group as 'Group', respectively.
#' @param type the type of otu table:  relative/absolute abundance otu-table
#' @param times the zoom scale of the points.
#' @export
#' @return cleaned data
#' @author Yilei
#' @examples
#' file_path <- system.file("extdata", package = "aggregatesR" )
#' otutab <- read.delim(paste(file_path, "otutab.txt", sep = "/"), header=T, row.names=1)
#' design <- read.delim(paste(file_path, "metadata.txt", sep = "/"), header=T, row.names=NULL)
#' design = design[ , c("SampleID","Group")]
#' otu_tern <- data_clean(otutab, design, type="absolute", threshold=0.001, times=100)
#' head(otu_tern,n=3)
#' p <- ggtern(data = otu_tern, aes(x = KO, y = OE, z = WT)) +
#' geom_point(aes(size = size), alpha = 0.8, show.legend = T) +
#' scale_size(range = c(0, 6)) + geom_mask() +
#' guides(colour = "none") + theme_bw() +
#' theme(axis.text = element_blank(), axis.ticks = element_blank())

data_clean <- function(otu, design, type = c("relative", "absolute"), threshold = 0.001, times = 100){

  if (type == "absolute"){
    otu_relative <- apply(otu, 2, function(x){x/sum(x)})
  }else {otu_relative <- otu}

  idx <- rowSums(otu_relative > threshold) >= 1
  otu_threshold <- as.data.frame(otu_relative[idx, ])
  otu_threshold$OTUs <- row.names(otu_threshold)


  otu_longer <- pivot_longer(data = otu_threshold,
                             cols = -OTUs,
                             names_to = "SampleID",
                             values_to = "value")


  merge_data <- merge(otu_longer, design, by ="SampleID")

  # otu <- subset(merge_data, select = -SampleID)
  otu <- subset(merge_data, select = c("Group", "OTUs", "value"))

  otu_mean <- otu %>% group_by(OTUs, Group) %>%
    summarise(value = mean(value))

  otu_tern <- otu_mean %>%
    group_by(Group, OTUs) %>%
    mutate(index = row_number()) %>%
    pivot_wider(names_from = Group,values_from = value) %>%
    select(-index)

  otu_tern$size <- (apply(otu_tern[2:4], 1, mean))*times
  return(otu_tern)
}


#' enrich data
#' @title enrich_data
#' @description the differential analysis of otu table
#' @param otu a tabbed text files in which OTUs are rows and samples are columns,
#' one entry in the table is a number of reads.
#' @param design a data table in which fixing the header of sample as 'sampleID' and the header of group as 'Group', respectively.
#' @export
#' @return cleaned data
#' @author Yilei
#' @examples
#' enrich_index <- enrich_data(otutab, design, p.value = 0.05)
#' plot_data <- merge(otu_tern, enrich_index, by = "OTUs", all.x = T)
#' p <- ggtern(data = plot_data, aes(x = KO, y = OE, z = WT)) +
#' geom_mask() +
#' geom_point(aes(size=size, color=enrich),alpha=0.8) +
#' guides(size="none") +theme_bw() +
#' theme(axis.text=element_blank(),
#' axis.ticks=element_blank())
#' @importFrom edgeR DGEList
#'

enrich_data <- function(otu, design, p.value = 0.05, adjust.method = "fdr"){
  dge_list <- DGEList(counts=otu, group=design$Group)
  # Remove the lower abundance/(cpm, rpkm)
  keep <- rowSums(dge_list$counts) >= 0
  dge_keep <- dge_list[keep, ,keep.lib.sizes=F]
  # scale the raw library sizes dgelist
  dge <- calcNormFactors(dge_keep)
  # fit the GLM
  design.mat <- model.matrix(~ 0 + dge$samples$group)
  d2 <- estimateGLMCommonDisp(dge, design.mat)
  d2 <- estimateGLMTagwiseDisp(d2, design.mat)
  fit <- glmFit(d2, design.mat)
  #######
  # if (missing(adjust.method))
  #   adjust.method="fdr"
  # if (missing(p.value))
  #   p.value=0.05
  group_index <- as.character(design$Group[!duplicated(design$Group)])
  # enrich groups
  lrt_1_2 <- glmLRT(fit, contrast=c(1, -1, 0))
  lrt_1_3 <- glmLRT(fit, contrast=c(1, 0, -1))
  de_1_2 <- decideTestsDGE(lrt_1_2, adjust.method=adjust.method,
                           p.value=p.value)
  de_1_3 <- decideTestsDGE(lrt_1_3, adjust.method=adjust.method,
                           p.value=p.value)

  rich_1 <- rownames(otu)[de_1_2 == 1 & de_1_3 == 1]
  enrich_1 <- data.frame(OTUs=rich_1,
                         enrich=rep(group_index[1], length(rich_1)))
  ###############################
  lrt_2_3 <- glmLRT(fit, contrast=c(0, 1, -1))
  lrt_2_1 <- glmLRT(fit, contrast=c(-1, 1, 0))

  de_2_3 <- decideTestsDGE(lrt_2_3, adjust.method=adjust.method,
                           p.value=p.value)
  de_2_1 <- decideTestsDGE(lrt_2_1, adjust.method=adjust.method,
                           p.value=p.value)

  rich_2 <- rownames(otu)[de_2_3 == 1 & de_2_1 == 1]
  enrich_2 <- data.frame(OTUs=rich_2,
                         enrich=rep(group_index[2], length(rich_2)))
  ###################
  lrt_3_1 <- glmLRT(fit, contrast=c(-1, 0, 1))
  lrt_3_2 <- glmLRT(fit, contrast=c(0, -1, 1))

  de_3_1 <- decideTestsDGE(lrt_3_1, adjust.method=adjust.method,
                           p.value=p.value)
  de_3_2 <- decideTestsDGE(lrt_3_2, adjust.method=adjust.method,
                           p.value=p.value)

  rich_3 <- rownames(otu)[de_3_1 == 1 & de_3_2 == 1]
  enrich_3 <- data.frame(OTUs=rich_3,
                         enrich=rep(group_index[3], length(rich_3)))
  enrich_index <- rbind(enrich_1, enrich_2, enrich_3)
  return(enrich_index)
}


