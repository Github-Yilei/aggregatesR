#' D.stat
#'
#' @description The function D.stat is for computing the ABBA and BABA proportions at each site of a 'freq.py' table and use these to compute the D atstistic.
#' @description Further details are available on: \url{https://github.com/simonhmartin/tutorials/blob/master/ABBA_BABA_whole_genome/README.md}
#' @description freq.py: \url{https://github.com/simonhmartin/genomics_general}
#' @parm p populations.
#' @return a D statistic (D varies from -1 to 1). D > 0 indicative of excess shared ancestry between P2 and P3.
#'
#' @examples
#' P1 <- "LR"
#' P2 <- "KQ"
#' P3 <- "MD"
#' D <- D.stat(freq[,P1], freq[,P2], freq[,P3])
#' print(paste("D =", round(D,4)))
#'
#' @export
#'

D.stat <- function(p1, p2, p3) {
    ABBA <- (1 - p1) * p2 * p3
    BABA <- p1 * (1 - p2) * p3
    (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
}


#' get_block_indices
#' @description The function get_block_indices in the jackknife script will  the blocks, and return the 'indices' (i.e. the rows in our frequencies table) corresponding to each block.
#' @param block_size the block size along with chromosome
#' @param positions the position for each site to be analyzed.
#' @param chromosomes the chromosomes.
#'
#' @return a jackknife block indices
#'
#' @examples
#' block_indices <- get_block_indices(block_size=1e6,
#' positions = freq_table$position,
#' chromosomes = freq_table$scaffold)
#' n_blocks <- length(block_indices)
#' print(paste("Genome divided into", n_blocks, "blocks."))
#'
#' @export
#'

get_block_indices <- function(block_size, positions, chromosomes = NULL){
    if (is.null(chromosomes) == TRUE) {
        block_starts <- seq(min(positions), max(positions), block_size)

        block_ends <- block_starts + block_size - 1

        lapply(1:length(block_starts), function(x) which(positions >= block_starts[x] &
                                                         positions <= block_ends[x]))
        }
    else {
        chrom_names <- unique(chromosomes)

        block_starts <- lapply(chrom_names, function(chrom_name) seq(min(positions[chromosomes==chrom_name]),
                                                                     max(positions[chromosomes==chrom_name]), block_size))

        block_chroms <- unlist(lapply(1:length(block_starts), function(x) rep(chrom_names[x], length(block_starts[[x]]))))

        block_starts <- unlist(block_starts)

        block_ends <- block_starts + block_size - 1

        lapply(1:length(block_starts), function(x) which(chromosomes == block_chroms[x] &
                                                     positions >= block_starts[x] &
                                                     positions <= block_ends[x]))
        }
    }

#' get_jackknife_sd
#'
#' @description this function runs the jackknife procedure by calculating pseudovalues by removing one block at a time
#' @param block_indices the blocks defined by get_block_indices.
#' @param FUN the D statistic function (D.stat) we created earlier
#' @return the standad deviation of the D statistic.
#'
#' @examples
#' D_sd <- get_jackknife_sd(block_indices = block_indices,
#' FUN=D.stat,
#' freq_table[,P1], freq_table[,P2], freq_table[,P3])
#'
#' @export
#'

get_jackknife_sd <- function(block_indices, FUN, ...){
    n_blocks <- length(block_indices)
    args = list(...)
    overall_mean <- FUN(...)
    if (is.null(dim(args[1])) == TRUE){
        return(sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]]]))*(n_blocks-1))))
        }
    else{
        return(sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]],]))*(n_blocks-1))))
        }
}

#' block_jackknife
#'
#' @export
#'

block_jackknife <- function(block_indices, FUN, ...){
    n_blocks <- length(block_indices)
    args = list(...)
    overall_mean <- FUN(...)
    if (is.null(dim(args[1])) == TRUE){
        pseudovalues <- sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]]]))*(n_blocks-1))
        }
    else{
        pseudovalues <- sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]],]))*(n_blocks-1))
        }

    mean <- mean(pseudovalues)

    std_dev <- sd(pseudovalues)

    list(mean=mean, std_dev=std_dev)
}

#' f.stat
#'
#' @description The D statistic provides a powerful test for introgression, but it does not quantify the proportion of the genome that has been shared. A related method has been developed to estimate f, the 'admixture proportion'.
#' @description The idea behind this approach is that we compare the observed excess of ABBA over BABA sites, to that which would be expected under complete admixture. To approximate the expectation under complete admixture we re-count ABBA and BABA but substituting a second population of the P3 species in the place of P2. If you lack a second population, you can simply split your P3 samples into two.
#' @examples
#' P1 <- "LR"
#' P2 <- "KQ"
#' P3a <- "MD1"
#' P3b <- "MD2"
#' f <- f.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3a], freq_table[,P3b])
#'
#' @export
#'

f.stat <- function(p1, p2, p3a, p3b) {
    ABBA_numerator <- (1 - p1) * p2 * p3a
    BABA_numerator <- p1 * (1 - p2) * p3a

    ABBA_denominator <- (1 - p1) * p3b * p3a
    BABA_denominator <- p1 * (1 - p3b) * p3a

    (sum(ABBA_numerator) - sum(BABA_numerator)) /
        (sum(ABBA_denominator) - sum(BABA_denominator))
}
