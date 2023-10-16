#' @name Lasso2plot
#'
#' @title Extract and Visualize Lasso Information from glmnet package.
#' @description
#' \code{Lasso2plot} Extract and Return a ggplot2 list of Lasso Information from glmnet package.
#'
#' @param data A model from  from glmnet package.
#' @param format A list from Lasso or cv Lasso fit.
#'
#' @importFrom magrittr %>%
#' @importFrom XML xmlParse xmlRoot xmlSize xmlToList
#' @importFrom ggplot2 aes ggplot geom_segment geom_text theme_bw
#' @importFrom dplyr select filter mutate
#' @importFrom tidyr pivot_longer
#'
#' @return Return a ggplot2 list
#'
#' @examples
#' Not run:
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' mod <- glmnet(x, y)
#' mod <- cv.glmnet(x, y)
#' p <- Lasso2plot(mod, format = "cv")
## End(**Not run**)
#' @export


Lasso2plot <- function(mod, format = "cv"){
  if (format == "cv"){
    cvfit <- mod
    cvdf <- data.frame(lamda = cvfit[["lambda"]] %>% log() %>% round(.,3),
                      cvm = cvfit[["cvm"]],
                      cvup = cvfit[["cvup"]],
                      cvlo = cvfit[["cvlo"]],
                      nzero = cvfit[["nzero"]])

    lambda_min <- log(cvfit$lambda.min)
    lambda_1se <- log(cvfit$lambda.1se)

    xbreaks <-  cvdf$lamda[seq(1,100,10)]
    xlabels <- as.character(cvdf$nzero[seq(1,100,10)])

    p <- ggplot(data =  cvdf) + geom_errorbar(aes(x = lamda, ymin = cvlo, ymax = cvup,color = nzero),
                                            linewidth = .5, width = 0.1, alpha = 0.5) +
      geom_point(aes(x = lamda,y = cvm, color = nzero)) +
      scale_x_continuous(sec.axis = sec_axis(~., breaks = xbreaks, labels = xlabels)) +
      geom_text(aes(x = lambda_min + 0.1, y = 1.2, label = 'lambda_min'), size = 4)+
      geom_text(aes(x = lambda_1se + 0.1, y = 1.2, label = 'lambda.1se'), size = 4) +
      labs(x = expression("Log"["10"]*"(lambda)"), y = 'Binomial Deviance', color = 'vars') +
      geom_vline(xintercept = c(lambda_min, lambda_1se), lty = 2, col = 'grey50') +
      theme_bw() +
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 14)) +
      guides(colour = 'none')

    return(p)
  } else {
    mod <- mod
    gl_df <- as_tibble(as.matrix(coef(mod)), rownames = "coef") %>%
      tidyr::pivot_longer(cols = -coef,
                   names_to = "step",
                   names_transform = list(step = parse_number),
                   values_to = "estimate") %>%
      #mutate(step = step)  %>%   #  lambda not 0
      group_by(step) %>%
    # lambda not start at 0
      mutate(lambda = mod$lambda[step + 1 ],
             dev.ratio = mod$dev.ratio[step + 1])  %>% 
      filter(coef != '(Intercept)') 

    p <- ggplot(gl_df, aes(x = log(lambda), y = estimate, group = coef, color= coef)) +
      geom_line(size = 1.2)+
      guides(color = 'none') +
      geom_hline(yintercept = 0)+
      labs(x = expression("Log"["10"]*"(lambda)"), y = "Coefficients") +
      theme_bw() +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14))
    scale_y_continuous(sec.axis = sec_axis( ~rescale(.,c(0,0.5)),name = "Categroy",labels=sprintf("%d%%",(0:5)*10)))
    return(p)

  }
}

