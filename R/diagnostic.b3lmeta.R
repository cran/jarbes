#' @title Diagnostic function for b3lmeta object in jarbes
#'
#' @description This function performers an approximated Bayesian cross-validation for a b3lmeta object
#'
#' @param object The object generated by the function b3lmeta.
#' @param post.p.value.cut Posterior p-value cut point to assess outliers.
#' @param study.names Character vector containing names of the studies used.
#' @param size.forest Size of the center symbol mark in the forest-plot lines
#' @param lwd.forest  Thickness of the lines in the forest-plot
#' @param shape.forest Type of symbol for the center mark in the forest-plot lines
#' @param ... \dots
#'
#'
#' @import ggplot2
#'
#' @export

diagnostic.b3lmeta = function(object,
                            post.p.value.cut = 0.05,
                            study.names = NULL,
                            size.forest = 0.4,
                            lwd.forest = 0.2,
                            shape.forest = 23,
                            ...) {

  x=y=ylo=yhi=NULL

  # Data preparation
  y.ghost = object$BUGSoutput$sims.list$y.ghost
  g.m = apply(y.ghost, 2, median)
  g.u = apply(y.ghost, 2, quantile, prob = 0.975)
  g.l = apply(y.ghost, 2, quantile, prob = 0.025)

  n.studies = length(g.m)

  TE = object$data$TE

  if (is.null(study.names)) {
    study.names = 1:n.studies
  }

  # Posterior p-values to detect outliers...
  p.vec = NULL
  for(i in 1:n.studies)
  {
    p1 = sum(y.ghost[,i]<TE[i])/length(y.ghost[,i])
    p2 = sum(y.ghost[,i]>TE[i])/length(y.ghost[,i])
    p.val = min(p1, p2)
    p.vec = c(p.vec, p.val)
  }

  p.col = ifelse(p.vec < post.p.value.cut, "red", "blue")

  data.plot = data.frame(
    x = study.names,
    TE = TE,
    g.m = g.m,
    ylo  = g.l,
    yhi  = g.u,
    p.vec = p.vec,
    p.col = p.col)


  p = ggplot(data.plot, aes(x = x, y = TE,
                            ymin = ylo, ymax = yhi)) +
    geom_pointrange(colour = p.col,
                    lwd = lwd.forest,        # Thickness of the lines
                    shape = shape.forest,
                    size = size.forest)+     # Point size
    coord_flip() +
    xlab("Study") +
    ylab("Posterior Predictive observation") +
    ggtitle("Bayesian Cross-Valdiation") +
    theme_bw()


  return(p)
}
