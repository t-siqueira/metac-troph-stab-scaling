#' Split Violin plot
#'
#' This function is to plot the split violin graph.
#' It is built on [ggplot2::geom_violin()],
#' and adopted from jan-glx's stackoverflow answer \url{https://stackoverflow.com/a/45614547}
#' A violin plot is a compact display of a continuous distribution.
#' It is a blend of [ggplot2::geom_boxplot()] and [ggplot2::geom_density()]:
#' A violin plot is a mirrored density plot displayed in the same way as a boxplot.
#' A split violin plot allows for odd violins to be plotted to the left,
#' and even violins to be plotted to the right
#'
#' @param mapping Set of aesthetic mappings created by `aes()` or `aes_()`.
#' If specified and `inherit.aes = TRUE` (the default),
#' it is combined with the default mapping at the top level of the plot.
#' You must supply mapping if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three options:
#'
#' If `NULL`, the default, the data is inherited from the plot data
#' as specified in the call to [ggplot2::ggplot()].
#'
#' A `data.frame`, or other object, will override the plot data.
#' All objects will be fortified to produce a data frame.
#' See [ggplot2::fortify()] for which variables will be created.
#'
#' A function will be called with a single argument, the plot data.
#' The return value must be a `data.frame`, and will be used as the layer data.
#' A function can be created from a formula (e.g. `~ head(.x, 10)`).
#' @param stat Use to override the default connection between `geom_density` and `stat_density`.
#' @param position Position adjustment, either as a string,
#' or the result of a call to a position adjustment function.
#' @param ... Other arguments passed on to `layer()`.
#' These are often aesthetics, used to set an aesthetic to a fixed value,
#'  like `colour = "red"` or `size = 3`.
#' They may also be parameters to the paired geom/stat.
#' @param draw_quantiles If `not(NULL)` (default), draw horizontal lines
#' at the given quantiles of the density estimate.
#' @param trim If `TRUE `(default), trim the tails of the violins to the range of the data.
#' If `FALSE`, don't trim the tails.
#' @param scale If `"area"` (default), all violins have the same area (before trimming the tails).
#' If `"count"`, areas are scaled proportionally to the number of observations.
#' If `"width"`, all violins have the same maximum width.
#'
#' @param na.rm If `FALSE`, the default, missing values are removed with a warning.
#' If `TRUE`, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends?
#' `NA`, the default, includes if any aesthetics are mapped.
#' `FALSE` never includes, and `TRUE` always includes.
#' It can also be a named logical vector to finely select the aesthetics to display.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics,
#' rather than combining with them.
#' This is most useful for helper functions that define both data and aesthetics
#' and shouldn't inherit behaviour from the default plot specification, e.g. [ggplot2::borders()].
#'
#' @export
geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      trim = trim, scale = scale, draw_quantiles = draw_quantiles,
      na.rm = na.rm, ...
    )
  )
}

#' @format NULL
#' @usage NULL
GeomSplitViolin <- ggplot2:::ggproto("GeomSplitViolin",
                                     ggplot2::GeomViolin,
                                     draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                       data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                       grp <- data[1, "group"]
                                       newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                       newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                       newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                                       if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                         stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                                   1))
                                         quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                         aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                         aesthetics$alpha <- rep(1, nrow(quantiles))
                                         both <- cbind(quantiles, aesthetics)
                                         quantile_grob <- GeomPath$draw_panel(both, ...)
                                         ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                       }
                                       else {
                                         ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                       }
                                     }
)
