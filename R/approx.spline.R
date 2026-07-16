#' Highcharts-style cubic spline interpolation
#'
#' Drop-in replacement for \code{stats::approx}.
#' Uses the Catmull-Rom spline (with the same smoothing
#' parameter as Highcharts \code{type = "spline"}) instead of straight-line
#' segments.  The call signature and return value are identical to
#' \code{stats::approx}, so it can be substituted transparently.
#'
#' @param x Numeric vector of abscissa values.
#' @param y Numeric vector of ordinate values, same length as \code{x}.
#' @param xout Optional numeric vector of points where the spline is
#'   evaluated.  If \code{NULL}, a regular grid of length \code{n} is used.
#' @param n Integer. Number of points to generate when \code{xout} is
#'   \code{NULL}.  Default is 50.
#' @param rule Handling of points outside the range of \code{x}.  \code{1}
#'   returns \code{NA}; \code{2} repeats the nearest endpoint value.  Matches
#'   the semantics of \code{stats::approx}.
#' @param log Logical. If \code{TRUE}, interpolation is performed in
#'   log-log space.  Default is \code{FALSE}.
#' @param smoothing Numeric tension parameter used in the Catmull-Rom to
#'   Bezier conversion.  Default is 1.5 (Highcharts default).
#' @param ... Additional arguments ignored; present only for full signature
#'   compatibility with \code{stats::approx}.
#'
#' @return A list with components \code{x} and \code{y}, just like
#'   \code{stats::approx}.
#'
#' @seealso \code{\link[stats]{approx}}
#'
#' @examples
#' x <- c(0.05, 0.10, 0.20, 0.30, 0.50, 1.0, 2.0, 3.0)
#' y <- exp(-x)
#' Td <- seq(0.05, 3.0, length.out = 100)
#' approx.spline(x, y, xout = Td, log = TRUE)$y
#'
#' @export
approx.spline <- function(x, y = NULL,
                          xout = NULL,
                          n = 50,
                          rule = 1,
                          log = FALSE,
                          smoothing = 1.5,
                          ...) {
    if (is.null(y)) stop("y must be supplied")
    stopifnot(
        length(x) == length(y),
        length(x) >= 2,
        is.numeric(x), is.numeric(y)
    )

    ord <- order(x)
    x <- x[ord]
    y <- y[ord]

    ## -- handle log-log option ----------------------------------------------
    HasZero <- FALSE
    if (log) {
        if (x[1L] == 0) { # drop leading Tn = 0, keep to re-insert
            HasZero <- TRUE
            x0 <- x[1L]
            y0 <- y[1L]
            x <- x[-1L]
            y <- y[-1L]
        }
        if (any(x <= 0) || any(y <= 0)) {
            stop("log = TRUE requires strictly positive x and y, except for a single leading Tn = 0.")
        }
        tx <- log(x)
        ty <- log(y)
        toX <- exp
        toY <- exp
    } else {
        tx <- x
        ty <- y
        toX <- identity
        toY <- identity
    }

    npts <- length(tx)
    if (npts < 2L) {
        stop("Need at least two positive-period points to spline.")
    }

    ## -- Bezier control-point helper ----------------------------------------
    controls <- function(i) {
        p <- max(1, i - 1)
        nx <- min(npts, i + 1)

        dx1 <- tx[i] - tx[p]
        dy1 <- ty[i] - ty[p]
        dx2 <- tx[nx] - tx[i]
        dy2 <- ty[nx] - ty[i]
        denom <- abs(dx1) + abs(dx2)
        if (denom == 0) denom <- 1

        fL <- (abs(dx1) / denom) / smoothing
        fR <- (abs(dx2) / denom) / smoothing

        list(
            P0 = c(tx[i], ty[i]),
            P1 = c(tx[i] - fL * dx1, ty[i] - fL * dy1),
            P2 = c(tx[i] + fR * dx2, ty[i] + fR * dy2),
            P3 = c(tx[nx], ty[nx])
        )
    }

    ## -- choose evaluation abscissae ----------------------------------------
    if (!is.null(xout)) {
        xo <- if (log) ifelse(xout <= 0, NaN, log(xout)) else xout
    } else {
        xo <- seq(min(tx), max(tx), length.out = n)
    }

    ## -- locate interval for each xo ----------------------------------------
    idx <- findInterval(xo, tx, all.inside = FALSE)
    yout <- rep(NA_real_, length(xo))

    if (rule == 2L) {
        yout[idx == 0] <- ty[1L]
        yout[idx >= npts] <- ty[npts]
    }

    inside <- which(idx > 0 & idx < npts)
    if (length(inside)) {
        for (k in inside) {
            i <- idx[k]
            C <- controls(i)
            t <- (xo[k] - C$P0[1]) / (C$P3[1] - C$P0[1])
            yout[k] <- (1 - t)^3 * C$P0[2] +
                3 * (1 - t)^2 * t * C$P1[2] +
                3 * (1 - t) * t^2 * C$P2[2] +
                t^3 * C$P3[2]
        }
    }

    yout <- toY(yout)
    if (HasZero) { # splice (0, Sa0) back in
        xo <- c(x0, xo)
        yout <- c(y0, yout)
    }

    list(
        x = if (log) toX(xo) else xo,
        y = yout
    )
}
