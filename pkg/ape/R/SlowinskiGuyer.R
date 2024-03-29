## SlowinskiGuyer.R (2011-03-10)

##   Tests of Diversification Shifts with Sister-Clades

## Copyright 2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

slowinskiguyer.test <- function(x, detail = FALSE)
{
    r <- x[, 1]
    n <- x[, 1] + x[, 2]
    pp <- (n - r)/(n - 1)
    chi <- -2 * sum(log(pp))
    df <- as.integer(2 * length(pp))
    pval <- pchisq(chi, df, lower.tail = FALSE)
    res <- data.frame("chisq" = chi, "df" = df, "P.val" = pval, row.names = "")
    if (detail)
        res <- list(res, individual_Pvalues = pp)
    res
}

mcconwaysims.test <- function(x)
{
    LRTp <- function(x) {
        f <- function(x) ifelse(x == 0, 0, x * log(x))
        n1 <- x[1]
        n2 <- x[2]
        1.629*(f(n1 - 1) - f(n1) + f(n2 - 1) - f(n2) - f(2) -
               f(n1 + n2 - 2) + f(n1 + n2))
    }
    chi <- sum(apply(x, 1, LRTp))
    pval <- pchisq(chi, df <- nrow(x), lower.tail = FALSE)
    data.frame("chisq" = chi, "df" = df, "P.val" = pval, row.names = "")
}

richness.yule.test <- function(x, t)
{
    n1 <- x[, 1]
    n2 <- x[, 2]
    n <- c(n1, n2)
    tb <- c(t, t)

    ## taken from /home/paradis/data/Projets/FISHLOSS/BDfit.R:
    .PrNt.Yule <- function(N, age, birth) {
        tmp <- exp(-birth * age)
        tmp * (1 - tmp)^(N - 1)
    }

    minusloglik0 <- function(l)
        -sum(log(.PrNt.Yule(n, tb, l)))

    minusloglika <- function(l)
        -sum(log(.PrNt.Yule(n1, t, l[1]))) - sum(log(.PrNt.Yule(n2, t, l[2])))

    out0 <- nlminb(0.01, minusloglik0, lower = 0, upper = 1)
    outa <- nlminb(rep(out0$par, 2), minusloglika,
                   lower = c(0, 0), upper = c(1, 1))
    chi <- 2 * (out0$objective - outa$objective)
    pval <- pchisq(chi, 1, lower.tail = FALSE)
    data.frame("chisq" = chi, "df" = 1, "P.val" = pval, row.names = "")
}

diversity.contrast.test <-
    function(x, method = "ratiolog", alternative = "two.sided", nrep = 0, ...)
{
    method <- match.arg(method, c("ratiolog", "proportion", "difference"))
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

    minmax <- t(apply(x, 1, sort)) # sort all rows
    DIFF <- x[, 1] - x[, 2]
    SIGN <- sign(DIFF)

    CONTRAST <- switch(method,
                       "ratiolog" = {
                           if (any(minmax == 1))
                               minmax <- minmax + 1 # prevent division by 0
                           ## Note: if min = max, no need to set the contrast
                           ## to zero since this is done with sign()
                           log(minmax[, 2]) / log(minmax[, 1])
                       },
                       "proportion" = minmax[, 2] / (minmax[, 2] + minmax[, 1]),
                       "difference" = abs(DIFF))

    y <- SIGN * CONTRAST # the signed contrasts

    if (nrep) {
        n <- length(SIGN)
        RND <-
            replicate(nrep,
                      sum(sample(c(-1, 1), size = n, replace = TRUE) * CONTRAST))
        cases <- switch(alternative,
                        "two.sided" = sum(abs(RND) > sum(y)),
                        "less" = sum(RND < sum(y)),
                        "greater" = sum(RND > sum(y)))
        cases/nrep
    } else wilcox.test(x = y, alternative = alternative, ...)$p.value
}
