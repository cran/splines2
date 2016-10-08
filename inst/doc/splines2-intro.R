## ----mSpline, fig.width=7, fig.height=4, fig.cap = "Quadratic M-splines with three internal knots."----
library(splines2)
knots <- c(0.3, 0.5, 0.6)
x <- seq(0, 1, 0.01)
msOut <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
library(graphics) # attach graphics just in case for plots
par(mar = c(2.5, 2.5, 0, 0), mgp = c(1.5, 0.5, 0))
matplot(x, msOut, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray") # mark internal knots

## ----iSpline, fig.width=7, fig.height=4, fig.cap = "I-splines of degree two with three internal knots."----
isOut <- iSpline(x, knots = knots, degree = 2, intercept = TRUE)
par(mar = c(2.5, 2.5, 0, 0), mgp = c(1.5, 0.5, 0))
matplot(x, isOut, type = "l", ylab = "y")
abline(h = 1, v = knots, lty = 2, col = "gray")

## ----msMat---------------------------------------------------------------
all.equal(msOut, attr(isOut, "msMat")) # 'msMat' attribute for M-splines

## ----cSpline-scaled, fig.width=7, fig.height=4, fig.cap = "C-splines of degree two with three internal knots."----
csOut1 <- cSpline(x, knots = knots, degree = 2, intercept = TRUE)
par(mar = c(2.5, 2.5, 0, 0), mgp = c(1.5, 0.5, 0))
matplot(x, csOut1, type = "l", ylab = "y")
abline(h = 1, v = knots, lty = 2, col = "gray")

## ----ibs, fig.width=7, fig.height=4, fig.cap="Piecewise linear B-splines (left) and their integrals (right)."----
ibsOut <- ibs(x, knots = knots, degree = 1, intercept = TRUE)
par(mar = c(2.5, 2.5, 0, 0), mgp = c(1.5, 0.5, 0), mfrow = c(1, 2))
matplot(x, attr(ibsOut, "bsMat"), type = "l", ylab = "y")
abline(v = knots, h = 1, lty = 2, col = "gray")
matplot(x, ibsOut, type = "l", ylab = "y")
abline(v = knots, h = c(0.15, 0.2, 0.25), lty = 2, col = "gray")

## ----bSpline, fig.width=7, fig.height=4, fig.cap="B-splines (left) and M-splines (right) of degree zero"----
bsOut0 <- bSpline(x, knots = knots, degree = 0, intercept = TRUE)
msOut0 <- mSpline(x, knots = knots, degree = 0, intercept = TRUE)
par(mar = c(2.5, 2.5, 0, 0), mgp = c(1.5, 0.5, 0), mfrow = c(1, 2))
matplot(x, bsOut0, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")
matplot(x, msOut0, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")

## ----predict-------------------------------------------------------------
predict(isOut, c(0.275, 0.525, 0.8))
