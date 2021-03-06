## ----setup, echo = FALSE------------------------------------------------------
knitr::opts_knit$set(global.par = TRUE)
knitr::opts_chunk$set(fig.width = 7, fig.height = 4)

## ----set-par, echo = FALSE----------------------------------------------------
library(graphics)
par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))

## ----bSpline, fig.cap="B-splines of degree one with two internal knots."------
library(splines2)
knots <- c(0.3, 0.5, 0.6)
x <- seq(0, 1, 0.01)
bsMat <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)
matplot(x, bsMat, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")

## ----ibs, fig.cap="Piecewise linear B-splines (left) and their integrals (right)."----
ibsMat <- ibs(x, knots = knots, degree = 1, intercept = TRUE)
par(mfrow = c(1, 2))
matplot(x, bsMat, type = "l", ylab = "y")
abline(v = knots, h = 1, lty = 2, col = "gray")
matplot(x, ibsMat, type = "l", ylab = "y")
abline(v = knots, h = c(0.15, 0.2, 0.25), lty = 2, col = "gray")

## ----dbs, fig.cap="Cubic B-spline (left) and their first derivative (right)."----
bsMat <- bSpline(x, knots = knots, intercept = TRUE)
dbsMat <- dbs(x, knots = knots, intercept = TRUE)
par(mfrow = c(1, 2))
matplot(x, bsMat, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")
matplot(x, dbsMat, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")

## ----dbsMat-------------------------------------------------------------------
is_equivalent <- function(a, b) {
    all.equal(a, b, check.attributes = FALSE)
}
stopifnot(is_equivalent(dbsMat, deriv(bsMat)))

## ----reset-par-mSpline, echo = FALSE------------------------------------------
par(mfrow = c(1, 1))

## ----mSpline, fig.cap = "Quadratic M-spline with three internal knots."-------
msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
matplot(x, msMat, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")

## ----mSpline-derivs-----------------------------------------------------------
dmsMat1 <- mSpline(x, knots = knots, degree = 2, intercept = TRUE, derivs = 1)
dmsMat2 <- deriv(msMat)
stopifnot(is_equivalent(dmsMat1, dmsMat2))

## ----pms-basis, fig.cap = "Cubic periodic M-splines."-------------------------
x1 <- seq.int(0, 3, 0.01)
pmsMat <- mSpline(x1, knots = knots, degree = 3, intercept = TRUE,
                  periodic = TRUE, Boundary.knots = c(0, 1))
matplot(x1, pmsMat, type = "l", xlab = "x", ylab = "Periodic Basis")
abline(v = seq.int(0, 3), lty = 2, col = "gray")

## ----pms-deriv, fig.cap = "The first derivatives of the periodic M-splines."----
dpmsMat <- deriv(pmsMat)
matplot(x1, dpmsMat, type = "l", xlab = "x", ylab = "The 1st derivatives")
abline(v = seq.int(0, 3), lty = 2, col = "gray")

## ----pms-integral, fig.cap = "The integrals of the periodic M-splines."-------
ipmsMat <- mSpline(x1, knots = knots, degree = 3, intercept = TRUE,
                   periodic = TRUE, Boundary.knots = c(0, 1), integral = TRUE)
matplot(x1, ipmsMat, type = "l", xlab = "x", ylab = "Integrals")
abline(v = seq.int(0, 3), h = seq.int(0, 3), lty = 2, col = "gray")

## ----iSpline, fig.cap = "I-splines of degree two with three internal knots."----
isMat <- iSpline(x, knots = knots, degree = 2, intercept = TRUE)
matplot(x, isMat, type = "l", ylab = "y")
abline(h = 1, v = knots, lty = 2, col = "gray")

## ----msMat--------------------------------------------------------------------
stopifnot(is_equivalent(msMat, deriv(isMat)))

## ----dmsMat-------------------------------------------------------------------
dmsMat3 <- deriv(isMat, 2)
stopifnot(is_equivalent(dmsMat1, dmsMat3))

## ----cSpline-scaled, fig.cap = "C-splines of degree two with three internal knots."----
csMat1 <- cSpline(x, knots = knots, degree = 2, intercept = TRUE)
matplot(x, csMat1, type = "l", ylab = "y")
abline(h = 1, v = knots, lty = 2, col = "gray")

## ----cSpline-not-scaled-------------------------------------------------------
csMat2 <- cSpline(x, knots = knots, degree = 2, intercept = TRUE, scale = FALSE)
stopifnot(is_equivalent(isMat, deriv(csMat2)))
stopifnot(is_equivalent(msMat, deriv(csMat2, 2)))
stopifnot(is_equivalent(msMat, deriv(deriv(csMat2))))

## ----bp-1, fig.cap = "Bernstein polynomials of degree 4 over [0, 1] (left) and the generalized version over [- 1, 1] (right)."----
x1 <- seq.int(0, 1, 0.01)
x2 <- seq.int(- 1, 1, 0.01)
bpMat1 <- bernsteinPoly(x1, degree = 4, intercept = TRUE)
bpMat2 <- bernsteinPoly(x2, degree = 4, intercept = TRUE)
par(mfrow = c(1, 2))
matplot(x1, bpMat1, type = "l", ylab = "y")
matplot(x2, bpMat2, type = "l", ylab = "y")

## ----bp-2, fig.height=6, fig.cap = "The integrals (upper panel) and the first derivatives (lower panel) of Bernstein polynomials of degree 4."----
ibpMat1 <- bernsteinPoly(x1, degree = 4, intercept = TRUE, integral = TRUE)
ibpMat2 <- bernsteinPoly(x2, degree = 4, intercept = TRUE, integral = TRUE)
dbpMat1 <- bernsteinPoly(x1, degree = 4, intercept = TRUE, derivs = 1)
dbpMat2 <- bernsteinPoly(x2, degree = 4, intercept = TRUE, derivs = 1)
par(mfrow = c(2, 2))
matplot(x1, ibpMat1, type = "l", ylab = "Integrals")
matplot(x2, ibpMat2, type = "l", ylab = "y")
matplot(x1, dbpMat1, type = "l", ylab = "y")
matplot(x2, dbpMat2, type = "l", ylab = "y")

## ----bp-deriv-----------------------------------------------------------------
stopifnot(is_equivalent(dbpMat1, deriv(bpMat1)))
stopifnot(is_equivalent(dbpMat2, deriv(bpMat2)))
stopifnot(is_equivalent(dbpMat1, deriv(ibpMat1, 2)))
stopifnot(is_equivalent(dbpMat2, deriv(ibpMat2, 2)))

## ----ns-basis, fig.cap = "Nonnegative natural cubic splines (left) and corresponding integrals (right)."----
nsMat <- naturalSpline(x, knots = knots, intercept = TRUE)
insMat <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
par(mfrow = c(1, 2))
matplot(x, nsMat, type = "l", ylab = "Basis")
matplot(x, insMat, type = "l", ylab = "Integrals")
stopifnot(is_equivalent(nsMat, deriv(insMat)))

## ----ns-deriv, fig.cap = "The derivatives of natural cubic splines."----------
d1nsMat <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
d2nsMat <- deriv(nsMat, 2)
par(mfrow = c(1, 2))
matplot(x, d1nsMat, type = "l", ylab = "The 1st derivatives")
matplot(x, d2nsMat, type = "l", ylab = "The 2nd derivatives")

## ----predict------------------------------------------------------------------
new_x <- c(0.275, 0.525, 0.8)
names(new_x) <- paste0("x=", new_x)
predict(isMat, new_x)

