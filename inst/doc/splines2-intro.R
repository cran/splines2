## ----setup, echo = FALSE------------------------------------------------------
knitr::opts_knit$set(global.par = TRUE)
knitr::opts_chunk$set(fig.width = 7, fig.height = 4)

## ----set-par, echo = FALSE----------------------------------------------------
library(graphics)
par(mar = c(2.5, 2.5, 0.5, 0.1), mgp = c(1.5, 0.5, 0))

## ----load-lib-----------------------------------------------------------------
library(splines2)
packageVersion("splines2")

## ----bSpline, fig.cap="B-splines of degree one with three internal knots placed at 0.3, 0.5, and 0.6."----
knots <- c(0.3, 0.5, 0.6)
x <- seq(0, 1, 0.01)
bsMat <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)
plot(bsMat, mark_knots = "all")

## ----ibs, fig.cap="Piecewise linear B-splines (left) and their integrals (right)."----
ibsMat <- ibs(x, knots = knots, degree = 1, intercept = TRUE)
op <- par(mfrow = c(1, 2))
plot(bsMat, mark_knots = "internal")
plot(ibsMat, mark_knots = "internal")
abline(h = c(0.15, 0.2, 0.25), lty = 2, col = "gray")

## ----dbs, fig.cap="Cubic B-spline (left) and their first derivative (right)."----
bsMat <- bSpline(x, knots = knots, intercept = TRUE)
dbsMat <- dbs(x, knots = knots, intercept = TRUE)
plot(bsMat, mark_knots = "internal")
plot(dbsMat, mark_knots = "internal")

## ----dbsMat-------------------------------------------------------------------
is_equivalent <- function(a, b) {
    all.equal(a, b, check.attributes = FALSE)
}
stopifnot(is_equivalent(dbsMat, deriv(bsMat)))

## ----pbs----------------------------------------------------------------------
px <- seq(0, 3, 0.01)
pbsMat <- bSpline(px, knots = knots, Boundary.knots = c(0, 1),
                  intercept = TRUE, periodic = TRUE)
ipMat <- ibs(px, knots = knots, Boundary.knots = c(0, 1),
             intercept = TRUE, periodic = TRUE)
dp1Mat <- deriv(pbsMat)
dp2Mat <- deriv(pbsMat, derivs = 2)
par(mfrow = c(1, 2))
plot(pbsMat, ylab = "Periodic B-splines", mark_knots = "boundary")
plot(ipMat, ylab = "The integrals", mark_knots = "boundary")
plot(dp1Mat, ylab = "The 1st derivatives", mark_knots = "boundary")
plot(dp2Mat, ylab = "The 2nd derivatives", mark_knots = "boundary")

## ----mSpline, fig.cap = "Quadratic M-spline with three internal knots placed at 0.3, 0.5, and 0.6."----
msMat <- mSpline(x, knots = knots, degree = 2, intercept = TRUE)
par(op)
plot(msMat, mark_knots = "all")

## ----mSpline-derivs-----------------------------------------------------------
dmsMat1 <- mSpline(x, knots = knots, degree = 2, intercept = TRUE, derivs = 1)
dmsMat2 <- deriv(msMat)
stopifnot(is_equivalent(dmsMat1, dmsMat2))

## ----pms-basis, fig.cap = "Cubic periodic M-splines."-------------------------
pmsMat <- mSpline(px, knots = knots, intercept = TRUE,
                  periodic = TRUE, Boundary.knots = c(0, 1))
plot(pmsMat, ylab = "Periodic Basis", mark_knots = "boundary")

## ----pms-deriv, fig.cap = "The first derivatives of the periodic M-splines."----
dpmsMat <- deriv(pmsMat)
plot(dpmsMat, ylab = "The 1st derivatives", mark_knots = "boundary")

## ----pms-integral, fig.cap = "The integrals of the periodic M-splines."-------
ipmsMat <- mSpline(px, knots = knots, intercept = TRUE,
                   periodic = TRUE, Boundary.knots = c(0, 1), integral = TRUE)
plot(ipmsMat, ylab = "Integrals", mark_knots = "boundary")
abline(h = seq.int(0, 3), lty = 2, col = "gray")

## ----iSpline, fig.cap = "I-splines of degree two with three internal knots placed at 0.3, 0.5, and 0.6."----
isMat <- iSpline(x, knots = knots, degree = 2, intercept = TRUE)
plot(isMat, mark_knots = "internal")

## ----msMat--------------------------------------------------------------------
stopifnot(is_equivalent(msMat, deriv(isMat)))

## ----dmsMat-------------------------------------------------------------------
dmsMat3 <- deriv(isMat, 2)
stopifnot(is_equivalent(dmsMat1, dmsMat3))

## ----cSpline-scaled, fig.cap = "C-splines of degree two with three internal knots placed at 0.3, 0.5, and 0.6."----
csMat1 <- cSpline(x, knots = knots, degree = 2, intercept = TRUE)
plot(csMat1)
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
plot(bpMat1)
plot(bpMat2)

## ----bp-2, fig.height=6, fig.cap = "The integrals (upper panel) and the first derivatives (lower panel) of Bernstein polynomials of degree 4."----
ibpMat1 <- bernsteinPoly(x1, degree = 4, intercept = TRUE, integral = TRUE)
ibpMat2 <- bernsteinPoly(x2, degree = 4, intercept = TRUE, integral = TRUE)
dbpMat1 <- bernsteinPoly(x1, degree = 4, intercept = TRUE, derivs = 1)
dbpMat2 <- bernsteinPoly(x2, degree = 4, intercept = TRUE, derivs = 1)
par(mfrow = c(2, 2))
plot(ibpMat1, ylab = "Integrals")
plot(ibpMat2, ylab = "Integrals")
plot(dbpMat1, ylab = "Derivatives")
plot(dbpMat2, ylab = "Derivatives")

## ----bp-deriv-----------------------------------------------------------------
stopifnot(is_equivalent(dbpMat1, deriv(bpMat1)))
stopifnot(is_equivalent(dbpMat2, deriv(bpMat2)))
stopifnot(is_equivalent(dbpMat1, deriv(ibpMat1, 2)))
stopifnot(is_equivalent(dbpMat2, deriv(ibpMat2, 2)))

## ----ns-basis, fig.cap = "Nonnegative natural cubic splines (left) and corresponding integrals (right)."----
nsMat <- naturalSpline(x, knots = knots, intercept = TRUE)
insMat <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
par(mfrow = c(1, 2))
plot(nsMat, ylab = "Basis")
plot(insMat, ylab = "Integrals")
stopifnot(is_equivalent(nsMat, deriv(insMat)))

## ----ns-deriv, fig.cap = "The derivatives of natural cubic splines."----------
d1nsMat <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
d2nsMat <- deriv(nsMat, 2)
matplot(x, d1nsMat, type = "l", ylab = "The 1st derivatives")
matplot(x, d2nsMat, type = "l", ylab = "The 2nd derivatives")

## ----nsk----------------------------------------------------------------------
nskMat <- nsk(x, knots = knots, intercept = TRUE)
par(op)
plot(nskMat, ylab = "nsk()", mark_knots = "all")
abline(h = 1, col = "red", lty = 3)

## ----update-nsk---------------------------------------------------------------
nskMat2 <- update(nskMat, knots = c(knots, 0.2, 0.4), trim = 0.025)
knots(nskMat2)
stopifnot(all.equal(quantile(x, c(0.025, 0.975), names = FALSE),
                    knots(nskMat2, "boundary")))

## ----predict------------------------------------------------------------------
new_x <- c(0.275, 0.525, 0.8)
names(new_x) <- paste0("x=", new_x)
(isMat2 <- predict(isMat, newx = new_x)) # the basis functions at new x
stopifnot(all.equal(predict(isMat, newx = new_x), update(isMat, x = new_x)))
## compute the first derivative of the I-spline function in different ways
beta <- seq(0.1, by = 0.1, length.out = ncol(isMat))
deriv_ispline1 <- predict(isMat, newx = new_x, coef = beta, derivs = 1)
deriv_ispline2 <- predict(update(isMat, x = new_x, derivs = 1), coef = beta)
deriv_ispline3 <- c(predict(deriv(isMat), newx = new_x) %*% beta)
stopifnot(all.equal(deriv_ispline1, deriv_ispline2))
stopifnot(all.equal(deriv_ispline2, deriv_ispline3))

## ----plot-coef----------------------------------------------------------------
beta <- seq.int(0.2, length.out = ncol(nskMat), by = 0.2)
plot(nskMat, ylab = "nsk()", mark_knots = "all", coef = beta)
abline(h = beta, col = seq_along(beta), lty = 3)

## ----formula-alias------------------------------------------------------------
b <- bSpline # create an alias for B-splines
mod1 <- lm(weight ~ b(height, degree = 1, df = 3), data = women)
iknots <- with(women, knots(bSpline(height, degree = 1, df = 3)))
mod2 <- lm(weight ~ bSpline(height, degree = 1, knots = iknots), data = women)
pred1 <- predict(mod1, head(women, 10))
pred2 <- predict(mod2, head(women, 10))
stopifnot(all.equal(pred1, pred2))

## ----formula-wrap-failed------------------------------------------------------
## generates quadratic spline basis functions based on log(x)
qbs <- function(x, ...) {
    splines2::bSpline(log(x), ..., degree = 2)
}
mod3 <- lm(weight ~ qbs(height, df = 5), data = women)
mod4 <- lm(weight ~ bsp(log(height), degree = 2, df = 5), data = women)
stopifnot(all.equal(unname(coef(mod3)), unname(coef(mod4)))) # the same coef
pred3 <- predict(mod3, head(women, 10))
pred4 <- predict(mod4, head(women, 10))
all.equal(pred3, pred4)
pred0 <- predict(qbs(women$height, df = 5),
                 newx = head(log(women$height), 10),
                 coef = coef(mod3)[- 1]) + coef(mod3)[1]
stopifnot(all.equal(pred0, pred4, check.names = FALSE))

## ----predict-qbs--------------------------------------------------------------
## generates quadratic spline basis functions based on log(x) with a new class
qbs <- function(x, ...) {
    res <- splines2::bSpline(log(x), ..., degree = 2)
    class(res) <- c("qbs", class(res))
    return(res)
}
## a utility to help model.frame() create the right matrices
makepredictcall.qbs <- function(var, call) {
    if (as.character(call)[1L] == "qbs" ||
    (is.call(call) && identical(eval(call[[1L]]), qbs))) {
        at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                                "periodic", "derivs", "integral")]
        call <- call[1L:2L]
        call[names(at)] <- at
    }
    call
}
## the same example
mod3 <- lm(weight ~ qbs(height, df = 5), data = women)
mod4 <- lm(weight ~ bsp(log(height), degree = 2, df = 5), data = women)
stopifnot(all.equal(unname(coef(mod3)), unname(coef(mod4)))) # the same coef
pred3 <- predict(mod3, head(women, 10))
pred4 <- predict(mod4, head(women, 10))
all.equal(pred3, pred4) # should be TRUE this time

## ----extract------------------------------------------------------------------
c(nskMat2$trim, attr(nskMat2, "trim"))

