library(splines2)

x <- seq.int(0, 1, 0.2)
knots <- c(0.3, 0.5, 0.6)
newx <- seq.int(0.1, 0.9, 0.2)

## Cubic B-spline basis functions
bs_mat <- bSpline(x, knots = knots)

## compute the B-spline basis functions at new x
predict(bs_mat, newx)

## compute the B-spline function for the specified coefficients
beta <- runif(ncol(bs_mat))
predict(bs_mat, coef = beta)

## compute the first derivative of the B-spline function
predict(bs_mat, coef = beta, derivs = 1)
## or equivalently
predict(deriv(bs_mat), coef = beta)

## compute the second derivative
predict(bs_mat, coef = beta, derivs = 2)
## or equivalently
predict(deriv(bs_mat, derivs = 2), coef = beta)

## compute the integral
predict(bs_mat, coef = beta, integral = TRUE)
## or equivalently
predict(update(bs_mat, integral = TRUE), coef = beta)

## visualize
op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.5, 0.1), mgp = c(1.5, 0.5, 0))
plot(bs_mat, coef = beta, ylab = "B-Spline Function", mark_knots = "all")
plot(deriv(bs_mat), coef = beta, ylab = "1st Derivative", mark_knots = "all")
plot(deriv(bs_mat, derivs = 2), coef = beta,
     ylab = "2nd Derivative", mark_knots = "all")
plot(update(bs_mat, integral = TRUE), coef = beta,
     ylab = "Integral", mark_knots = "all")
par(op)
