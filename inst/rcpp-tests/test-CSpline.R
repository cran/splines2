Rcpp::sourceCpp("test-CSpline.cpp")

x <- seq.int(0, 10, 0.02)
inter_knots <- c(2.4, 3.5, 5.2, 8)
bound_knots <- c(- 1, 12)
degree <- 4

foo <- function(...) {
    mat <- cSpline(..., intercept = TRUE)
    d1mat <- deriv(mat)
    d2mat <- deriv(d1mat)
    d3mat <- deriv(d2mat)
    list(basis = unclass(mat),
         d1 = unclass(d1mat),
         d2 = unclass(d2mat),
         d3 = unclass(d3mat),
         degree = attr(mat, "degree"),
         internal_knots = knots(mat),
         boundary_knots = attr(mat, "Boundary.knots"))
}

res <- foo(x = x, knots = inter_knots, degree = degree,
           Boundary.knots = bound_knots)

## default constructors with setter methods
res00 <- rcpp_cspline00(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res00)

res01 <- rcpp_cspline01(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res01)

res02 <- rcpp_cspline02(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res02)

res03 <- rcpp_cspline03(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res03)

res04 <- rcpp_cspline04(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res04)

res05 <- rcpp_cspline05(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res05)

## non-default constructor 1
res1 <- rcpp_cspline1(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res1)

## non-default constructor 2
res2 <- rcpp_cspline2(x, 10, degree, bound_knots)
res20 <- foo(x = x, degree = degree, df = 10,
             Boundary.knots = bound_knots)
expect_equivalent(res20, res2)

## non-default constructor 3: simple knot sequence
knot_seq <- sort(c(rep(bound_knots, each = degree + 1), inter_knots))
res31 <- rcpp_cspline3(x, degree, knot_seq)
expect_equivalent(res, res31)

## non-default constructor 3: extended knot sequence
knot_seq <- sort(c(seq.int(0, 10, 1), 1, rep(4, 2), rep(7, 2)))
res32 <- rcpp_cspline3(x, degree, knot_seq)
expect_equivalent(
    lapply(res32, function(a) {
        tmp <- dim(a)
        if (is.null(tmp)) {
            tmp <- length(a)
        }
        tmp
    }),
    {
        tmp <- c(length(x), length(knot_seq) - degree - 1)
        c(rep(list(tmp), 4), 1, length(knot_seq) - 2 * (degree + 1), 2)
    }
)

## non-default constructor 4
res4 <- rcpp_cspline4(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res4)

## conversion from BernsteinPoly
res5 <- rcpp_cspline5(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res5)

## conversion from PeriodicCspline
res6 <- rcpp_cspline6(x, inter_knots, degree, bound_knots)
expect_equivalent(res, res6)
