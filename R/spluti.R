AddInterceptPP <- function (PPObj, intercept) {
  if (!inherits(PPObj, "PP")) {
    mstop("'PPObj' does not inherit \"PP\" class!")
  }
  if (!(PPObj$shift)) {
    mstop("\"PP\" or \"PPA\" object in unshifted form is not supported!")
  }
  PPCoef <- PPObj$coef
  if (missing(intercept)) {
    mstop("You need to provide intercept values!")
  }
  intercept <- c(intercept)
  if (is.matrix(PPCoef)) {
    if (length(intercept) == 1) {
      PPCoef[1, ] <- PPCoef[1, ] + intercept
    } else {
      PPCoef <- array(PPCoef, c(dim(PPCoef), length(intercept)))
      PPCoef[1, , ] <- PPCoef[1, , ] + rep(intercept, each = dim(PPCoef)[2])
    }
  } else {
    if (length(intercept) == 1) {
      PPCoef[1, , ] <- PPCoef[1, , ] + intercept
    } else if (length(intercept) == dim(PPCoef)[3]) {
      PPCoef[1, , ] <- PPCoef[1, , ] + rep(intercept, each = dim(PPCoef)[2])
    } else {
      mstop("Incorrect number of intercept values!",
            "Either set one common intercept for all splines or one per each.")
    }
  }
  lst <- list(coef = PPCoef, knots = PPObj$knots, shift = PPObj$shift)
  if (is.matrix(PPCoef)) structure(lst, class = "PP")
  else structure(lst, class = c("PP", "PPA"))
}

AddPP <- function (a, b) {
  if (!(inherits(a, "PP") && inherits(b, "PP"))) {
    mstop("'a' and 'b' must both inherit \"PP\" class!")
  }
  if (!(a$shift && b$shift)) {
    mstop("\"PP\" or \"PPA\" object in unshifted form is not supported!")
  }
  if (!identical(a$knots, b$knots)) {
    new.knots <- sort(union(a$knots, b$knots))
    a <- MapPP(a, new.knots)
    b <- MapPP(b, new.knots)
  }
  if (inherits(a, "PPA")) {
    if (inherits(b, "PPA")) PPA_add_PPA(a, b)
    else PP_add_PPA(b, a)
  } else {
    if (inherits(b, "PPA")) PP_add_PPA(a, b)
    else PP_add_PP(a, b)
  }
}

PP_add_PP <- function (a, b) {
  PPCoef1 <- a$coef  ## matrix
  PPCoef2 <- b$coef  ## matrix
  ord1 <- dim(PPCoef1)[1]
  ord2 <- dim(PPCoef2)[1]
  PPObj <- b
  if (ord1 == ord2) {
    PPCoef <- PPCoef1 + PPCoef2
  } else if (ord1 < ord2) {
    PPCoef <- PPCoef2
    PPCoef[1:ord1, ] <- PPCoef[1:ord1, ] + PPCoef1
  } else {
    PPCoef <- PPCoef1
    PPCoef[1:ord2, ] <- PPCoef[1:ord2, ] + PPCoef2
  }
  PPObj$coef <- PPCoef
  PPObj
}

PP_add_PPA <- function (a, b) {
  PPCoef1 <- a$coef  ## matrix
  PPCoef2 <- b$coef  ## array
  ord1 <- dim(PPCoef1)[1]
  ord2 <- dim(PPCoef2)[1]
  PPObj <- b
  if (ord1 == ord2) {
    PPCoef <- c(PPCoef1) + PPCoef2
  } else if (ord1 < ord2) {
    PPCoef <- PPCoef2
    PPCoef[1:ord1, , ] <- PPCoef[1:ord1, , ] + c(PPCoef1)
  } else {
    PPCoef <- array(PPCoef1, dim = c(ord1, dim(PPCoef2)[2:3]))
    PPCoef[1:ord2, , ] <- PPCoef[1:ord2, , ] + PPCoef2
  }
  PPObj$coef <- PPCoef
  PPObj
}

PPA_add_PPA <- function (a, b) {
  PPCoef1 <- a$coef  ## array
  PPCoef2 <- b$coef  ## array
  if (dim(PPCoef1)[3] != dim(PPCoef2)[3]) {
    mstop("'a' and 'b' do not have the same number of splines!")
  }
  ord1 <- dim(PPCoef1)[1]
  ord2 <- dim(PPCoef2)[1]
  PPObj <- b
  if (ord1 == ord2) {
    PPCoef <- PPCoef1 + PPCoef2
  } else if (ord1 < ord2) {
    PPCoef <- PPCoef2
    PPCoef[1:ord1, , ] <- PPCoef[1:ord1, , ] + PPCoef1
  } else {
    PPCoef <- PPCoef1
    PPCoef[1:ord2, , ] <- PPCoef[1:ord2, , ] + PPCoef2
  }
  PPObj$coef <- PPCoef
  PPObj
}

BSpl2PP <- function (x, fun = "bs", df = NULL, degree = 3, intercept = FALSE,
                     knots = NULL, Boundary.knots = range(x), coef) {
  if (missing(coef)) {
    mstop("Basis coefficients must be provided via argument 'coef'!")
  }
  BSplCall <- MakeBSplCall(x, fun, df, degree, intercept, knots, Boundary.knots)
  PPCoef <- BSpl2PP_Kernel(BSplCall, coef)
  lst <- list(coef = PPCoef, knots = BSplKnots(BSplCall), shift = TRUE)
  if (is.matrix(PPCoef)) structure(lst, class = "PP")
  else structure(lst, class = c("PPA", "PP"))
}


MakeBSplCall <- function (x, fun = "bs",
                          df = NULL, degree = 3, intercept = FALSE,
                          knots = NULL, Boundary.knots = range(x)) {
  fun_set <- c("bs", "ns", "pbs")
  if (!(fun %in% fun_set)) {
    mstop(sprintf("'fun' can only be: %s", toString(sprintf("'%s'", fun_set))))
  }
  if (fun == "ns") {
    degree <- 3
    min_k <- 0
    min_df <- 1 + intercept
  } else if (fun == "pbs") {
    min_k <- degree
    min_df <- intercept
  } else {
    min_k <- 0
    min_df <- degree + intercept
  }
  if (any(is.infinite(x))) mstop("'x' should not contain Inf or -Inf!")
  if (anyNA(x)) x <- x[!is.na(x)]
  min_x <- min(x)
  max_x <- max(x)
  if (min_x == max_x) mstop("'x' should contain more unique values!")
  Boundary.knots <- sort(Boundary.knots)
  lbnd <- Boundary.knots[1]
  rbnd <- Boundary.knots[2]
  if (!is.finite(lbnd) || !is.finite(rbnd) || lbnd == rbnd) {
    mstop("The 2 boundary knots should be non NA/NaN/Inf/-Inf and different!")
  }
  if (any(x < lbnd | x > rbnd)) mstop("x-value outside boundary is not allowed here!")
  if (is.null(knots)) {
    if (is.null(df)) {
      mstop("'df' and 'knots' can't both be NULL!")
    } else {
      k <- df - min_df
      if (k < min_k) mstop(sprintf("'df' should be %d at minimum!", min_df + min_k))
      prob <- seq.int(0, 1, length.out = k + 2)
      prob <- prob[2:(k + 1)]
      knots <- quantile(x, prob, names = FALSE)
    }
  } else {
    k <- length(knots)
    if (k < min_k) mstop(sprintf("At least %d interior knots are required!", min_k))
    knots <- sort(knots)
    if (anyDuplicated(knots)) mstop("Interior knots should not contain duplicates!")
    if ((knots[1] <= min_x) || (knots[k] >= max_x)) {
      mstop("Interior knots must be within range of data 'x'!")
    }
    df_implied <- k + min_df
    if (is.null(df)) {
      df <- df_implied
    } else {
      if (df != df_implied) {
        mstop("Inconsistency between specified 'df' and 'knots'",
              "You may consider setting either to NULL!")
      }
    }
  }
  if (fun == "ns") {
    call("ns", quote(x), df = df, knots = knots, intercept = intercept,
         Boundary.knots = Boundary.knots)
  } else {
    call(fun, quote(x), df = df, knots = knots, degree = degree, intercept = intercept,
         Boundary.knots = Boundary.knots)
  }
}

BSplKnots <- function (BSplCall) {
  lst <- as.list(BSplCall)
  interior_knots <- lst$knots
  lbnd_knots <- lst$Boundary.knots[1]
  rbnd_knots <- lst$Boundary.knots[2]
  c(lbnd_knots, interior_knots, rbnd_knots)
}

BSpl2PP_Kernel <- function (BSplCall, BSplCoef) {
  if (is.matrix(BSplCoef)) n_replicates <- ncol(BSplCoef)
  else n_replicates <- 1
  x <- BSplKnots(BSplCall)
  degree <- BSplCall$degree
  if (is.null(degree)) degree <- 3  ## for ns()
  ord <- degree + 1
  df <- BSplCall$df
  if ((n_replicates > 1) && (nrow(BSplCoef) != df)) {
      mstop(sprintf("Number of coefficients (rows) mismatch number of basis: %d!", df))
    }
  if ((n_replicates == 1) && (length(BSplCoef) != df)) {
      mstop(sprintf("Number of coefficients mismatch number of basis: %d!", df))
    }
  n_pieces <- length(x) - 1
  xg <- .Call("C_EvalGrid", x, ord)
  BSplBasis <- eval(BSplCall, list(x = xg))
  yg <- BSplBasis %*% BSplCoef  ## (ord x n_pieces) x n_replicates
  yg <- matrix(yg, ord, n_pieces * n_replicates)
  P <- outer(seq.int(0, 1, length.out = ord), 0:degree, "^")
  PPCoef <- solve.default(P, yg)
  h <- diff.default(x)
  H <- rep(1 / h, each = ord) ^ rep.int(0:degree, length(n_pieces))
  PPCoef <- H * PPCoef
  if (n_replicates > 1) PPCoef <- array(PPCoef, dim = c(ord, n_pieces, n_replicates))
  PPCoef
}

DerivPoly <- function (pc, deriv = 0) {
  degree <- length(pc) - 1
  if (deriv > degree) return(0)
  if (deriv > 0) {
    pc[-seq_len(deriv)] * choose(deriv:degree, deriv) * factorial(deriv)
  } else {
    pc
  }
}

DerivPP <- function (PPObj, deriv = 0) {
  if (!inherits(PPObj, "PP")) {
    mstop("'PPObj' does not inherit \"PP\" class!")
  }
  if (!(PPObj$shift)) {
    mstop("\"PP\" or \"PPA\" object in unshifted form is not supported!")
  }
  PPCoef <- PPObj$coef
  degree <- dim(PPCoef)[1] - 1
  if (deriv > degree) {
    mstop(sprintf("Piecewise polynomial in this object has degree %d.", degree),
          sprintf("deriv <= %d expected, or polynomial becomes 0 after differentiation!", degree))
  }
  if (deriv > 0) {
    MulFctr <- choose(deriv:degree, deriv) * factorial(deriv)
    if (is.matrix(PPCoef)) {
      PPCoef <- PPCoef[-seq_len(deriv), , drop = FALSE] * MulFctr
    } else {
      PPCoef <- PPCoef[-seq_len(deriv), , , drop = FALSE] * MulFctr
    }
    PPObj$coef <- PPCoef
  }
  PPObj
}

General2PP <- function (degree, knots, PredFun, ...) {
  
}


ICSpl2PP <- function (x, y, method = "natural") {
  method_set <- c("fmm", "natural", "periodic", "monoH.FC", "hyman")
  if (!(method %in% method_set)) {
    mstop(sprintf("'method' can only be: %s!", toString(sprintf("'%s'", method_set))))
  }
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]
    y <- y[ind]
  }
  if (anyDuplicated(x)) mstop("Tied 'x' values are not allowed!")
  if ((method == "periodic") && (y[1] != y[length(y)])) 
    mstop("y-values are not periodic!")
  if (method %in% c("monoH.FC", "hyman")) {
    dy <- diff.default(y)
    if (!(all(dy <= 0) || all(dy >= 0))) mstop("y-values are not monotonic!")
  }
  if (TRUE) ICSpl2PP_v1(x, y, method)
  else ICSpl2PP_v2(x, y, method)
}

ICSpl2PP_v1 <- function (x, y, method = "natural") {
  if (method == "monoH.FC") mstop("'monoH.FC' method is not supported!")
  ICSpl <- stats::splinefun(x, y, method)
  info <- environment(ICSpl)$z
  ind <- seq_len(length(info$x) - 1)
  PPCoef <- with(info, rbind(y[ind], b[ind], c[ind], d[ind], deparse.level = 0))
  structure(list(coef = PPCoef, knots = info$x, shift = TRUE), class = "PP")
}

ICSpl2PP_v2 <- function (x, y, method = "natural") {
  ICSpl <- stats::splinefun(x, y, method)
  xl <- x[-length(x)]
  y0 <- y[-length(y)]
  y1 <- ICSpl(xl, deriv = 1)
  y2 <- ICSpl(xl, deriv = 2)
  xr <- x[-1]
  xm <- (xl + xr) / 2
  y3 <- ICSpl(xm, deriv = 3)
  a <- y0
  b <- y1
  c <- y2 / 2
  d <- y3 / 6
  PPCoef <- rbind(a, b, c, d, deparse.level = 0)
  structure(list(coef = PPCoef, knots = x, shift = TRUE), class = "PP")
}



MapPP <- function (PPObj, new.knots) {
  if (!inherits(PPObj, "PP")) {
    mstop("'PPObj' does not inherit \"PP\" class!")
  }
  if (!(PPObj$shift)) {
    mstop("\"PP\" or \"PPA\" object in unshifted form is not supported!")
  }
  knots <- PPObj$knots
  if (identical(knots, new.knots)) return(PPObj)
  if (!all(knots %in% new.knots)) {
    mstop("'new.knots' is not more refined than 'PPObj$knots'!",
          "'PPObj$knots' needs be inclusive in 'new.knots' for this function to work!")
  }
  PPCoef <- PPObj$coef
  aug_ind <- findInterval(new.knots[-length(new.knots)], knots, rightmost.closed = TRUE)
  aug_knots <- knots[aug_ind]
  if (is.matrix(PPCoef)) aug_PPCoef <- PPCoef[, aug_ind]
  else aug_PPCoef <- PPCoef[, aug_ind, ]
  degree <- dim(PPCoef)[1] - 1
  power <- 0:degree
  BinomCoefMat <- outer(power, power, choose)
  n_pieces <- length(new.knots) - 1
  if (is.matrix(PPCoef)) {
    for (i in 1:n_pieces) {
      a <- aug_knots[i]
      b <- new.knots[i]
      if (b != a) {
        PowerMat <- toeplitz((b - a) ^ power)
        aug_PPCoef[, i] <- base::crossprod(BinomCoefMat * PowerMat, aug_PPCoef[, i])
      }
    }
  } else {
    for (i in 1:n_pieces) {
      a <- aug_knots[i]
      b <- new.knots[i]
      if (b != a) {
        PowerMat <- toeplitz((b - a) ^ power)
        aug_PPCoef[, i, ] <- base::crossprod(BinomCoefMat * PowerMat, aug_PPCoef[, i, ])
      }
    }
  }
  lst <- list(coef = aug_PPCoef, knots = new.knots, shift = TRUE)
  if (is.matrix(PPCoef)) structure(lst, class = "PP")
  else structure(lst, class = c("PPA", "PP"))
}

mstop <- function (...) {
  stop(sprintf("\n  %s", unlist(list(...))), call. = FALSE)
}

ListTerms <- function (termsObject) {
  tl <- attr(termsObject, "term.labels")  ## tl <- labels(termsObject)
  cat("Available terms are:\n", sprintf(" * %s\n", tl))
}

ShiftPoly <- function (pc, a, b) {
  degree <- length(pc) - 1
  power <- 0:degree
  BinomCoefMat <- outer(power, power, choose)
  x0 <- b - a
  PowerMat <- toeplitz(x0 ^ power)
  c(base::crossprod(BinomCoefMat * PowerMat, pc))
}

Poly2PP <- function (pc, knots) {
  knots <- sort(knots)
  if (anyDuplicated(knots)) mstop("Knots should not contain duplicates!")
  n_pieces <- length(knots) - 1
  if (n_pieces < 1) mstop("At least 2 knots are expected!")
  if (is.matrix(pc)) degree <- nrow(pc) - 1
  else degree <- length(pc) - 1
  if (degree < 1) mstop("The polynomial should be at least linear!")
  power <- 0:degree
  BinomCoefMat <- outer(power, power, choose)
  if (is.matrix(pc)) {
    PPCoef <- array(0, dim = c(degree + 1, n_pieces, ncol(pc)))
    for (i in 1:n_pieces) {
      PowerMat <- toeplitz(knots[i] ^ power)
      PPCoef[, i, ] <- base::crossprod(BinomCoefMat * PowerMat, pc)
    }
  } else {
    PPCoef <- matrix(0, degree + 1, n_pieces)
    for (i in 1:n_pieces) {
      PowerMat <- toeplitz(knots[i] ^ power)
      PPCoef[, i] <- base::crossprod(BinomCoefMat * PowerMat, pc)
    }
  }
  lst <- list(coef = PPCoef, knots = knots, shift = TRUE)
  if (is.matrix(PPCoef)) structure(lst, class = "PP")
  else structure(lst, class = c("PPA", "PP"))
}

PolyEqn <- function (pc, shift, a) {
  
  if (length(pc) == 1) return(sprintf("%.3g", pc))
  
  pc0 <- pc[1]       ## coefficient for 1
  pc1 <- pc[2]       ## coefficient for x
  pc2 <- pc[-c(1, 2)]  ## coefficients for higher power term; may be numeric(0)
  
  
  rev_a_sgn <- if (a > 0) " - " else " + "
  a <- abs(a)
  
  xterm0 <- sprintf("%.3g", pc0)
  coef_sgn <- if (pc1 > 0) " + " else " - "
  if (shift) {
    xterm1 <- sprintf("%s%.3g * (x%s%.3g)", coef_sgn, abs(pc1), rev_a_sgn, a)
  } else {
    xterm1 <- sprintf("%s%.3g * x", coef_sgn, abs(pc1))
  }
  
  if (length(pc2)) {
    sgn <- rep.int(" - ", length(pc2))
    sgn[pc2 > 0] <- " + "
    if (shift) {
      xterm2 <- sprintf("%s%.3g * (x%s%.3g) ^ %d", sgn, abs(pc2), rev_a_sgn, a, seq_along(pc2) + 1)
    } else {
      xterm2 <- sprintf("%s%.3g * x ^ %d", sgn, abs(pc2), seq_along(pc2) + 1)
    }
    xterm2 <- paste0(xterm2, collapse = "")
  } else xterm2 <- ""
  
  paste0(xterm0, xterm1, xterm2, collapse = "")
}

polydiv <- function (f, g) {
  if (length(f) < length(g)) {
    mstop("f's order is smaller than g's")
  }
  out <- .Call("C_polydiv", f, g)
  list(quotient = out[length(g):length(f)], remainder = out[seq_len(length(g) - 1)])
}

multiplicity <- function (f, g, tol) {
  if (!is.matrix(g)) g <- matrix(g)
  if (length(f) < nrow(g)) {
    mstop("f's order is smaller than g's")
  }
  .Call("C_multiplicity", f, g, tol)
}

polysturm <- function (pc) {
  n <- length(pc)
  sturm <- vector("list", n)
  f <- pc
  g <- pc[-1] * seq_len(n - 1)  ## coefficients of differentiated polynomial
  f <- f / abs(f[n])
  g <- g / abs(g[n - 1])
  sturm[[1]] <- f
  sturm[[2]] <- g
  i <- 3
  while (n > 2) {
    r <- -1 * polydiv(f, g)$remainder
    if (max(abs(r)) < 1e-8) break
    r <- r / abs(r[n - 2])
    sturm[[i]] <- r
    i <- i + 1
    n <- n - 1
    f <- g
    g <- r
  }
  sturm[lengths(sturm) > 0]
}

polyzero <- function (f) {
  sturm <- polysturm(f)
  sgn2 <- unlist(lapply(sturm, tail, 1))
  sgn1 <- sgn2 * (-1) ^ seq(from = length(f) - 1, by = -1, length.out = length(sturm))
  v1 <- sum(diff.default(sgn1) != 0)
  v2 <- sum(diff.default(sgn2) != 0)
  n_distinct_real <- v1 - v2
  GCD <- sturm[[length(sturm)]]
  f_deflated <- f
  if (length(GCD) > 1) f_deflated <- polydiv(f, GCD)$quotient
  Roots <- polyroot(f_deflated)
  RePart <- Re(Roots)
  ImPart <- Im(Roots)
  Im2Re <- abs(ImPart) / abs(RePart)
  tol <- 2 ^ (-26)  ## sqrt(.Machine$double.eps)
  isReal <- (Im2Re < tol)
  while (sum(isReal) < n_distinct_real) {
    tol <- tol / 2
    isReal <- (Im2Re < tol)
  }
  while (sum(isReal) > n_distinct_real) {
    tol <- tol * 2
    isReal <- (Im2Re < tol)
  }
  isComplex <- !isReal
  RealRoots <- RePart[isReal]
  ComplexRoots <- Roots[isComplex]
  ComplexRoots <- ComplexRoots[Im(ComplexRoots) > 0]
  n_distinct_complex <- length(ComplexRoots)  ## don't divide 2
  if (n_distinct_real > 1)
    RealRoots <- RealRoots[order(abs(RealRoots), decreasing = TRUE)]
  if (n_distinct_complex > 1)
    ComplexRoots <- ComplexRoots[order(Mod(ComplexRoots), decreasing = TRUE)]
  RealRootsMultiplicity <- integer(n_distinct_real)
  if (n_distinct_real) {
    g <- rbind(-RealRoots, 1)
    RealRootsMultiplicity <- multiplicity(f, g, tol)
  }
  ComplexRootsMultiplicity <- integer(n_distinct_complex)
  if (n_distinct_complex) {
    g <- rbind(Mod(ComplexRoots) ^ 2, -2 * Re(ComplexRoots), 1)
    ComplexRootsMultiplicity <- multiplicity(f, g, tol)
  }
  list(real = RealRoots, mreal = RealRootsMultiplicity,
       complex = ComplexRoots, mcomplex = ComplexRootsMultiplicity)
}

PPAEntry <- function (PPAObj, i) {
  if (!inherits(PPAObj, "PPA")) {
    mstop("'PPAObj' is not a 'PPA' object!")
  }
  if (length(i) > 1) mstop("Only a single value for 'i' is allowed!")
  if ((i < 1) || (i > dim(PPAObj$coef)[3])) mstop("Out-of-bound 'i' value!")
  PPCoef <- PPAObj$coef[, , i]
  if (!is.matrix(PPCoef)) PPCoef <- matrix(PPCoef, nrow = 1)
  structure(list(coef = PPCoef, knots = PPAObj$knots, shift = PPAObj$shift),
            class = "PP")
}


summary.PPA <- function (object, nhead = NULL, ...) {
  n <- dim(object$coef)[3]
  EqnLst <- vector("list", n)
  for (i in 1:n) {
    PPCoef <- object$coef[, , i]
    if (!is.matrix(PPCoef)) PPCoef <- matrix(PPCoef, nrow = 1)
    PPObj <- with(object, list(coef = PPCoef, knots = knots, shift = shift))
    EqnLst[[i]] <- summary.PP(PPObj, nhead)  ## skip S3 method dispatching
  }
  EqnLst
}

print.PPA <- function (x, ...) {
  cat("----------\n")
  cat(str(x))
  cat("----------\n")
  cat(sprintf("An array of %d splines, each with %d polynomial pieces of degree %d.\n",
              dim(x$coef)[3], dim(x$coef)[2], dim(x$coef)[1] - 1))
  cat("No preview of polynomial equations is done here.\n")
  cat("But you can use 'summary' to export all equation strings.\n")
  cat("----------\n")
}

solve.PPA <- function (a, b = 0, deriv = 0, ...) {
  if (!(a$shift)) {
    mstop("This function does not support \"PP\" objects in unshifted form!")
  }
  n <- dim(a$coef)[3]
  xr_lst <- vector("list", n)
  for (i in 1:n) {
    PPCoef <- a$coef[, , i]
    if (!is.matrix(PPCoef)) PPCoef <- matrix(PPCoef, nrow = 1)
    PPObj <- with(a, list(coef = PPCoef, knots = knots, shift = shift))
    xr_lst[[i]] <- solve.PP(PPObj, b, deriv)  ## skip S3 method dispatching
  }
  xr_lst
}



summary.PP <- function (object, nhead = NULL, ...) {
  PPObj <- object
  PPCoef <- PPObj$coef
  shift <- PPObj$shift
  knots <- PPObj$knots
  n_pieces <- dim(PPCoef)[2]
  if (is.null(nhead)) nhead <- n_pieces
  else nhead <- min(nhead, n_pieces)
  EqnVec <- character(nhead)
  for (i in 1:nhead) EqnVec[i] <- PolyEqn(PPCoef[, i], shift, knots[i])
  EqnVec
}

print.PP <- function (x, ...) {
  cat("----------\n")
  cat(str(x))
  preview <- summary.PP(x, 6)
  cat("----------\n")
  cat(sprintf("%d polynomial pieces of degree %d in this spline.\n",
              length(x$knots) - 1, dim(x$coef)[1] - 1))
  cat(sprintf("The first %d are printed below for preview.\n", length(preview)))
  cat("You can use 'summary' to export all equations strings.\n")
  cat("----------\n")
  cat("", paste0(preview, "\n"))
  cat("----------\n")
}

plot.PP <- function (x, spread = 3, deriv = 0, show.knots = FALSE, ylim.hint = NULL, ...) {
  PPObj <- x
  if (!(PPObj$shift)) {
    mstop("This function does not support \"PP\" objects in unshifted form!")
  }
  PPObj <- DerivPP(PPObj, deriv)
  PPCoef <- PPObj$coef
  n_pieces <- dim(PPCoef)[2]
  degree <- dim(PPCoef)[1] - 1
  x <- PPObj$knots
  power <- 0:degree
  ng <- spread * (degree + 1)
  xg <- numeric(ng * n_pieces)
  if (is.matrix(PPCoef)) {
    yg <- numeric(length(xg))
    ind <- 1:ng
    for (i in 1:n_pieces) {
      xg[ind] <- xg_i <- seq.int(x[i], x[i + 1], length.out = ng)
      yg[ind] <- outer(xg_i - x[i], power, "^") %*% PPCoef[, i]
      ind <- ind + ng
    }
  } else {
    yg <- matrix(0, length(xg), dim(PPCoef)[3])
    ind <- 1:ng
    for (i in 1:n_pieces) {
      xg[ind] <- xg_i <- seq.int(x[i], x[i + 1], length.out = ng)
      yg[ind, ] <- outer(xg_i - x[i], power, "^") %*% PPCoef[, i, ]
      ind <- ind + ng
    }
  }
  ylim <- range(yg, ylim.hint)
  if (is.matrix(yg)) {
    matplot(xg, yg, type = "l", lty = 1, xlab = "x", ylab = "y", ylim = ylim)
  } else {
    plot.default(xg, yg, type = "l", xlab = "x", ylab = "y", ylim = ylim)
  }
  if (show.knots) abline(v = x, lty = 3, col = 8)
  xgyg <- list(x = xg, y = yg)
}

predict.PP <- function (object, newx, deriv = 0, ...) {
  PPObj <- object
  if (!(object$shift)) {
    mstop("This function does not support \"PP\" objects in unshifted form!")
  }
  if (missing(newx)) mstop("No 'newx' provided!")
  PPCoef <- PPObj$coef
  degree <- dim(PPCoef)[1] - 1
  if (deriv > degree) {
    mstop(sprintf("Piecewise polynomial in this object has degree %d.", degree),
          sprintf("deriv <= %d expected, or polynomial becomes 0 after differentiation!", degree))
  }
  power <- 0:(degree - deriv)
  x <- PPObj$knots
  piece <- findInterval(newx, x, rightmost.closed = TRUE)
  if (any(piece == 0 | piece == length(x))) {
    mstop("Out-of-boundary prediction is not allowed!")
  }
  ind <- split.default(seq_len(length(newx)), piece)
  piece_with_newx <- as.integer(names(ind))
  n_pieces_to_predict <- length(piece_with_newx)
  if (is.matrix(PPCoef)) {
    y <- numeric(length(newx))
    MulFctr <- choose(deriv:degree, deriv) * factorial(deriv)
    for (i in 1:n_pieces_to_predict) {
      ii <- piece_with_newx[i]
      xi <- newx[ind[[i]]] - x[ii]
      pc <- PPCoef[, ii]
      if (deriv > 0) {
        pc <- pc[-seq_len(deriv)] * MulFctr
      }
      y[ind[[i]]] <- outer(xi, power, "^") %*% pc
    }
  } else {
    y <- matrix(0, length(newx), dim(PPCoef)[3])
    for (i in 1:n_pieces_to_predict) {
      ii <- piece_with_newx[i]
      xi <- newx[ind[[i]]] - x[ii]
      pc <- PPCoef[, ii, ]
      if (deriv > 0) {
        pc <- pc[-seq_len(deriv), ] * MulFctr
      }
      y[ind[[i]], ] <- outer(xi, power, "^") %*% pc
    }
  }
  y
}

solve.PP <- function (a, b = 0, deriv = 0, ...) {
  PPObj <- a
  y <- b
  if (!(PPObj$shift)) {
    mstop("This function does not support \"PP\" objects in unshifted form!")
  }
  PPCoef <- PPObj$coef
  n_pieces <- dim(PPCoef)[2]
  degree <- dim(PPCoef)[1] - 1
  x <- PPObj$knots
  if (deriv >= degree) {
    mstop(sprintf("Piecewise polynomial in this object has degree %d.", degree),
          sprintf("deriv < %d is expected in the inverse problem!", degree))
  }
  xr <- vector("list", n_pieces)
  MulFctr <- choose(deriv:degree, deriv) * factorial(deriv)
  for (i in 1:n_pieces) {
    pc <- PPCoef[, i]
    if (deriv > 0) {
      pc <- pc[-seq_len(deriv)] * MulFctr
    }
    pc[1] <- pc[1] - y
    if (TRUE) { 
      croots <- base::polyroot(pc)
      rroots <- Re(croots)[round(Im(croots), 10) == 0]
      rroots <- base::unique(rroots)
    } else {
      zeros <- polyzero(pc)
      rroots <- zeros$real
    }
    rroots <- rroots + x[i]
    xr[[i]] <- rroots[(rroots >= x[i]) & (rroots <= x[i + 1])]
    }
  unlist(xr)
}


ExtractBSplTerm <- function (termsObject, BSplTerm) {
  tl <- attr(termsObject, "term.labels")  ## tl <- labels(termsObject)
  if (!(BSplTerm %in% tl)) {
    h <- "\n  Required BSplTerm not found! Available terms are:\n"
    ol <- sprintf("    * %s\n", tl)
    stop(h, ol, call. = FALSE)
  }
  is_bs <- grepl("^bs\\(", BSplTerm)
  is_ns <- grepl("^ns\\(", BSplTerm)
  if (!(is_bs || is_ns)) {
    mstop("Provided BSplTerm is neither 'bs()' nor 'ns()' from package 'splines'!")
  }
  pos <- match(BSplTerm, tl)
  predvars <- attr(termsObject, "predvars")  ## try: as.list(predvars)
  BSplCall <- predvars[[2 + pos]]  ## try: as.list(BSplCall)
  BSplCall[[2]] <- quote(x)
  if (is_bs) {
    degree <- BSplCall[[3]]
    interior_knots <- unname(BSplCall[[4]])
    boundary_knots <- BSplCall[[5]]
    BSplCall[[4]] <- interior_knots
  } else {
    degree <- 3
    interior_knots <- unname(BSplCall[[3]])
    boundary_knots <- BSplCall[[4]]
    BSplCall[[3]] <- interior_knots
  }
  x <- c(boundary_knots[1], interior_knots, boundary_knots[2])
  df <- length(interior_knots) + degree + 1 - (1 - BSplCall$intercept) - 2 * is_ns
  BSplCall$df <- df
  list(pos = pos, knots = x, call = BSplCall)
}


RegBSpl2PP <- function (RegModel, BSplTerm) {
  UseMethod("RegBSpl2PP")
}

RegBSpl2PP.lm <- function (RegModel, BSplTerm) {
  RegBSpl <- ExtractBSplTerm(terms(RegModel), BSplTerm)
  pos <- RegBSpl$pos
  RegBSplCoef <- with(RegModel, coefficients[assign == pos])
  PPCoef <- BSpl2PP_Kernel(RegBSpl$call, RegBSplCoef)
  lst <- list(coef = PPCoef, knots = RegBSpl$knots, shift = TRUE)
  structure(lst, class = "PP")
}

RegBSpl2PP.mlm <- function (RegModel, BSplTerm) {
  RegBSpl <- ExtractBSplTerm(terms(RegModel), BSplTerm)
  pos <- RegBSpl$pos
  RegBSplCoef <- with(RegModel, coefficients[assign == pos, , drop = FALSE])
  PPCoef <- BSpl2PP_Kernel(RegBSpl$call, RegBSplCoef)
  lst <- list(coef = PPCoef, knots = RegBSpl$knots, shift = TRUE)
  structure(lst, class = c("PPA", "PP"))
}

RegBSpl2PP.lme <- function (RegModel, BSplTerm) {
  RegBSpl <- ExtractBSplTerm(terms(RegModel), BSplTerm)
  ind <- attr(RegModel$fixDF, "assign")[[BSplTerm]]
  RegBSplCoef <- RegModel$coefficients$fixed[ind]
  PPCoef <- BSpl2PP_Kernel(RegBSpl$call, RegBSplCoef)
  lst <- list(coef = PPCoef, knots = RegBSpl$knots, shift = TRUE)
  structure(lst, class = "PP")
}



SmSpl2PP <- function (SmSpl) {
  if (!inherits(SmSpl, "smooth.spline")) {
    mstop("'SmSpl' is not a 'smooth.spline' object!")
  }
  kx <- with(SmSpl$fit, knot * range + min)
  kx <- kx[4:(length(kx) - 3)]
  ky <- predict(SmSpl, kx)[[2]]  ## see ?predict.smooth.spline
  ICSpl2PP(kx, ky, method = "natural")
}

GenToyReg <- function () {
  set.seed(0)
  n_data <- 101
  n_basis <- 6
  n_replicates <- 4
  x <- seq(0, 10, length = n_data)
  B <- bs(x, df = n_basis, intercept = TRUE)
  b <- matrix(rnorm(n_basis * n_replicates), n_basis) 
  y_true <- B %*% b
  y_true[, 2] <- y_true[, 3] + 0.2 * x - 0.5
  y_true[, 3] <- y_true[, 3] + 0.1 * x - 1
  y_true[, 4] <- y_true[, 4] - 0.25 * x
  matplot(x, y_true, type = "l", lty = 1)  ## looks good!
  e <- matrix(rnorm(n_data * n_replicates, 0, 0.2 * sd(y_true)), nrow = n_data)
  y <- round(y_true + e, 2)
  matplot(x, y, type = "p", lty = 1, pch = 20)
  wtoyreg <- data.frame(x = x, y1 = y[, 1], y2 = y[, 2], y3 = y[, 3], y4 = y[, 4])
  save(wtoyreg, file = "wtoyreg.rda", compress = "xz")
  ind1 <- sort(sample.int(100, size = 50)); x1 <- x[ind1]; y1 <-y[ind1, 1]
  ind2 <- sort(sample.int(100, size = 50)); x2 <- x[ind2]; y2 <-y[ind2, 2]
  ind3 <- sort(sample.int(100, size = 50)); x3 <- x[ind3]; y3 <-y[ind3, 3]
  ind4 <- sort(sample.int(100, size = 50)); x4 <- x[ind4]; y4 <-y[ind4, 4]
  ltoyreg <- data.frame(x = c(x1, x2, x3, x4), y = c(y1, y2, y3, y4), id = gl(4, 50))
  save(ltoyreg, file = "ltoyreg.rda", compress = "xz")
  with(ltoyreg, plot(x, y, pch = 20, col = as.integer(id)))
}

GenToyInterp <- function () {
  set.seed(0)
  x <- 1:10 + runif(10, -0.1, 0.1)
  y1 <- rnorm(10, 3, 1)
  y2 <- y1[c(1:9, 1)]
  y3 <- sort(y1)
  toyinterp <- data.frame(x, y1, y2, y3)
  save(toyinterp, file = "toyinterp.rda", compress = "xz")
}

UnshiftPP <- function (PPObj) {
  if (!inherits(PPObj, "PP")) {
    mstop("'PPObj' does not inherit \"PP\" class!")
  }
  if (PPObj$shift) {
    PPCoef <- PPObj$coef
    x <- PPObj$knots
    degree <- dim(PPCoef)[1] - 1
    power <- 0:degree
    n_pieces <- length(x) - 1
    NewPPCoef <- PPCoef
    BinomCoefMat <- outer(power, power, choose)
    if (is.matrix(PPCoef)) {
      for (i in 1:n_pieces) {
        PowerMat <- toeplitz((-x[i]) ^ power)
        NewPPCoef[, i] <- base::crossprod(BinomCoefMat * PowerMat, PPCoef[, i])
      }
    } else {
      for (i in 1:n_pieces) {
        PowerMat <- toeplitz((-x[i]) ^ power)
        NewPPCoef[, i, ] <- base::crossprod(BinomCoefMat * PowerMat, PPCoef[, i, ])
      }
    }
    PPObj$coef <- NewPPCoef
    PPObj$shift <- FALSE
  }
  warning("Use of unshifted form is deprecated in {spluti}.\n", call. = FALSE)
  PPObj
}

