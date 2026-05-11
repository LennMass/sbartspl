# Helpers: preprocessing for inputs/outputs
# ----------------------------------------------------
# Small utilities for turning a data.frame into a numeric model matrix,
# rank-transforming covariates onto [0, 1], and scaling outcomes to
# [-0.5, 0.5] (and back). Used as a thin preprocessing layer in front of
# BART / SoftBART fits.


# Turn a data.frame into a numeric model matrix (no intercept) and return
# a 0-based group index per column, useful when grouped variable selection
# is needed downstream (e.g. one group per original factor variable).
#
# Args:
#   X : data.frame of covariates (any mix of numeric/factor/character).
#
# Returns a list with:
#   X     : numeric model matrix
#   group : integer vector, one entry per column of X, indexing the
#           original variable each column came from (0-based)
preprocess_df <- function(X) {
  stopifnot(is.data.frame(X))
  
  X <- model.matrix(~.-1, data = X)
  group <- attr(X, "assign") - 1
  
  return(list(X = X, group = group))
  
}


# Tied-rank transform onto [0, 1]. Ties get the maximum rank
# (so the largest value maps to 1, ties are not broken arbitrarily).
trank <- function(x) {
  x_unique <- unique(x)
  x_ranks <- rank(x_unique, ties.method = "max")
  tx <- x_ranks[match(x,x_unique)] - 1
  
  tx <- tx / length(unique(tx))
  tx <- tx / max(tx)
  
  return(tx)
}


# Column-wise rank transform of a matrix or data.frame via trank().
quantile_normalize_bart <- function(X) {
  apply(X = X, MARGIN = 2, trank)
}


# Min-max scale a numeric vector to [-0.5, 0.5]. Standard outcome
# rescaling for BART; pair with unnormalize_bart() to reverse.
normalize_bart <- function(y) {
  a <- min(y)
  b <- max(y)
  y <- (y - a) / (b - a) - 0.5
  return(y)
}

# Inverse of normalize_bart(). Caller must keep the original a = min(y)
# and b = max(y) from the training data and pass them back in here.
#
# Args:
#   z : normalised vector on [-0.5, 0.5]
#   a : original min(y) before normalisation
#   b : original max(y) before normalisation
unnormalize_bart <- function(z, a, b) {
  y <- (b - a) * (z + 0.5) + a
  return(y)
}
