
preprocess_df <- function(X) {
  stopifnot(is.data.frame(X))
  
  X <- model.matrix(~.-1, data = X)
  group <- attr(X, "assign") - 1
  
  return(list(X = X, group = group))
  
}

trank <- function(x) {
  x_unique <- unique(x)
  x_ranks <- rank(x_unique, ties.method = "max")
  tx <- x_ranks[match(x,x_unique)] - 1
  
  tx <- tx / length(unique(tx))
  tx <- tx / max(tx)
  
  return(tx)
}

quantile_normalize_bart <- function(X) {
  apply(X = X, MARGIN = 2, trank)
}

normalize_bart <- function(y) {
  a <- min(y)
  b <- max(y)
  y <- (y - a) / (b - a) - 0.5
  return(y)
}

unnormalize_bart <- function(z, a, b) {
  y <- (b - a) * (z + 0.5) + a
  return(y)
}
