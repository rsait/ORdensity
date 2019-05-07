density <- function(dx, label,  K=10)
{
  # densidad de Falsos Positivos en un area de K vecinos
  # Input: 
  # dx0: distances from x0 to the rest of individual
  # label: label indicating whether suspuciuos (1) or null (0)
  # Output: densidad
  o <- order(dx)[1:K]
  r <- (dx[o])[K]
  # p <- sum(label[o]==0)/K
  p <- sum(label[o]==0)
  f <- p/r
  out <- c(p, f, r)
  out
}
