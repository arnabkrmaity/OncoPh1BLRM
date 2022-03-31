Beta.ms2ab <- function (m, s, probs = c(0.025, 0.5, 0.975)) 
{
  names(m) = NULL
  names(s) = NULL
  n = m * (1 - m)/s/s - 1
  a = n * m
  b = n * (1 - m)
  qntls <- do.call(rbind, mapply(qbeta, shape1 = a, shape2 = b, 
                                 MoreArgs = list(p = probs), SIMPLIFY = FALSE))
  colnames(qntls) = paste(round(100 * probs, 1), "%", sep = "")
  return(list(a = a, b = b, n = n, qntls = qntls))
}