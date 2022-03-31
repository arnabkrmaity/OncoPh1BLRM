Graph2Combo <- function (P = NULL, Pcat = NULL, Pcat.labels = NULL, est = "mean",
          low = "2.5%", high = "97.5%", int.col = NULL, int.crit = 0.5,
          drug1.labels = NULL, drug2.labels = NULL, data = NULL, Pcutoffs = NULL,
          Pcutoffs.lty = 3, est.pch = 16, v.space = 0.2, h.space = 0.1,
          cex.01 = 0.75, cex.lab = 1, cex.main = 1, cex.axis1 = 1,
          cex.axis2 = 1, cex.data = 0.8, back.col = "white", title = "",
          title.cex = 1, cex.Pcat.labels = 0.8, adj.Pcat.labels = 0.2,
          pos.Pcat.labels = c(1, 1), lab.drug2 = "Agent 2", lab.drug1 = "Agent 1")
{
  if (length(dim(P)) != 3)
    stop("P must be of dimension d1 x d2 x #summaries")
  if (length(dim(Pcat)) != 3)
    stop("Pcat must be of dimension d1 x d2 x #intervals")
  d1 = dim(P)[1]
  d2 = dim(P)[3]
  Nint = dim(Pcat)[2]
  est = P[, est, ]
  est.low = P[, low,]
  est.high = P[, high,]
  est = matrix(est, d1, d2)
  est.low = matrix(est.low, d1, d2)
  est.high = matrix(est.high, d1, d2)
  if (is.null(drug1.labels))
    drug1.labels = dimnames(P)[[1]]
  if (is.null(drug2.labels))
    drug2.labels = dimnames(P)[[3]]
  xx = seq(0, by = 10 * (1 + h.space), length = d2)
  yy = seq(0, by = 10 * (1 + v.space), length = d1)
  Nint = dim(Pcat)[2]
  if (length(int.crit) == 1)
    int.crit = rep(int.crit, Nint)
  if (is.null(int.col))
    int.col = lapply(1:Nint, function(e) c("white", "grey"))
  int.col = lapply(int.col, function(e) {
    if (length(e) == 1)
      return(c(e, e))
    if (length(e) == 2)
      return(e)
  })
  xcord = seq(0, 10, length = Nint + 2)
  plot(c(0, (1 + h.space) * d2 * 10), c(0, (1 + v.space) *d1 * 10), axes = FALSE, type = "n", xlab = lab.drug2,
       ylab = lab.drug1, main = title, cex.main = title.cex,
       cex.lab = cex.lab, cex.main = cex.main)
  for (j in 1:d2) {
    for (jj in 1:d1) {
      rect(xx[j], yy[jj], xx[j] + 10, yy[jj] + 10, col = back.col)
      for (j2 in 1:1) {
        lines(rep(xx[j] + xcord[j2 + 1], 2), c(yy[jj], yy[jj] + 10))
      }
      if (!is.null(data))
        text((xx[j] + xx[j] + xcord[2])/2, yy[jj] + 9, data[jj, j], adj = 0.5, cex = cex.data)
      points(xx[j] + xcord[2]/2, yy[jj] + est[jj, j] *10, pch = est.pch)
      if (!is.null(est.low) & !is.null(est.high))
        lines(rep(xx[j] + xcord[2]/2, 2), yy[jj] + 10 *c(est.low[jj, j], est.high[jj, j]), lwd = 3)
      if (!is.null(Pcutoffs)) {
        for (jjj in 1:length(Pcutoffs)) lines(c(xx[j], xx[j] + xcord[2]), rep(yy[jj] + 10 * Pcutoffs[jjj],2), lty = Pcutoffs.lty)
      }
      for (j2 in 1:Nint) {
        rect(xx[j] + xcord[j2 + 1] + 0.2 * (xcord[j2 + 2] - xcord[j2 + 1]), yy[jj], xx[j] + xcord[j2 +1] + 0.8 * (xcord[j2 + 2] - xcord[j2 + 1]),
             yy[jj] + 10*Pcat[jj, j2, j], col = ifelse(Pcat[jj,j2, j] <= int.crit[j2], int.col[[j2]][1], int.col[[j2]][2]))

        if (!is.null(Pcat.labels)) {
          if (pos.Pcat.labels[1] == 0 & pos.Pcat.labels[2] == 0)
            text(xx[j] + xcord[j2 + 1] + 0.2 * (xcord[j2 + 2] - xcord[j2 + 1]), yy[jj] + 9, Pcat.labels[j2],
                 adj = adj.Pcat.labels, cex = cex.Pcat.labels)
          if (pos.Pcat.labels[1] == j & pos.Pcat.labels[2] == jj)
            text(xx[j] + xcord[j2 + 1] + 0.2 * (xcord[j2 + 2] - xcord[j2 + 1]), yy[jj] + 9, Pcat.labels[j2],
                 adj = adj.Pcat.labels, cex = cex.Pcat.labels)
        }
        lines(c(xx[j] + xcord[j2 + 1], xx[j] + xcord[j2 + 2]), rep(yy[jj] + 10 * int.crit[j2], 2), lty = 3)
      }
    }
  }
  text(rep(0, d2 * 2), c(yy, yy + 10), labels = c(rep(0, d1),
                                                  rep(1, d1)), adj = 1, cex = cex.01)
  axis(1, at = xx + 5, drug2.labels, tick = FALSE, cex.axis = cex.axis1)
  axis(2, at = yy + 5, drug1.labels, tick = FALSE, cex.axis = cex.axis2)
  return(NULL)
}
