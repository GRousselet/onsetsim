# true onset = 160 ms, F=81, max at F=126
true_onset <- 160
Xf <- seq(0, 500, 2)
Nf <- length(Xf)
temp1 <- vector(mode = "numeric", length = Nf)
erp <- dnorm(seq(-1.5,1.5,length.out=47), 0, 1)
erp <- erp - min(erp)
erp <- erp / max(erp)
temp2.half <- c(rep(0, 79), erp, rep(0, 79+46))
