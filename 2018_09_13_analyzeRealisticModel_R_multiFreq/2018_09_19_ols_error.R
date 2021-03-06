# c·½ϋaZo
RSS <- numeric(6)
for (i in 1:6) {
  x_tmp = t(rate_EMCLs[i,])
  y_tmp = t(aveSOS2[i,])
  lm_tmp = lm(y_tmp ~ x_tmp)
  RSS[i]<- deviance(lm_tmp)
}

frq2 <- frq/1e3
# c·½ϋaΖόgΜvbg
windows()
par(ps=20,mex=1.5)
plot(RSS ~ frq2[,1],
     xlab = 'όg [kHz]',
     ylab = 'c·½ϋa@@[(m/s)^2]',
     log = 'x',
     type = 'o',
     ylim = c(0 ,30))
savePlot(filename = "//Azlab-fs01/€Ί/Βlwork/|ΰ(Π)/result/2018_09_13_analyzeRealisticModel_R_multiFreq/RSS.png",type="png")
dev.off()

# c·vbg
windows()
plot(lm_tmp)

# c·Zo
residuals <- matrix(nrow = 6,ncol = 25)
for (i in 1:6) {
  x_tmp = t(rate_EMCLs[i,])
  y_tmp = t(aveSOS2[i,])
  lm_tmp = lm(y_tmp ~ x_tmp)
  residuals[i,]<- residuals(lm_tmp)
}

# c·ΖB_YΜvbg
windows()
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  x_tmp <- t(varDiff_Point[i,])
  y_tmp <- residuals[i,]
  plot(y_tmp ~ x_tmp,
       xlab = 'oHΟ»ͺΜͺU',
       ylab = 'c·[m/s]',
       main = sprintf('Sόg=%d kHz',frq2[i,1]),
       ylim=c(-2.2,2.2))
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
savePlot(filename = "//Azlab-fs01/€Ί/Βlwork/|ΰ(Π)/result/2018_09_13_analyzeRealisticModel_R_multiFreq/residuals_vs_displaceFAS.png",type="png")
dev.off()

# c·ΖΌlSΜvbg
windows()
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  x_tmp <- t(fwhm_cell[i,])
  y_tmp <- residuals[i,]
  plot(y_tmp ~ x_tmp,
       xlab = 'ΌlS[fq]',
       ylab = 'c·[m/s]',
       main = sprintf('Sόg=%d kHz',frq2[i,1]),
       ylim=c(-2.2,2.2))
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
savePlot(filename = "//Azlab-fs01/€Ί/Βlwork/|ΰ(Π)/result/2018_09_13_analyzeRealisticModel_R_multiFreq/residuals_vs_FWHM.png",type="png")
dev.off()



# oHΟ»ͺΜͺUΖεtFpΜvbg
windows()
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  x_tmp <- t(refractive_dispersion[i,])
  y_tmp <- t(varDiff_Point[i,])
  plot(y_tmp ~ x_tmp,
       xlab = 'εtFp[rad]',
       ylab = 'oHΟ»ͺΜͺU',
       main = sprintf('Sόg=%d kHz',frq2[i,1])
  )
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
savePlot(filename = "//Azlab-fs01/€Ί/Βlwork/|ΰ(Π)/result/2018_09_13_analyzeRealisticModel_R_multiFreq/var_diff_point_vs_refractive_dispersion.png",type="png")
dev.off()


# ΌlSΖtFpΜvbg
windows()
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  x_tmp <- t(refractive_dispersion[i,])
  y_tmp <- t(fwhm_cell[i,])
  plot(y_tmp ~ x_tmp,
       xlab = 'tFp[rad]',
       ylab = 'ΌlS',
       main = sprintf('Sόg=%d kHz',frq2[i,1])
  )
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
savePlot(filename = "//Azlab-fs01/€Ί/Βlwork/|ΰ(Π)/result/2018_09_13_analyzeRealisticModel_R_multiFreq/FWHM_vs_refractive_dispersion.png",type="png")
dev.off()