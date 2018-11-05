# load(file = "data/workspace.RData")
# 単回帰プロット
windows()
r_single_reg_multifrq <- matrix(nrow = 6,ncol = 1)
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  x_tmp <- t(rate_EMCLs[i,])
  y_tmp <- t(aveSOS2[i,])
  plot(y_tmp ~ x_tmp,
       xlab = 'EMCL[%]',
       ylab = '平均音速値[m/s]',
       main = sprintf('中心周波数=%d MHz',frq[i,1]/1e3),
       ylim = c(1525 ,1565))
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
savePlot(filename = "//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/linear_regression_multi_frq.png",type="png")
dev.off()

# load(file = "data/workspace.RData")
# 単回帰プロット
windows()
r_single_reg_multifrq <- matrix(nrow = 6,ncol = 1)
par(mfrow=c(2,3),ps=14,mex=0.7,pty='s')
for (i in 1:6) {
  x_tmp <- t(rate_EMCLs[i,])
  y_tmp <- t(aveSOS2[i,])
  plot(y_tmp ~ x_tmp,
       xlab = 'IMAT[%]',
       ylab = 'Area Average SOS[m/s]',
       main = sprintf('Center frequency =%d kHz',frq[i,1]/1e3),
       ylim = c(1525 ,1565))
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
savePlot(filename = "H:/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/linear_regression_multi_frq_eng.png",type="png")
dev.off()

# 単回帰プロット
# 500 kHzで経路平均音速と平面波送信(単一ch送受信)とのIMAT推定精度比較を行なう．
windows()
r_single_reg_multifrq <- matrix(nrow = 6,ncol = 1)
par(ps=20,mex=1.0,pty='s')
for (i in 3) {
  # each path
  x_tmp <- t(rate_EMCLs[i,])
  y_tmp <- t(aveSOS2[i,])
  plot(y_tmp ~ x_tmp,
       xlab = 'IMAT[%]',
       ylab = 'Area Average SOS[m/s]',
       main = sprintf('Center frequency =%d kHz',frq[i,1]/1e3),
       ylim = c(1530 ,1570),
       col = 'blue')
  abline(lm(y_tmp ~ x_tmp))
  result <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result, interval="confidence", level = 0.95)
  conf.predict <- predict(result, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "blue")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "blue",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "blue",lty=2)
  # single transimit-receive [plane wave]
  par(new=T)
  x_tmp <- t(rate_EMCLs[i,])
  y_tmp <- t(aveSOS_plwv[i,])
  plot(y_tmp ~ x_tmp, ann = F, ylim = c(1530 ,1570),col='darkgreen')
  abline(lm(y_tmp ~ x_tmp))
  result2 <- lm(y_tmp ~ x_tmp)
  conf.interval <- predict(result2, interval="confidence", level = 0.95)
  conf.predict <- predict(result2, interval="prediction", level = 0.95)
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 2], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.interval[order(x_tmp), 3], col = "darkgreen")
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 2], col = "darkgreen",lty=2)
  lines(x_tmp[order(x_tmp)], conf.predict[order(x_tmp), 3], col = "darkgreen",lty=2)
}
legend("bottomleft",legend = c("proposed","single-channel"),col = c("blue","darkgreen"),pch = c(1,1))
savePlot(filename = "H:/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/linear_regression_multi_frq_eng_vs_plwv.png",type="png")
dev.off()

#残差ヒストグラム表示
windows()
residuals_multifrq <- matrix(nrow = 6,ncol = 25)
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  x_tmp = t(rate_EMCLs[i,])
  y_tmp = t(aveSOS2[i,])
  result <- lm(y_tmp ~ x_tmp)
  for (j in 1:25){
    residuals_multifrq[i,j] <- as.vector(result$residuals[j]) 
  }
  hist(residuals_multifrq[i,],breaks ='Sturges',freq = FALSE,
       xlab = '残差[m/s]',
       main = sprintf('中心周波数=%d MHz',frq[i,1]/1e3),
       xlim = c(-2.5, 2.5),
       ylim = c(0,1.0)
  )

}
savePlot(filename = "//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/residuals_multi_frq.png",type="png")
dev.off()

#残差プロット
windows()
par(ps=18,mex=1.2,pty='s')
plot(result,which=1)
savePlot(filename = "//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/residuals_fitted_500kHz.png",type="png")
dev.off()

windows()
par(ps=18,mex=1.2,pty='s')
plot(result,which=2)
savePlot(filename = "//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/normal_q-q_500kHz.png",type="png")
dev.off()

#平均音速値推定誤差
x <- seq(-50,50,0.05)
windows()
par(mfrow=c(2,3),ps=18,mex=1.2,pty='s')
for (i in 1:6) {
  a <- t(v_error[,i])
  hist(a,breaks = 'Scott',freq = FALSE,
       main = sprintf('中心周波数=%d kHz',frq[i,1]/1e3),
       xlab='経路平均音速推定誤差[m/s]',xlim=c(-60,60),ylim=c(0,0.12))
  m <- mean(a)	# 標本平均列の平均
  s <- sd(a)	# 標本平均列の標準偏差
  curve(dnorm(x, mean=m, sd=s), type="l",add=T, col="red")
}
savePlot(filename = "//Azlab-fs01/東研究室/個人work/竹内(ひ)/result/2018_10_09_again_analysis_with_correct_conditions/2018_09_13_analyzeRealisticModel_R_multiFreq/v_estimation_errors2.png",type="png")
dev.off()

#素子ペアごとの平均音速推定値誤差の評価
mean_v_error <- matrix(nrow = 6,ncol = 1)
std_v_error <- matrix(nrow = 6,ncol = 1)
for (i in 1:6) {
  mean_v_error[i,1] <- mean(v_error[,i])
  std_v_error[i,1] <- sqrt(var(v_error[,i]))*(length(v_error[,i])-1)/(length(v_error[,i]))
}
t.test(v_error[,1])
